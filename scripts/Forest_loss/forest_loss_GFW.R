#-------------------------------------------------------------------------------
# This script measures forest loss by subregion using GFW data
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Initial set up
#-------------------------------------------------------------------------------
# Set time out for large downloads
options(timeout = max(60000000, getOption("timeout")))

# Load all libraries
packages <- list(
  "rjson", "dplyr", "tidyr", "purrr", "terra", "stars", "sf", "readr",
  "geodata", "gdalcubes", "rredlist", "stringr", "tmaptools", "ggplot2", "rstac",
  "lubridate", "RCurl", "codetools"
)
suppressPackageStartupMessages(lapply(packages, library, character.only = TRUE)) # Load libraries - packages

# Grab user inputs
input <- biab_inputs()
print("Inputs: ")
print(input)

# Load additional helper functions
#outputFolder <- "/home/jurietheron/Projects/bon-in-a-box-pipelines/output/Forest_loss"
if (!dir.exists(file.path(outputFolder))) {
  dir.create(outputFolder, recursive = TRUE, showWarnings = FALSE)
} else {
  print("dir exists")
}

#-------------------------------------------------------------------------------
# Prepare user inputs for analysis
#-------------------------------------------------------------------------------
# spatial resolution
#spat_res <- 50
spat_res <- ifelse(is.null(input$spat_res), 500, input$spat_res)

# Define SRS
#srs <- paste0("EPSG", ":", "31529")
srs <- paste0(input$crs$CRS$authority, ":", input$crs$CRS$code)
check_srs <- grepl("^[[:digit:]]+$", srs)
sf_srs <- if (check_srs) st_crs(as.numeric(srs)) else st_crs(srs)
srs_cube <- suppressWarnings(if (check_srs) {
  authorities <- c("EPSG", "ESRI", "IAU2000", "SR-ORG")
  auth_srid <- paste(authorities, srs, sep = ":")
  auth_srid_test <- map_lgl(auth_srid, ~ !"try-error" %in% class(suppressWarnings(try(st_crs(.x), silent = TRUE))))
  if (sum(auth_srid_test) != 1) print("--- Please specify authority name or provide description of the SRS ---") else auth_srid[auth_srid_test]
} else {
  srs
})

use_subregions <- !is.null(input$region_boundaries) && nzchar(trimws(input$region_boundaries))

# define time range
#t_0 <- 2010
#t_n <- 2020
t_0 <- input$t_0
t_n <- input$t_n # should be larger than t_0

# Load the GFW rasters produced by the separate download step
r_GFW_TC <- terra::rast(input$tree_cover)
r_loss_before_t0 <- terra::rast(input$loss_before_t0)
r_year_loss <- terra::rast(input$loss_period)
r_GFW_gain_aoi <- terra::rast(input$forest_gain)

#-------------------------------------------------------------------------------
# Conditional checks
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------------------------------------------
# 1. Load inputs
#-------------------------------------------------------------------------------------------------------------------
# Limit the analysis to the downloaded footprint and the selected area.
download_extent <- terra::ext(r_GFW_TC)
analysis_coverage <- sf::st_as_sfc(sf::st_bbox(
  c(
    xmin = unname(download_extent$xmin),
    ymin = unname(download_extent$ymin),
    xmax = unname(download_extent$xmax),
    ymax = unname(download_extent$ymax)
  ),
  crs = sf::st_crs(terra::crs(r_GFW_TC))
)) |>
  sf::st_transform(sf_srs)

if (identical(input$area_type, "Polygon")) {
  if (is.null(input$area_file) || !nzchar(input$area_file)) {
    biab_error_stop("Select a polygon area file.")
  }
  area_polygon <- sf::st_read(input$area_file, quiet = TRUE) |>
    sf::st_transform(sf_srs) |>
    sf::st_geometry() |>
    sf::st_union()
  analysis_coverage <- sf::st_intersection(analysis_coverage, area_polygon)
} else if (!use_subregions && identical(input$area_type, "Raster")) {
  if (is.null(input$area_file) || !nzchar(input$area_file)) {
    biab_error_stop("Select a raster area file.")
  }
  area_raster <- terra::rast(input$area_file)
  raster_crs <- sf::st_crs(terra::crs(area_raster))
  if (is.na(raster_crs)) {
    biab_error_stop("The area raster must have a coordinate reference system.")
  }
  raster_footprint <- terra::as.polygons(terra::ext(area_raster), crs = terra::crs(area_raster)) |>
    sf::st_as_sf() |>
    sf::st_geometry() |>
    sf::st_transform(sf_srs)
  analysis_coverage <- sf::st_intersection(analysis_coverage, raster_footprint)
} else if (!use_subregions && identical(input$area_type, "Country or region") &&
           length(input$crs$bbox) == 4L) {
  bbox_values <- setNames(input$crs$bbox, c("xmin", "ymin", "xmax", "ymax"))
  selected_bbox <- sf::st_as_sfc(sf::st_bbox(bbox_values, crs = sf_srs))
  analysis_coverage <- sf::st_intersection(analysis_coverage, selected_bbox)
}

if (length(analysis_coverage) == 0 || !any(as.numeric(sf::st_area(analysis_coverage)) > 0)) {
  biab_error_stop("The selected area does not overlap the downloaded GFW coverage.")
}
if (use_subregions) {
  threshold_text <- trimws(input$region_thresholds)
  if (length(threshold_text) == 0 || !nzchar(threshold_text)) {
    biab_error_stop("Enter the threshold table with the header region, threshold and one row per subregion.")
  }
  # Accept backslashes as row separators when the form uses a single-line text field.
  if (!grepl("\n", threshold_text, fixed = TRUE)) {
    threshold_text <- gsub("\\", "\n", threshold_text, fixed = TRUE)
  } else {
    threshold_text <- gsub("\\\\\\r?\\n", "\n", threshold_text)
  }
  thresholds <- readr::read_csv(
    I(threshold_text),
    col_types = readr::cols(.default = readr::col_character()),
    show_col_types = FALSE
  )
  if (!identical(names(thresholds), c("region", "threshold")) ||
      nrow(thresholds) == 0 || nrow(readr::problems(thresholds)) > 0) {
    biab_error_stop("Enter the header region, threshold followed by one name and percentage per row.")
  }
  thresholds$region <- trimws(thresholds$region)
  if (anyNA(thresholds$region) || any(!nzchar(thresholds$region))) {
    biab_error_stop("Every threshold row must have a subregion name.")
  }
  duplicate_regions <- unique(thresholds$region[duplicated(thresholds$region)])
  if (length(duplicate_regions) > 0) {
    biab_error_stop(paste("Duplicate threshold rows for:", paste(duplicate_regions, collapse = ", ")))
  }
  thresholds$threshold <- suppressWarnings(as.numeric(thresholds$threshold))
  invalid_thresholds <- !is.finite(thresholds$threshold) |
    thresholds$threshold < 0 | thresholds$threshold > 100
  if (any(invalid_thresholds)) {
    biab_error_stop(paste(
      "Minimum tree cover must be a number from 0 to 100 for:",
      paste(thresholds$region[invalid_thresholds], collapse = ", ")
    ))
  }

  regions <- sf::st_read(input$region_boundaries, quiet = TRUE)
  region_fields <- names(sf::st_drop_geometry(regions))
  if (is.null(input$region_field) || !input$region_field %in% region_fields) {
    biab_error_stop(paste(
      "The selected boundary layer does not contain the requested subregion name field.",
      "Available fields:", paste(region_fields, collapse = ", ")
    ))
  }
  if (nrow(regions) == 0 || is.na(sf::st_crs(regions)) ||
      !all(sf::st_geometry_type(regions) %in% c("POLYGON", "MULTIPOLYGON"))) {
    biab_error_stop("The subregion layer must contain polygons and have a defined CRS.")
  }
  region_names <- trimws(as.character(regions[[input$region_field]]))
  if (anyNA(region_names) || any(!nzchar(region_names))) {
    biab_error_stop("Every subregion polygon must have a name.")
  }
  regions <- sf::st_sf(region = region_names, geometry = sf::st_geometry(regions)) |>
    sf::st_transform(sf_srs) |>
    sf::st_make_valid() |>
    sf::st_collection_extract("POLYGON")

  regions <- sf::st_intersection(regions, analysis_coverage)
  # Exclude empty intersections and regions that only touch the footprint.
  regions <- regions[
    !sf::st_is_empty(regions) & as.numeric(sf::st_area(regions)) > 0,
  ]
  if (nrow(regions) == 0) {
    biab_error_stop("The analysis area does not overlap any supplied subregion within the downloaded GFW coverage.")
  }
  regions <- regions |>
    sf::st_collection_extract("POLYGON") |>
    dplyr::group_by(region) |>
    dplyr::summarise(.groups = "drop")
  threshold_rows <- match(regions$region, thresholds$region)
  if (anyNA(threshold_rows)) {
    biab_error_stop(paste(
      "Missing tree-cover thresholds for:",
      paste(regions$region[is.na(threshold_rows)], collapse = ", ")
    ))
  }
  regions$min_tree_cover_pct <- thresholds$threshold[threshold_rows]
} else {
  whole_area_threshold <- input$whole_area_threshold
  if (is.null(whole_area_threshold)) {
    whole_area_threshold <- 15
  }
  if (length(whole_area_threshold) != 1L || !is.numeric(whole_area_threshold) ||
      !is.finite(whole_area_threshold) || whole_area_threshold < 0 || whole_area_threshold > 100) {
    biab_error_stop("Enter a whole-area minimum tree-cover percentage from 0 to 100.")
  }
  regions <- sf::st_sf(
    region = "Analysis area total",
    min_tree_cover_pct = whole_area_threshold,
    geometry = sf::st_union(analysis_coverage)
  )
}

#-------------------------------------------------------------------------------------------------------------------
# 3. Perform analysis per subregion
#-------------------------------------------------------------------------------------------------------------------
habitat_change_map <- vector("list", nrow(regions))
for (region_index in seq_len(nrow(regions))) {
  # Select the region already clipped to the shared analysis area
  shape <- regions[region_index, ]
  print(paste0("========== Processing: ", shape$region, " =========="))

  # Subset the tree threshold values
  max_threshold <- 100
  min_threshold <- shape$min_tree_cover_pct

  #-------------------------------------------------------------------------------------------------------------------
  # 3.1 Apply per-region threshold to base forest layer
  #-------------------------------------------------------------------------------------------------------------------
  print("========== Processing base forest layer ==========")

  r_GFW_TC_threshold <- terra::clamp(
    r_GFW_TC,
    lower = min_threshold,
    upper = max_threshold,
    values = FALSE
  )
  r_GFW_TC_threshold <- terra::classify(r_GFW_TC_threshold, rcl = cbind(NA, 0))

  r_GFW_TC_threshold <- terra::ifel(r_GFW_TC_threshold > 0, 1, 0)

  # Rebase forest layer to t_0 by removing pre-t_0 loss (uses r_loss_before_t0 computed once above)
  if (t_0 != 2000) {
    r_GFW_TC_threshold <- terra::classify(r_GFW_TC_threshold - r_loss_before_t0, rcl = cbind(-1, 0))
    print("Rebased forest layer to t_0")
  }

  print("========== Base forest layer processed ==========")

  #-------------------------------------------------------------------------------------------------------------------
  # 3.2 Mask loss and TC layers to this region's shape
  #-------------------------------------------------------------------------------------------------------------------
  r_GFW_TC_threshold_mask <- terra::mask(r_GFW_TC_threshold, shape)
  r_year_loss_masked <- terra::mask(r_year_loss, r_GFW_TC_threshold_mask, maskvalues = 1, inverse = TRUE)
  r_year_loss_mask <- terra::classify(r_year_loss_masked, rcl = cbind(0, NA))

  # Clip gain layer to this region's shape
  r_GFW_gain_shape <- terra::crop(r_GFW_gain_aoi, shape, mask = TRUE)

  #-------------------------------------------------------------------------------------------------------------------
  # 3.3 Postprocess: reclassify and merge layers
  #-------------------------------------------------------------------------------------------------------------------
  r_GFW_TC_threshold_mask[r_GFW_TC_threshold_mask > 0] <- 1 # no change = 1
  r_year_loss_mask[r_year_loss_mask > 0] <- 2               # loss = 2
  r_GFW_gain_shape[r_GFW_gain_shape > 0] <- 3               # gain = 3

  # Merge loss, no-change, and gain into one raster and clip to shape
  v1 <- merge(r_year_loss_mask, r_GFW_TC_threshold_mask)
  v2 <- merge(r_GFW_gain_shape, v1)
  v3_clip <- terra::crop(
    terra::classify(v2, rcl = cbind(0, NA)),
    shape,
    mask = TRUE
  )

  # Store in list
  habitat_change_map[[region_index]] <- v3_clip

  print("========== Layer of forest loss generated ==========")
}

print("========== Outputting results ==========")

# Merge into single forest loss layer
habitat_change_map_merge <- if (length(habitat_change_map) == 1L) {
  habitat_change_map[[1L]]
} else {
  terra::merge(terra::sprc(habitat_change_map))
}

# Save
habitat_change_map_path <- file.path(outputFolder,"GFW_loss.tiff")
dir.create(dirname(habitat_change_map_path), recursive = TRUE, showWarnings = FALSE)
suppressWarnings(
  terra::writeRaster(
    habitat_change_map_merge,
    habitat_change_map_path,
    gdal = c("COMPRESS=DEFLATE", "TFW=YES"),
    filetype = "COG",
    overwrite = TRUE
  )
)

# Outputting result
biab_output("habitat_change_map", habitat_change_map_path)

#-------------------------------------------------------------------------------------------------------------------
# 4. Barplot of forest loss percentage by region and for the analysis area
#-------------------------------------------------------------------------------------------------------------------
print("========== Building forest change barplot ==========")

# Count the merged raster for the total so shared boundary cells are counted once.
plot_rasters <- habitat_change_map
plot_regions <- regions$region
if (use_subregions) {
  plot_rasters <- c(plot_rasters, list(habitat_change_map_merge))
  plot_regions <- make.unique(c(plot_regions, "Analysis area total"))
}
total_label <- tail(plot_regions, 1L)
region_levels <- plot_regions
plot_width <- max(8, 1.25 * length(plot_regions))
category_codes <- c("1", "2", "3")

change_counts <- purrr::map_dfr(seq_along(plot_rasters), function(i) {
  freq_tab <- terra::freq(plot_rasters[[i]])
  code_indices <- match(category_codes, freq_tab$value, nomatch = 0L)
  data.frame(
    region = factor(plot_regions[i], levels = region_levels),
    code = category_codes,
    count = c(0, freq_tab$count)[code_indices + 1L]
  )
})

# Forest loss % = loss / baseline forest, where baseline forest = no change + loss (gain excluded)
loss_pct_df <- change_counts |>
  dplyr::group_by(region) |>
  dplyr::summarise(
    baseline = sum(count[code %in% c("1", "2")]),
    loss = sum(count[code == "2"]),
    .groups = "drop"
  ) |>
  dplyr::mutate(
    loss_pct = dplyr::if_else(baseline > 0, 100 * loss / baseline, NA_real_)
  )

print(loss_pct_df)

# Barplot: one loss % bar per region
p_change <- ggplot2::ggplot(
  loss_pct_df,
  ggplot2::aes(x = region, y = loss_pct)
) +
  ggplot2::geom_col(fill = "#d95f02", width = 0.7) +
  ggplot2::geom_text(
    ggplot2::aes(label = sprintf("%.1f%%", loss_pct)),
    vjust = -0.3, size = 3
  ) +
  ggplot2::labs(
    title = paste0("Forest loss by region (", t_0, "-", t_n, ")"),
    subtitle = "Loss as % of baseline forest (no change + loss, gain excluded)",
    x = NULL, y = "Forest loss (%)"
  ) +
  ggplot2::theme_minimal() +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 30, hjust = 1))

# Save and register output
change_barplot_path <- file.path(outputFolder, "forest_change_barplot.png")
ggplot2::ggsave(change_barplot_path, p_change, width = plot_width, height = 5, dpi = 150, limitsize = FALSE)
biab_output("forest_change_barplot", change_barplot_path)

#-------------------------------------------------------------------------------------------------------------------
# 5. Barplot of area (ha) per change category by region and for the analysis area
#-------------------------------------------------------------------------------------------------------------------
print("========== Building forest area barplot ==========")

cat_labels <- c("1" = "No change", "2" = "Loss", "3" = "Gain")
cell_area_ha <- (spat_res^2) / 10000 # projected units assumed metres -> hectares

change_area_df <- change_counts |>
  dplyr::filter(count > 0 | region == total_label) |>
  dplyr::mutate(
    area_ha = count * cell_area_ha,
    category = factor(cat_labels[code], levels = unname(cat_labels))
  )

print(change_area_df)

# Grouped barplot: one cluster of bars per region, coloured by category
p_area <- ggplot2::ggplot(
  change_area_df,
  ggplot2::aes(x = region, y = area_ha, fill = category)
) +
  ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.8), width = 0.7) +
  ggplot2::geom_text(
    ggplot2::aes(label = formatC(round(area_ha), format = "d", big.mark = ",")),
    position = ggplot2::position_dodge(width = 0.8),
    vjust = -0.3, size = 2.5
  ) +
  ggplot2::scale_fill_manual(
    values = c("No change" = "#1b9e77", "Loss" = "#d95f02", "Gain" = "#7570b3")
  ) +
  ggplot2::labs(
    title = paste0("Forest change area by region (", t_0, "-", t_n, ")"),
    x = NULL, y = "Area (ha)", fill = NULL
  ) +
  ggplot2::theme_minimal() +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 30, hjust = 1))

# Save and register output
area_barplot_path <- file.path(outputFolder, "forest_area_barplot.png")
ggplot2::ggsave(area_barplot_path, p_area, width = plot_width, height = 5, dpi = 150, limitsize = FALSE)
biab_output("forest_area_barplot", area_barplot_path)

print("========== Analysis completed ==========")
