#-------------------------------------------------------------------------------
# This script measures forest loss across Guinea using GFW data
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
#path_script <- "scripts"
path_script <- Sys.getenv("SCRIPT_LOCATION")
source(file.path(path_script, "data/filterCubeRangeFunc.R"), echo = TRUE)
source(file.path(path_script, "data/loadCubeFunc.R"), echo = TRUE)

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

# AOI for the analysis to run within (optional)
#v_path_bbox_analysis <- "/home/jurietheron/Projects/bon-in-a-box-pipelines/scripts/Forest_loss/pnmb.gpkg"
v_path_bbox_analysis <- if (is.null(input$sf_bbox)) NA else input$sf_bbox

# Min forest threshold for GFW (level of forest for the species)
min_forest_Forest_Guinea <- if (is.null(input$min_forest_Forest_Guinea)) {
  NA
} else {
  input$min_forest_Forest_Guinea
}
min_forest_Maritime_Guinea <- if (is.null(input$min_forest_Maritime_Guinea)) {
  NA
} else {
  input$min_forest_Maritime_Guinea
}
min_forest_Middle_Guinea <- if (is.null(input$min_forest_Middle_Guinea)) {
  NA
} else {
  input$min_forest_Middle_Guinea
}
min_forest_Upper_Guinea <- if (is.null(input$min_forest_Upper_Guinea)) {
  NA
} else {
  input$min_forest_Upper_Guinea
}
#min_forest_Forest_Guinea <- 25
#min_forest_Maritime_Guinea <- 28
#min_forest_Middle_Guinea <- 23
#min_forest_Upper_Guinea <- 24
min_forest <- c(
  "Forest Guinea" = min_forest_Forest_Guinea,
  "Maritime Guinea" = min_forest_Maritime_Guinea,
  "Middle Guinea" = min_forest_Middle_Guinea,
  "Upper Guinea" = min_forest_Upper_Guinea
)

# define time range
#t_0 <- 2010
#t_n <- 2020
t_0 <- input$t_0
t_n <- input$t_n # should be larger than t_0

# Convert years to 2-digit lossyear codes
t_n_code <- as.numeric(substr(t_n, start = 3, stop = 4))

#-------------------------------------------------------------------------------
# Conditional checks
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------------------------------------------
# 1. Load inputs
#-------------------------------------------------------------------------------------------------------------------
# Load Guinea forest groupings polygon
#gin_for_shape <- st_read("/home/jurietheron/Projects/bon-in-a-box-pipelines/scripts/Forest_loss/gin_admbnda_adm1_ocha.gpkg") |>
#  st_transform(sf_srs)
gin_for_shape <- st_read(file.path(path_script, "Forest_loss/gin_admbnda_adm1_ocha.gpkg")) |>
  st_transform(sf_srs)
print(gin_for_shape)

# Load optional AOI — if provided, clip Guinea regions and downloads to it; otherwise use full Guinea extent
if (!is.na(v_path_bbox_analysis)) {
  sf_bbox_analysis <- st_read(v_path_bbox_analysis) |> st_transform(sf_srs)
  sf_ext_srs <- sf_bbox_analysis |> st_bbox()
  gin_for_shape_intersect <- gin_for_shape[
    unlist(sf::st_intersects(sf_bbox_analysis, gin_for_shape)),
  ]
} else {
  sf_bbox_analysis <- NULL
  sf_ext_srs <- gin_for_shape |> st_bbox()
  gin_for_shape_intersect <- gin_for_shape
}
print(sf_ext_srs)

#-------------------------------------------------------------------------------------------------------------------
# 2. Download all GFW layers once (loop-invariant)
#-------------------------------------------------------------------------------------------------------------------

# Download raw treecover2000 cube — thresholding is per-region and applied inside the loop
print("========== Downloading base forest layer ==========")
cube_GFW_TC <- load_cube(
  stac_path = "https://stac.geobon.org/",
  limit = 1000,
  collections = c("gfw-treecover2000"),
  bbox = sf_ext_srs,
  spatial.res = spat_res,
  srs.cube = srs_cube,
  temporal.res = "P1Y",
  t0 = "2000-01-01",
  t1 = "2000-12-31",
  resampling = "bilinear"
)
print("========== Base forest layer downloaded ==========")

# Download forest loss cube and derive the period loss raster
print("========== Downloading and processing forest loss maps ==========")
cube_GFW_loss <- load_cube(
  stac_path = "https://stac.geobon.org/",
  limit = 1000,
  collections = c("gfw-lossyear"),
  bbox = sf_ext_srs,
  srs.cube = srs_cube,
  spatial.res = spat_res,
  temporal.res = "P1Y",
  t0 = "2000-01-01",
  t1 = "2000-12-31",
  resampling = "mode",
  aggregation = "first"
)

if (t_0 != 2000) {
  t_0_code <- as.numeric(substr(t_0, start = 3, stop = 4))
  cube_loss_before_t0 <- funFilterCube_range(
    cube = cube_GFW_loss, min = 1, type_min = 1, max = t_0_code, type_max = 1, value = FALSE
  )
  r_loss_before_t0 <- suppressWarnings(cube_to_raster(cube_loss_before_t0, format = "terra"))
  r_loss_before_t0 <- terra::classify(r_loss_before_t0, rcl = cbind(NA, 0))
  # Loss between t_0 and t_n (exclusive of t_0, inclusive of t_n)
  cube_loss_period <- funFilterCube_range(
    cube = cube_GFW_loss, min = t_0_code, type_min = 2, max = t_n_code, type_max = 1, value = FALSE
  )
} else {
  # Loss from 2001 to t_n
  cube_loss_period <- funFilterCube_range(
    cube = cube_GFW_loss, min = 1, type_min = 1, max = t_n_code, type_max = 1, value = FALSE
  )
}

r_year_loss <- suppressWarnings(cube_to_raster(cube_loss_period, format = "terra"))
r_year_loss <- terra::classify(r_year_loss, rcl = cbind(NA, 0))
print("========== Forest loss layer downloaded and processed ==========")

# Download forest gain cube and rasterize
print("========== Downloading and processing forest gain maps ==========")
cube_GFW_gain <- load_cube(
  stac_path = "https://io.biodiversite-quebec.ca/stac",
  limit = 1000,
  collections = c("gfw-gain"),
  bbox = sf_ext_srs,
  srs.cube = srs_cube,
  spatial.res = spat_res,
  temporal.res = "P1Y",
  t0 = "2000-01-01",
  t1 = "2000-12-31",
  resampling = "near"
)
r_GFW_gain <- cube_to_raster(cube_GFW_gain, format = "terra")
r_GFW_gain_aoi <- if (!is.null(sf_bbox_analysis)) {
  terra::classify(terra::mask(r_GFW_gain, sf_bbox_analysis), rcl = cbind(0, NA))
} else {
  terra::classify(r_GFW_gain, rcl = cbind(0, NA))
}
print("========== Forest gain layer downloaded and processed ==========")

#-------------------------------------------------------------------------------------------------------------------
# 3. Perform analysis per Guinea forest grouping
#-------------------------------------------------------------------------------------------------------------------
#gin_group_index = 1
habitat_change_map <- list()
for (gin_group_index in 1:length(gin_for_shape_intersect)) {
  # Select geometry, clipped to the AOI if provided
  shape <- gin_for_shape_intersect[gin_group_index, ]
  if (!is.null(sf_bbox_analysis)) {
    shape <- st_intersection(shape, sf_bbox_analysis)
  }
  print(paste0("========== Processing: ", shape$group, " =========="))

  # Subset the tree threshold values
  max_threshold <- 100
  min_threshold <- min_forest[shape$group]

  #-------------------------------------------------------------------------------------------------------------------
  # 3.1 Apply per-region threshold to base forest layer
  #-------------------------------------------------------------------------------------------------------------------
  print("========== Processing base forest layer ==========")

  cube_GFW_TC_threshold <<- funFilterCube_range(
    cube_GFW_TC,
    min = min_threshold,
    max = max_threshold,
    value = FALSE
  )
  r_GFW_TC_threshold <- suppressWarnings(cube_to_raster(cube_GFW_TC_threshold, format = "terra"))
  r_GFW_TC_threshold <- terra::classify(r_GFW_TC_threshold, rcl = cbind(NA, 0))

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
  habitat_change_map[[gin_group_index]] <- v3_clip

  print("========== Layer of forest loss generated ==========")
}

print("========== Outputting results ==========")

# Merge into single forest loss layer
rsrc <- terra::sprc(habitat_change_map)
habiatat_change_map_merge <- terra::merge(rsrc)

# Save
habitat_change_map_path <- file.path(outputFolder,"GFW_loss.tiff")
dir.create(dirname(habitat_change_map_path), recursive = TRUE, showWarnings = FALSE)
suppressWarnings(
  terra::writeRaster(
    habiatat_change_map_merge,
    habitat_change_map_path,
    gdal = c("COMPRESS=DEFLATE", "TFW=YES"),
    filetype = "COG",
    overwrite = TRUE
  )
)

# Outputting result
biab_output("habitat_change_map", habitat_change_map_path)

#-------------------------------------------------------------------------------------------------------------------
# 4. Barplot of % per change category, per Guinea region
#-------------------------------------------------------------------------------------------------------------------
print("========== Building forest change barplot ==========")

# Category codes from the per-region rasters: 1 = no change, 2 = loss, 3 = gain
# Forest loss % = loss / baseline forest, where baseline forest = no change + loss (gain excluded)
loss_pct_df <- purrr::map_dfr(seq_along(habitat_change_map), function(i) {
  r <- habitat_change_map[[i]]
  freq_tab <- terra::freq(r) # columns: layer, value, count
  counts <- setNames(freq_tab$count, as.character(freq_tab$value))
  n_nochange <- ifelse(is.na(counts["1"]), 0, counts["1"])
  n_loss     <- ifelse(is.na(counts["2"]), 0, counts["2"])
  baseline   <- n_nochange + n_loss
  data.frame(
    region   = gin_for_shape_intersect$group[i],
    loss_pct = if (baseline > 0) 100 * n_loss / baseline else NA_real_,
    stringsAsFactors = FALSE
  )
})

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
ggplot2::ggsave(change_barplot_path, p_change, width = 8, height = 5, dpi = 150)
biab_output("forest_change_barplot", change_barplot_path)

#-------------------------------------------------------------------------------------------------------------------
# 5. Barplot of area (ha) per change category, per Guinea region
#-------------------------------------------------------------------------------------------------------------------
print("========== Building forest area barplot ==========")

cat_labels <- c("1" = "No change", "2" = "Loss", "3" = "Gain")
cell_area_ha <- (spat_res^2) / 10000 # projected units assumed metres -> hectares

# Tally area (ha) per category for each stored region raster
change_area_df <- purrr::map_dfr(seq_along(habitat_change_map), function(i) {
  r <- habitat_change_map[[i]]
  freq_tab <- terra::freq(r) # columns: layer, value, count
  data.frame(
    region  = gin_for_shape_intersect$group[i],
    code    = as.character(freq_tab$value),
    area_ha = freq_tab$count * cell_area_ha,
    stringsAsFactors = FALSE
  )
}) |>
  dplyr::filter(!is.na(code), code %in% names(cat_labels)) |>
  dplyr::mutate(category = factor(cat_labels[code], levels = unname(cat_labels)))

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
ggplot2::ggsave(area_barplot_path, p_area, width = 8, height = 5, dpi = 150)
biab_output("forest_area_barplot", area_barplot_path)

print("========== Analysis completed ==========")
