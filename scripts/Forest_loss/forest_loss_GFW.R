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
#path_script <- "scripts"
#outputFolder <- "/home/jurietheron/Projects/bon-in-a-box-pipelines/output/Forest_loss/"
if (!dir.exists(file.path(outputFolder))) {
  dir.create(outputFolder, recursive = TRUE, showWarnings = FALSE)
} else {
  print("dir exists")
}
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

# AOI for the analysis to run within
#v_path_bbox_analysis <- "/home/jurietheron/Projects/bon-in-a-box-pipelines/scripts/Forest_loss/pnmb.gpkg"
v_path_bbox_analysis <- if (is.null(input$sf_bbox)) {
  NA
} else {
  input$sf_bbox
}

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
#min_forest_Forest_Guinea <- 30
#min_forest_Maritime_Guinea <- 25
#min_forest_Middle_Guinea <- 20
#min_forest_Upper_Guinea <- 35
min_forest <- c(
  "Forest Guinea" = min_forest_Forest_Guinea,
  "Maritime Guinea" = min_forest_Maritime_Guinea,
  "Middle Guinea" = min_forest_Middle_Guinea,
  "Upper Guinea" = min_forest_Upper_Guinea
)

# Max forest threshold for GFW (level of forest for the species)
max_forest_Forest_Guinea <- if (is.null(input$max_forest_Forest_Guinea)) {
  NA
} else {
  input$max_forest_Forest_Guinea
}
max_forest_Maritime_Guinea <- if (is.null(input$max_forest_Maritime_Guinea)) {
  NA
} else {
  input$max_forest_Maritime_Guinea
}
max_forest_Middle_Guinea <- if (is.null(input$max_forest_Middle_Guinea)) {
  NA
} else {
  input$max_forest_Middle_Guinea
}
max_forest_Upper_Guinea <- if (is.null(input$max_forest_Upper_Guinea)) {
  NA
} else {
  input$max_forest_Upper_Guinea
}
#max_forest_Forest_Guinea <- 100
#max_forest_Maritime_Guinea <- 100
#max_forest_Middle_Guinea <- 100
#max_forest_Upper_Guinea <- 100
max_forest <- c(
  "Forest Guinea" = max_forest_Forest_Guinea,
  "Maritime Guinea" = max_forest_Maritime_Guinea,
  "Middle Guinea"   = max_forest_Middle_Guinea,
  "Upper Guinea"    = max_forest_Upper_Guinea
)

# define time range
#t_0 <- 2010
#t_n <- 2020
t_0 <- input$t_0
t_n <- input$t_n # should be larger than t_0

#-------------------------------------------------------------------------------
# Conditional checks
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------------------------------------------
# 1. Load inputs
#-------------------------------------------------------------------------------------------------------------------
# Load AOI and create BBOX
sf_bbox_analysis <- st_read(v_path_bbox_analysis) |> 
  st_transform(sf_srs)
sf_ext_srs <- sf_bbox_analysis |> st_bbox()
print(sf_ext_srs)

# Load Guinea forest groupings polygon
gin_for_shape <- st_read("/home/jurietheron/Projects/bon-in-a-box-pipelines/scripts/Forest_loss/gin_admbnda_adm1_ocha.gpkg") |> 
  st_transform(sf_srs)
gin_for_shape <- st_read(file.path(path_script, "Forest_loss/gin_admbnda_adm1_ocha.gpkg")) |> 
  st_transform(sf_srs)
print(gin_for_shape)

# See how many gin shapes intersect with the AOI
gin_for_shape_intersect <- gin_for_shape[
  unlist(
    sf::st_intersects(
      sf_bbox_analysis,
      gin_for_shape
    )
  ), 
]

# Perform analysis per Guinea forest grouping
habitat_change_map_path <- c()
for(gin_group_index in 1:length(gin_for_shape_intersect)){
  # Select geometry
  shape <- gin_for_shape_intersect[gin_group_index,]
  print(paste0("========== Processing: ",shape$group, " =========="))

  # Subset the tree threshold values
  max_threshold <- max_forest[shape$group]
  min_threshold <- min_forest[shape$group]

  #-------------------------------------------------------------------------------------------------------------------
  # 2. GFW base forest layer processing
  #-------------------------------------------------------------------------------------------------------------------

  print("========== Downloadeding and processing the base forest layer ==========")

  # Load forest base map from 2000
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

  # Apply thresholding
  cube_GFW_TC_threshold <<- funFilterCube_range(
    cube_GFW_TC,
    min = min_threshold,
    max = max_threshold,
    value = FALSE
  )

  # Load into memory / download
  r_GFW_TC_threshold <- suppressWarnings(
    cube_to_raster(
      cube_GFW_TC_threshold,
      format = "terra"
    )
  )

  # Reclassify
  r_GFW_TC_threshold <- terra::classify(
    r_GFW_TC_threshold,
    rcl = cbind(NA, 0)
  )

  print("========== Base forest layer downloaded and proccessed ==========")

  #-------------------------------------------------------------------------------------------------------------------
  # 3. GFW forest loss processing
  #-------------------------------------------------------------------------------------------------------------------

  print("========== Downloadeding and processing forest loss maps ==========")

  # Load forest loss maps
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

  # Convert years to 2-digit lossyear codes
  t_n_code <- as.numeric(substr(t_n, start = 3, stop = 4))

  # If t_0 != 2000, rebase the forest layer by removing loss that occurred before t_0
  if (t_0 != 2000) {
    t_0_code <- as.numeric(substr(t_0, start = 3, stop = 4))
    cube_loss_before_t0 <- funFilterCube_range(
      cube = cube_GFW_loss, min = 1, type_min = 1, max = t_0_code, type_max = 1, value = FALSE
    )
    r_loss_before_t0 <- suppressWarnings(cube_to_raster(cube_loss_before_t0, format = "terra"))
    r_loss_before_t0 <- terra::classify(r_loss_before_t0, rcl = cbind(NA, 0))
    r_GFW_TC_threshold <- terra::classify(r_GFW_TC_threshold - r_loss_before_t0, rcl = cbind(-1, 0))
    print("Rebased forest layer to t_0")

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

  # Convert loss cube to raster
  r_year_loss <- suppressWarnings(cube_to_raster(cube_loss_period, format = "terra"))
  r_year_loss <- terra::classify(r_year_loss, rcl = cbind(NA, 0))

  # Mask to AOI and then to forested area
  r_GFW_TC_threshold_mask <- terra::mask(r_GFW_TC_threshold, sf_bbox_analysis)
  r_year_loss_masked <- terra::mask(r_year_loss, r_GFW_TC_threshold_mask, maskvalues = 1, inverse = TRUE)

  # Reclassify: keep only loss pixels
  r_year_loss_mask <- terra::classify(r_year_loss_masked, rcl = cbind(0, NA))

  print("========== Forest loss layer downloaded and proccessed ==========")

  #-------------------------------------------------------------------------------------------------------------------
  # 4. GFW forest gain processing
  #-------------------------------------------------------------------------------------------------------------------

  print("========== Downloadeding and processing forest gain maps ==========")

  # Load forest gain maps
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

  # Turn cube to raster
  r_GFW_gain <- cube_to_raster(
    cube_GFW_gain,
    format = "terra"
  )

  # Reclassify and mask
  r_GFW_gain_mask <- terra::classify(
    terra::mask(
      r_GFW_gain,
      sf_bbox_analysis
    ),
    rcl = cbind(0, NA)
  )

  print("========== Forest gain layer downloaded and proccessed ==========")

  #-------------------------------------------------------------------------------------------------------------------
  # 5. Postprocess GFW data
  #-------------------------------------------------------------------------------------------------------------------

  # Post process GFW
  r_GFW_TC_threshold_mask[r_GFW_TC_threshold_mask > 0] <- 1 # turn no change to 1
  r_year_loss_mask[r_year_loss_mask > 0] <- 2 # turn loss value to 2
  r_GFW_gain_mask[r_GFW_gain_mask > 0] <- 3 # turn gain value to 3

  # Put no change, loss, and gain together in one raster
  v1 <- merge(r_year_loss_mask, r_GFW_TC_threshold_mask) # merge loss and no change
  v2 <- merge(r_GFW_gain_mask, v1) # merge gain
  v3 <- terra::classify(v2, rcl = cbind(0, NA)) # turn 0 to NA

  # Recategorize
  # Add colours
  # Add labels


  # Save
  habitat_change_map_path[gin_group_index] <- file.path(outputFolder, shape$group, paste0(shape$group, "_GFW_loss.tiff"))
  dir.create(habitat_change_map_path[gin_group_index], recursive = TRUE, showWarnings = FALSE)
  habitat_change_map <- suppressWarnings(
    terra::writeRaster(
      v3,
      habitat_change_map_path[gin_group_index],
      gdal = c("COMPRESS=DEFLATE", "TFW=YES"),
      filetype = "COG",
      overwrite = TRUE
    )
  )

  print("========== Map of forest loss generated ==========")
}

print("========== Outputting results ==========")





# Outputting result
biab_output("habitat_change_map", habitat_change_map_path)

print("========== Analysis completed ==========")