#-------------------------------------------------------------------------------
# This script downloads Global Forest Watch (GFW) layers for the analysis area
# and time period.
#-------------------------------------------------------------------------------

options(timeout = max(60000000, getOption("timeout")))

packages <- list(
  "rjson", "dplyr", "tidyr", "purrr", "terra", "stars", "sf", "readr",
  "geodata", "gdalcubes", "rredlist", "stringr", "tmaptools", "ggplot2", "rstac",
  "lubridate", "RCurl", "codetools"
)
suppressPackageStartupMessages(lapply(packages, library, character.only = TRUE)) # Load libraries - packages

input <- biab_inputs()
print("Inputs: ")
print(input)

if (!dir.exists(outputFolder)) {
  dir.create(outputFolder, recursive = TRUE, showWarnings = FALSE)
}

path_script <- Sys.getenv("SCRIPT_LOCATION")
source(file.path(path_script, "data/filterCubeRangeFunc.R"), echo = TRUE)
source(file.path(path_script, "data/loadCubeFunc.R"), echo = TRUE)

# Write larger cubes directly to disk before opening them with terra. This
# avoids converting a stars proxy in cube_to_raster(), which fails for large
# areas with recent versions of stars/abind.
cube_to_terra <- function(cube, prefix) {
  cube_dir <- file.path(outputFolder, "gdalcubes")
  dir.create(cube_dir, recursive = TRUE, showWarnings = FALSE)
  tif_prefix <- basename(tempfile(pattern = paste0(prefix, "_"), tmpdir = cube_dir))
  tif_paths <- gdalcubes::write_tif(
    cube,
    dir = cube_dir,
    prefix = tif_prefix,
    COG = TRUE,
    creation_options = list(COMPRESS = "DEFLATE")
  )
  if (length(tif_paths) == 0) {
    biab_error_stop(paste0("No GeoTIFF was produced for ", prefix, "."))
  }
  terra::rast(tif_paths)
}

#-------------------------------------------------------------------------------
# Prepare user inputs
#-------------------------------------------------------------------------------
spat_res <- ifelse(is.null(input$spat_res), 500, input$spat_res)
t_0 <- input$t_0
t_n <- input$t_n
t_n_code <- as.numeric(substr(t_n, start = 3, stop = 4))

srs <- paste0(input$crs$CRS$authority, ":", input$crs$CRS$code)
sf_srs <- sf::st_crs(srs)
srs_cube <- srs

area_type <- input$area_type

if (area_type == "Country or region") {
  bbox_values <- input$crs$bbox
  if (is.null(bbox_values) || length(bbox_values) != 4) {
    biab_error_stop("The selected country or region does not contain a valid bounding box.")
  }
  names(bbox_values) <- c("xmin", "ymin", "xmax", "ymax")
  area_bbox <- sf::st_as_sfc(sf::st_bbox(
    bbox_values,
    crs = srs
  ))
} else if (area_type == "Polygon") {
  if (is.null(input$area_file) || !nzchar(input$area_file)) {
    biab_error_stop("Select a polygon area file.")
  }
  area_bbox <- sf::st_read(input$area_file, quiet = TRUE)
} else if (area_type == "Raster") {
  if (is.null(input$area_file) || !nzchar(input$area_file)) {
    biab_error_stop("Select a raster area file.")
  }
  raster_area <- terra::rast(input$area_file)
  raster_crs <- terra::crs(raster_area, proj = TRUE)
  if (is.na(raster_crs) || !nzchar(raster_crs)) {
    biab_error_stop("The area raster must have a coordinate reference system.")
  }
  raster_extent <- terra::ext(raster_area)
  area_bbox <- sf::st_as_sfc(sf::st_bbox(
    c(
      xmin = raster_extent$xmin,
      ymin = raster_extent$ymin,
      xmax = raster_extent$xmax,
      ymax = raster_extent$ymax
    ),
    crs = raster_crs
  ))
} else {
  biab_error_stop("Area type must be Country or region, Polygon, or Raster.")
}

if (is.na(sf::st_crs(area_bbox))) {
  biab_error_stop("The selected download area must have a coordinate reference system.")
}

# The STAC data are downloaded for the bounding box of the selected area.
sf_ext_srs <- area_bbox |>
  sf::st_transform(sf_srs) |>
  sf::st_bbox()

print(sf_ext_srs)

#-------------------------------------------------------------------------------
# Download all GFW layers once
#-------------------------------------------------------------------------------
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
r_GFW_TC <- suppressWarnings(cube_to_terra(cube_GFW_TC, "treecover2000"))
print("========== Base forest layer downloaded ==========")

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
    cube = cube_GFW_loss, min = 1, type_min = 1,
    max = t_0_code, type_max = 1, value = FALSE
  )
  r_loss_before_t0 <- suppressWarnings(
    cube_to_terra(cube_loss_before_t0, "loss_before_t0")
  )
  r_loss_before_t0 <- terra::classify(r_loss_before_t0, rcl = cbind(NA, 0))
  cube_loss_period <- funFilterCube_range(
    cube = cube_GFW_loss, min = t_0_code, type_min = 2,
    max = t_n_code, type_max = 1, value = FALSE
  )
} else {
  cube_loss_period <- funFilterCube_range(
    cube = cube_GFW_loss, min = 1, type_min = 1,
    max = t_n_code, type_max = 1, value = FALSE
  )
  # Not used for a 2000 start, but keep the interface between steps constant.
  r_loss_before_t0 <- NULL
}

r_year_loss <- suppressWarnings(cube_to_terra(cube_loss_period, "loss_period"))
r_year_loss <- terra::classify(r_year_loss, rcl = cbind(NA, 0))
if (is.null(r_loss_before_t0)) {
  r_loss_before_t0 <- terra::init(r_year_loss, 0)
}
print("========== Forest loss layer downloaded and processed ==========")

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
r_GFW_gain <- cube_to_terra(cube_GFW_gain, "gain")
r_GFW_gain <- terra::classify(r_GFW_gain, rcl = cbind(0, NA))
print("========== Forest gain layer downloaded and processed ==========")

#-------------------------------------------------------------------------------
# Save outputs for the forest-loss analysis step
#-------------------------------------------------------------------------------
write_gfw_raster <- function(raster, filename) {
  output_path <- file.path(outputFolder, filename)
  suppressWarnings(terra::writeRaster(
    raster,
    output_path,
    gdal = c("COMPRESS=DEFLATE", "TFW=YES"),
    filetype = "COG",
    overwrite = TRUE
  ))
  output_path
}

tree_cover_path <- write_gfw_raster(r_GFW_TC, "GFW_treecover2000.tif")
loss_before_t0_path <- write_gfw_raster(r_loss_before_t0, "GFW_loss_before_t0.tif")
loss_period_path <- write_gfw_raster(r_year_loss, "GFW_loss_period.tif")
gain_path <- write_gfw_raster(r_GFW_gain, "GFW_gain.tif")

biab_output("tree_cover", tree_cover_path)
biab_output("loss_before_t0", loss_before_t0_path)
biab_output("loss_period", loss_period_path)
biab_output("forest_gain", gain_path)

print("========== GFW download completed ==========")
