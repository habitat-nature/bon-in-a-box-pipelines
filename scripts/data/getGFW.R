#-------------------------------------------------------------------------------
# This script dowloads Global Forest Watch (GFW) layers for the analysis area and time period.
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Initial set up
#-------------------------------------------------------------------------------
# Set time out for large downloads
options(timeout = max(60000000, getOption("timeout")))

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
