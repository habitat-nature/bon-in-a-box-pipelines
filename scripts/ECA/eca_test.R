## Libraries
library("rjson")
library("dplyr")
library("sf")
library("ggplot2")


## Reading inputs
input <- biab_inputs()
print("Inputs: ")
print(input)


## Read AOI
aoi <- sf::st_read(
     input$aoi_polygon,
     quiet = TRUE
)

# Error/ condition checking
if (nrow(aoi) == 0){
     biab_error_stop("AOI has no geometries!")
}


## Save results
aoi_output_path <- file.path(outputFolder,"aoi.gpkg")
sf::st_write(
     aoi,
     aoi_output_path,
     quiet = TRUE
)
biab_output("aoi", aoi_output_path)


## Make a plot
test_plot_aoi <- ggplot2::ggplot() +
  ggplot2::geom_sf(data = aoi, fill = "blue", color = "blue", alpha = 0.5) +
  ggplot2::theme_minimal()


## Save results
plot_output_path <- file.path(outputFolder, "test_plot.png")
ggplot2::ggsave(plot_output_path, test_plot_aoi)
biab_output("test_plot", plot_output_path)


## Info
biab_info("Some information message")
