_Author(s): Maria Isabel Arce-Plata, Guillaume Larocque, Jurie Theron, and Samara Manzin_

## Introduction

The Measure Forest Change pipeline uses Global Forest Watch (GFW) data to map forest loss, forest gain, and baseline forest with no detected change. Run it for the whole selected area with one minimum tree-cover threshold, or supply subregion boundaries and set a separate threshold for each subregion.

## Uses

Use the pipeline to locate forest loss, compare the share of baseline forest lost across subregions, and summarize forest change in hectares for monitoring, conservation planning, or reporting.

## Pipeline limitations

- GFW tree-cover loss records canopy removal. It can reflect permanent deforestation, harvesting, fire, or other disturbances.
- The selected years determine the loss analysis. The forest-gain layer covers the fixed 2000–2012 source period.
- The source data resolution and classification accuracy limit how reliably the pipeline detects small or fragmented changes. Choosing a smaller output cell size does not improve source detail.
- Results cover the selected area within available GFW coverage. When you supply subregions, results cover their overlap with that area.

## Before you start

Choose a projected coordinate reference system (CRS) that uses metres. Prepare a georeferenced area file if you plan to use a polygon or raster. For subregion reporting, prepare a GeoPackage with one polygon layer and an attribute field containing the names. Use non-overlapping subregions. The pipeline combines features with the same name.

## Running the pipeline

### Pipeline inputs

| Input | What to enter |
| --- | --- |
| Area type | Choose **Country or region** to use the selected bounding box, **Polygon** to clip to an uploaded boundary, or **Raster** to use the uploaded raster's rectangular extent. |
| Country, region, bounding box, and output CRS | Select a projected CRS that uses metres. In **Country or region** mode, select a country or region or draw a bounding box. |
| Area file | Supply a polygon or georeferenced raster for the corresponding area type. Leave this blank for **Country or region**. |
| Output cell size (metres) | Set the side length of each square output cell, for example 50 metres. |
| Start year and End year | Choose a start year of 2000 or later and an end year after the start year, up to 2020. Forest loss is measured after the start year through the end year. |
| Subregion boundaries (optional) | Choose the subregion GeoPackage, or leave this blank to analyze the selected area as a whole. |
| Subregion name field | When using subregions, enter the exact attribute field name from the GeoPackage, such as `group`. |
| Minimum tree-cover thresholds | When using subregions, paste the header `region, threshold` followed by one row per subregion. Match the polygon names and capitalization, and enter a percentage from 0 to 100. If the field accepts only one line, separate rows with a backslash (`\`). |
| Whole-area minimum tree cover (%) (optional) | When Subregion boundaries is blank, enter a percentage from 0 to 100 for the whole area. The default is 15%. |

For the Guinea configuration, choose **Country or region**, select Guinea, and select the bundled `gin_admbnda_adm1_ocha.gpkg` as **Subregion boundaries**. Enter `group` as the **Subregion name field** and paste the following into **Minimum tree-cover thresholds**. [See what the subregions look like in the Guinea example](https://github.com/user-attachments/assets/e6ccac41-c645-4de0-88a1-ef24a94a6c8a).

```text
region, threshold
Forest Guinea, 25
Maritime Guinea, 15
Middle Guinea, 15
Upper Guinea, 15
```

### Pipeline steps

#### 1. Download Global Forest Watch data

The pipeline retrieves 2000 tree-cover percentages, annual forest loss, and the forest-gain layer for the selected area. It prepares separate loss rasters for loss through the start year and loss after the start year through the end year. [GEO BON STAC catalog](https://stac.geobon.org/).

#### 2. Measure forest change

The pipeline applies the minimum tree-cover threshold to the 2000 tree-cover raster. For a start year after 2000, it removes forest lost through that year to establish baseline forest. It combines the baseline with the loss and gain rasters, clips the results to the analysis boundary, and summarizes them across the whole area and by subregion when supplied.

### Pipeline outputs

| Output | How to read it |
| --- | --- |
| Forest change map (raster) | GeoTIFF with values `1` for no detected change in baseline forest, `2` for forest loss, and `3` for forest gain. |
| Forest change area (ha) | Bar chart showing hectares in each change category. |
| Forest loss (% of baseline forest) | Bar chart showing loss area divided by the sum of loss area and unchanged baseline forest area. Gain is excluded from the denominator. A missing bar can indicate that the area had no baseline forest. |
| GFW tree cover in 2000 | Source raster with tree-cover percentages for 2000. |
| GFW forest gain | Source raster showing detected gain during 2000–2012. |
| GFW loss during analysis period | Source raster showing detected canopy loss after the start year through the end year. |

The **Analysis area total** in the charts summarizes the merged map. With subregions, the charts also show results for each analyzed subregion.

## Example

## Troubleshooting

- **The selected boundary layer does not contain the requested subregion name field.** Check the available fields in the error message. For the Guinea example, select `gin_admbnda_adm1_ocha.gpkg` and enter `group`.
- **Missing tree-cover thresholds for a subregion.** Check the names and capitalization against the GeoPackage and provide exactly one row per analyzed subregion.
- **The threshold table is rejected.** Paste the table directly into the text field. Use the header `region, threshold` and enter numbers from 0 to 100 without percent signs.
- **The CRS is in latitude and longitude or does not use metres.** Select a projected CRS with metre units.
- **The analysis area does not overlap any supplied subregion.** Check that the selected country, bounding box, polygon, or raster extent overlaps the supplied subregion polygons.

## References

- Hansen, M. C., et al. (2013). High-Resolution Global Maps of 21st-Century Forest Cover Change. _Science, 342_(6160), 850–853. [DOI](https://doi.org/10.1126/science.1244693).
- Hansen/UMD/Google/USGS/NASA. Global Forest Change data. [Global Forest Watch](https://www.globalforestwatch.org/map/).
