# WorldClim Phase 1 Inputs

The local WorldClim metadata and raster headers were inspected directly. Each period has BIO1, and BIO1 values are stored as integer °C × 10 with NoData `-32768`. The pipeline masks NoData, converts BIO1 to °C, and aligns past rasters to the present BIO1 transform using bilinear geospatial reprojection.

| Period | BIO1 path | Format | CRS | Dimensions | Nominal resolution | Bounds |
|---|---|---|---|---|---|---|
| Present baseline (1960–1990) | `data/raw/worldclim/present/bio/bio_1/` | ESRI Grid | EPSG:4326 | 2160 columns × 900 rows | 10 arc-minutes | approximately −180° to 180° longitude, −60° to 90° latitude |
| LGM CCSM4 (~22 ka) | `data/raw/worldclim/lgm/cclgmbi1.tif` | GeoTIFF | EPSG:4326 | 2160 columns × 900 rows | 10 arc-minutes | −180° to 180° longitude, −60° to 90° latitude |
| Mid-Holocene CCSM4 (~6 ka) | `data/raw/worldclim/mid_holocene/ccmidbi1.tif` | GeoTIFF | EPSG:4326 | 2160 columns × 900 rows | 10 arc-minutes | −180° to 180° longitude, −60° to 90° latitude |

Present BIO1 has transform resolution approximately 0.166666675359° and origin (−180°, 90.000007823°). LGM and mid-Holocene BIO1 have resolution 1/6° and origin (−180°, 90°). Although dimensions and nominal coverage match, their transforms differ slightly; the present raster is the authoritative target grid.

The corresponding BIO12 files are present in all three datasets (`bio_12/w001001.adf`, `cclgmbi12.tif`, `ccmidbi12.tif`) and encode annual precipitation in mm/year. Phase 1 currently analyzes BIO1 only.