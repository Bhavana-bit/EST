"""
Palaeoclimate Stability and Modern Biodiversity
Project Pipeline & Data Downloader Script
"""

import os
import sys
import urllib.request
import zipfile
import pandas as pd
import numpy as np

def print_dataset_summary():
    summary = """
================================================================================
 project: Palaeoclimate Stability and Modern Biodiversity
 Dataset Overview & Access Protocols
================================================================================

1. PaleoClim Data (http://www.paleoclim.org/)
   - Purpose: Palaeoclimate reconstructions (LGM ~21ka, Mid-Holocene ~6ka, LIG ~130ka, Pliocene ~3.2Ma).
   - Key Variables: BIO1 (Annual Mean Temp), BIO12 (Annual Precip), Temp & Precip Seasonality.
   - Access: Direct zip downloads or `paleoclim` python/R API.

2. WorldClim Palaeoclimate (https://www.worldclim.org/)
   - Purpose: High-resolution present (2.1) and downscaled past climate (LGM/CMIP5 models).
   - Key Variables: bio_1 to bio_19.
   - Access: HTTP raster downloads (GeoTIFF, 2.5', 5', 10' resolution).

3. GBIF Occurrence Data (https://www.gbif.org/)
   - Purpose: Global georeferenced species observations for calculating observed species richness.
   - Access: `pygbif` Python library or GBIF REST API (`https://api.gbif.org/v1/occurrence/search`).

4. IUCN Red List Spatial Data (https://www.iucnredlist.org/resources/spatial-data-download)
   - Purpose: Expert-validated species geographical range polygons.
   - Access: GeoPackage / Shapefiles for Mammals, Amphibians, Reptiles, Birds (BirdLife).

5. Copernicus DEM (GLO-30 / GLO-90) (https://registry.opendata.aws/copernicus-dem/)
   - Purpose: Global elevation & topographic heterogeneity (slope, roughness) controls.
   - Access: AWS S3 Open Data (`s3://copernicus-dem-30m`) or OpenTopography / EarthEngine API.

================================================================================
"""
    print(summary)

if __name__ == "__main__":
    print_dataset_summary()
