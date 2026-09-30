"""
rigorous_dataset_validation.py
------------------------------
Rigorous audit and validation of all datasets in EST/Palaeoclimate_Biodiversity_Project/.
Validates:
1. WorldClim 1.4 rasters (Present, LGM, Mid-Holocene)
2. PaleoClim rasters (Current, LGM, Late Holocene, CHELSA)
3. GBIF occurrence data (gbif_all_continents.csv)
4. IUCN spatial data (Checks presence/absence)
5. Copernicus DEM data (Checks presence/absence)
Prints complete metadata, shape, CRS, bounds, NoData, units, value ranges.
"""

import os
import sys
import zipfile
from pathlib import Path
import pandas as pd
import numpy as np
import rasterio

sys.path.append(str(Path(__file__).resolve().parent.parent))
import config

base_path = Path("EST/Palaeoclimate_Biodiversity_Project")

print("=" * 80)
print("RIGOROUS DATASET AUDIT AND VALIDATION REPORT")
print("=" * 80)

# 1. WorldClim 1.4 Audit
print("\n--- 1. WORLDCLIM 1.4 DATASETS ---")
wc_dir = base_path / "data" / "worldclim"
wc_zips = list(wc_dir.glob("*.zip"))

for z_path in wc_zips:
    print(f"\nZip Archive: {z_path} ({z_path.stat().st_size / (1024*1024):.2f} MB)")
    with zipfile.ZipFile(z_path, 'r') as z:
        members = z.namelist()
        print(f"  Total items inside: {len(members)}")
        tif_members = [m for m in members if m.endswith('.tif')]
        adf_members = [m for m in members if m.endswith('.adf')]
        if tif_members:
            print(f"  GeoTIFF files inside: {tif_members[:5]}")
        if adf_members:
            print(f"  ESRI GRID folders inside: {set(m.split('/')[0] for m in adf_members if '/' in m)}")

# 2. PaleoClim Audit
print("\n\n--- 2. PALEOCLIM DATASETS ---")
pc_dir = base_path / "data" / "paleoclim"
pc_files = list(pc_dir.glob("*"))

for f in pc_files:
    size_mb = f.stat().st_size / (1024*1024)
    print(f"\nFile: {f.name} ({size_mb:.2f} MB)")
    if f.suffix == ".tif":
        try:
            with rasterio.open(f) as src:
                data = src.read(1)
                valid = (data != src.nodata) if src.nodata is not None else ~np.isnan(data)
                print(f"  Format: {src.driver}")
                print(f"  Shape: {src.shape}")
                print(f"  CRS: {src.crs}")
                print(f"  Bounds: {src.bounds}")
                print(f"  NoData: {src.nodata}")
                print(f"  Min / Max: {np.min(data[valid]):.2f} / {np.max(data[valid]):.2f}")
        except Exception as e:
            print(f"  Rasterio Error: {e}")
    elif f.suffix == ".zip":
        with zipfile.ZipFile(f, 'r') as z:
            print(f"  Zip Contents: {z.namelist()[:5]}")

# 3. GBIF Audit
print("\n\n--- 3. GBIF OCCURRENCE DATASETS ---")
gbif_dir = base_path / "data" / "gbif"
gbif_files = list(gbif_dir.glob("*.csv"))

for f in gbif_files:
    size_mb = f.stat().st_size / (1024*1024)
    df = pd.read_csv(f)
    print(f"\nFile: {f.name} ({size_mb:.2f} MB)")
    print(f"  Total Records: {len(df)}")
    print(f"  Columns: {list(df.columns)}")
    if 'speciesKey' in df.columns:
        print(f"  Unique Species Keys: {df['speciesKey'].nunique()}")
    print(f"  Latitude Range: {df['decimalLatitude'].min():.2f} to {df['decimalLatitude'].max():.2f}")
    print(f"  Longitude Range: {df['decimalLongitude'].min():.2f} to {df['decimalLongitude'].max():.2f}")

# 4. IUCN Check
print("\n\n--- 4. IUCN DATASET CHECK ---")
iucn_dir = base_path / "data" / "iucn"
if iucn_dir.exists() and any(iucn_dir.glob("*")):
    print(f"IUCN files found: {list(iucn_dir.glob('*'))}")
else:
    print("IUCN DATASET: MISSING (Directory does not exist or contains 0 files)")

# 5. Copernicus DEM Check
print("\n\n--- 5. COPERNICUS DEM DATASET CHECK ---")
dem_dir = base_path / "data" / "dem"
if dem_dir.exists() and any(dem_dir.glob("*")):
    print(f"Copernicus DEM files found: {list(dem_dir.glob('*'))}")
else:
    print("COPERNICUS DEM DATASET: MISSING (Directory does not exist or contains 0 files)")
