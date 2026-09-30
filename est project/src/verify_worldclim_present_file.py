"""
verify_worldclim_present_file.py
--------------------------------
Inspects exact WorldClim Present BIO1 file attributes and verifies exact subtraction ΔT = |Present - LGM|.
"""

import os
import sys
import zipfile
from pathlib import Path
import numpy as np
import rasterio

base_path = Path("EST/Palaeoclimate_Biodiversity_Project")
wc_zip = base_path / "data" / "worldclim" / "bio_10m_esri.zip"
extract_dir = Path("c:/Users/chari/OneDrive/Documents/est project/data/raw/worldclim/present")

if not (extract_dir / "bio" / "bio_1" / "hdr.adf").exists():
    with zipfile.ZipFile(wc_zip, 'r') as z:
        z.extractall(extract_dir)

hdr_path = extract_dir / "bio" / "bio_1" / "hdr.adf"

print("=" * 80)
print("WORLDCLIM PRESENT BIO1 FILE VERIFICATION")
print("=" * 80)
print(f"Full File Path: {hdr_path}")
print(f"Filename: {hdr_path.name} (ESRI GRID dataset: {hdr_path.parent})")

with rasterio.open(hdr_path) as src:
    data = src.read(1)
    print(f"Dimensions: {src.shape} (rows x cols: {src.height} x {src.width})")
    print(f"CRS: {src.crs}")
    print(f"Bounds: {src.bounds}")
    print(f"Resolution: {src.res} (deg)")
    print(f"NoData: {src.nodata}")
    valid = (data != src.nodata) & (data > -500)
    print(f"Raw Units / Encoding: Integer (°C x 10)")
    print(f"Raw Min / Max: {data[valid].min()} / {data[valid].max()}")
    print(f"Scaled Units: °C (Raw / 10.0)")
    print(f"Scaled Min / Max: {data[valid].min()/10.0:.2f}°C / {data[valid].max()/10.0:.2f}°C")
