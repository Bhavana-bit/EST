"""
audit_workspace_and_system.py
-----------------------------
Audits existing files in the workspace and user directories to locate WorldClim 1.4,
PaleoClim, GBIF, IUCN, and DEM rasters/archives.
Inspects raster metadata (dimensions, CRS, transform, NoData, dtype, value ranges) using rasterio if found.
"""

import os
from pathlib import Path
import rasterio

USER_PROFILE = Path(os.environ.get("USERPROFILE", "C:/Users/chari"))
WORKSPACE = Path("c:/Users/chari/OneDrive/Documents/est project")

search_paths = [
    WORKSPACE,
    USER_PROFILE / "Downloads",
    USER_PROFILE / "Desktop",
    USER_PROFILE / "OneDrive" / "Documents",
    USER_PROFILE / "Documents"
]

raster_extensions = {".tif", ".bil", ".hdr", ".asc", ".nc", ".adf"}
zip_extensions = {".zip", ".tar", ".gz"}

found_rasters = []
found_zips = []

print("=" * 80)
print("AUDITING WORKSPACE AND LOCAL SYSTEM FOR DATASETS")
print("=" * 80)

for base in search_paths:
    if not base.exists():
        continue
    print(f"Scanning directory: {base}")
    for root, dirs, files in os.walk(base):
        root_path = Path(root)
        # Skip hidden or system folders to avoid infinite loops
        if ".git" in root_path.parts or ".antigravity" in root_path.parts or "AppData" in root_path.parts:
            continue
        
        # Check for ESRI GRID folders (directories containing w001001.adf)
        if (root_path / "w001001.adf").exists():
            found_rasters.append({
                "type": "ESRI GRID Directory",
                "path": root_path,
                "file": "w001001.adf"
            })
            
        for f in files:
            ext = Path(f).suffix.lower()
            if ext in raster_extensions and f != "w001001.adf":
                found_rasters.append({
                    "type": "Raster File",
                    "path": root_path / f,
                    "file": f
                })
            elif ext in zip_extensions:
                if any(k in f.lower() for k in ["wc", "bio", "ccsm", "lgm", "holocene", "paleo", "gbif", "iucn", "dem", "10m"]):
                    found_zips.append({
                        "type": "Archive",
                        "path": root_path / f,
                        "file": f
                    })

print(f"\n--- FINDINGS SUMMARY ---")
print(f"Found {len(found_rasters)} potential raster files/directories.")
print(f"Found {len(found_zips)} potential dataset archives.")

if found_rasters:
    print("\nINSpectING FOUND RASTERS:")
    for item in found_rasters[:15]:
        p = item["path"]
        print(f"\nPath: {p}")
        try:
            with rasterio.open(p) as src:
                print(f"  Shape: {src.shape} (rows x cols)")
                print(f"  CRS: {src.crs}")
                print(f"  Bounds: {src.bounds}")
                print(f"  NoData: {src.nodatavals}")
                print(f"  Data Types: {src.dtypes}")
                # Read sample window to check scaling
                sample = src.read(1, window=((0, 100), (0, 100)))
                print(f"  Sample Min/Max: {sample.min()} / {sample.max()}")
        except Exception as e:
            print(f"  Rasterio Inspect Error: {e}")

if found_zips:
    print("\nFOUND DATASET ARCHIVES:")
    for z in found_zips:
        print(f"  - {z['path']} ({z['path'].stat().st_size / (1024*1024):.2f} MB)")
