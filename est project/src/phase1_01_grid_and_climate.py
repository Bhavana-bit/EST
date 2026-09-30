"""
phase1_01_grid_and_climate.py
-----------------------------
Phase 1 Step 1: Common Grid & WorldClim Climate Processing.
Strict Enforcement: NO SYNTHETIC FALLBACKS.
Fails loudly if official WorldClim raster files are missing.
"""

import sys
from pathlib import Path
import numpy as np
import rasterio

sys.path.append(str(Path(__file__).resolve().parent.parent))
import config

def process_phase1_climate():
    print("=" * 80)
    print("PHASE 1: WORLDCLIM CLIMATE PROCESSING & ΔT STABILITY METRIC")
    print("=" * 80)
    
    present_path = config.WORLDCLIM_DIR / "present" / "bio1.bil"
    lgm_path = config.WORLDCLIM_DIR / "lgm" / "bio1.bil"
    
    # Check for alternative GeoTIFF extension if BIL is absent
    if not present_path.exists():
        present_path = config.WORLDCLIM_DIR / "present" / "bio1.tif"
    if not lgm_path.exists():
        lgm_path = config.WORLDCLIM_DIR / "lgm" / "bio1.tif"
        
    if not present_path.exists() or not lgm_path.exists():
        raise FileNotFoundError(
            f"[FATAL ERROR] Official WorldClim rasters missing!\n"
            f"Expected Present BIO1 at: {present_path}\n"
            f"Expected LGM BIO1 at:     {lgm_path}\n"
            f"Pipeline cannot continue without real verified raster datasets."
        )
        
    print(f"Opening Present BIO1: {present_path}")
    with rasterio.open(present_path) as src_pres:
        bio1_present_raw = src_pres.read(1)
        nodata_pres = src_pres.nodata
        meta_pres = src_pres.meta
        
    print(f"Opening LGM BIO1: {lgm_path}")
    with rasterio.open(lgm_path) as src_lgm:
        bio1_lgm_raw = src_lgm.read(1)
        nodata_lgm = src_lgm.nodata
        
    # Temperature scaling correction (°C x 10 -> °C)
    bio1_present_degc = np.where(bio1_present_raw == nodata_pres, np.nan, bio1_present_raw / 10.0)
    bio1_lgm_degc = np.where(bio1_lgm_raw == nodata_lgm, np.nan, bio1_lgm_raw / 10.0)
    
    # Compute absolute temperature difference ΔT (°C)
    delta_t_lgm = np.abs(bio1_present_degc - bio1_lgm_degc)
    
    # Save outputs
    np.save(config.PROCESSED_DATA_DIR / "bio1_present.npy", bio1_present_degc)
    np.save(config.PROCESSED_DATA_DIR / "bio1_lgm.npy", bio1_lgm_degc)
    np.save(config.PROCESSED_DATA_DIR / "delta_t_lgm.npy", delta_t_lgm)
    
    print(f"[SUCCESS] Processed real WorldClim rasters successfully.")
    print(f"  Shape: {meta_pres['height']} rows x {meta_pres['width']} cols")
    print(f"  CRS: {meta_pres['crs']}")

if __name__ == "__main__":
    process_phase1_climate()
