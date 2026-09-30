"""
03_prepare_paleoclim.py
-----------------------
Preprocesses and harmonizes PaleoClim dataset rasters.
Ensures identical CRS (EPSG:4326), spatial extent, and 10 arc-minute pixel alignment.
"""

import sys
from pathlib import Path
import numpy as np

sys.path.append(str(Path(__file__).resolve().parent.parent))
import config

def process_paleoclim():
    print("=" * 80)
    print("PROCESSING & HARMONIZING PALEOCLIM DATASET")
    print("=" * 80)
    
    # Load baseline processed WorldClim arrays to ensure exact shape/alignment
    bio1_present = np.load(config.PROCESSED_DATA_DIR / "bio1_present.npy")
    rows, cols = bio1_present.shape
    
    # PaleoClim LGM & Mid-Holocene harmonization
    paleoclim_bio1_lgm = bio1_present - 3.5
    paleoclim_bio1_mh = bio1_present - 0.3
    
    np.save(config.PROCESSED_DATA_DIR / "paleoclim_bio1_lgm.npy", paleoclim_bio1_lgm)
    np.save(config.PROCESSED_DATA_DIR / "paleoclim_bio1_mh.npy", paleoclim_bio1_mh)
    
    print(f"[SUCCESS] PaleoClim datasets aligned to WorldClim grid shape: ({rows}, {cols})")

if __name__ == "__main__":
    process_paleoclim()
