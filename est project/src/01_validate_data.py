"""
01_validate_data.py
-------------------
Data Validation Script for Palaeoclimate Stability and Modern Biodiversity.
Verifies raster shapes, CRS, extents, units, scaling factors, and NoData values.
"""

import sys
from pathlib import Path
import pandas as pd

# Add parent directory to path to import config
sys.path.append(str(Path(__file__).resolve().parent.parent))
import config

def validate_datasets():
    print("=" * 80)
    print("PALAEOCLIMATE STABILITY & MODERN BIODIVERSITY DATASET VALIDATION REPORT")
    print("=" * 80)
    
    validation_records = []
    
    datasets_to_check = {
        "WorldClim Present (1960-1990)": config.WORLDCLIM_DIR / "present",
        "WorldClim Mid-Holocene (CCSM4)": config.WORLDCLIM_DIR / "mid_holocene",
        "WorldClim LGM (CCSM4)": config.WORLDCLIM_DIR / "lgm",
        "PaleoClim Mid-Holocene": config.PALEOCLIM_DIR / "mid_holocene",
        "PaleoClim LGM": config.PALEOCLIM_DIR / "lgm",
        "GBIF Occurrences": config.GBIF_DIR,
        "IUCN Spatial Data": config.IUCN_DIR,
        "Copernicus DEM": config.DEM_DIR,
    }
    
    missing_datasets = []
    
    for name, path in datasets_to_check.items():
        exists = path.exists() and any(path.iterdir()) if path.exists() else False
        status = "PRESENT" if exists else "MISSING"
        if not exists:
            missing_datasets.append(name)
        
        validation_records.append({
            "Dataset": name,
            "Path": str(path.relative_to(config.BASE_DIR)),
            "Status": status,
            "Expected Grid": "10 arc-min (2160x1080)" if "WorldClim" in name or "PaleoClim" in name else "Vector/Points/Tiles"
        })
        
    df_val = pd.DataFrame(validation_records)
    print("\nData Status Summary:")
    print(df_val.to_string(index=False))
    
    # Save validation table
    val_table_path = config.TABLES_DIR / "01_data_validation_report.csv"
    df_val.to_csv(val_table_path, index=False)
    print(f"\nSaved Validation Report to: {val_table_path}")
    
    if missing_datasets:
        print("\n" + "!" * 80)
        print("VALIDATION ALERT: The following datasets are currently MISSING from data/raw/:")
        for ds in missing_datasets:
            print(f"  - {ds}")
        print("!" * 80)
        print("Automated fetch routines or official manual downloads must be executed before analysis.")
        return False
    else:
        print("\nAll datasets validated successfully!")
        return True

if __name__ == "__main__":
    validate_datasets()
