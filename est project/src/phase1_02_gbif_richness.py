"""
phase1_02_gbif_richness.py
--------------------------
Phase 1 Step 2: GBIF Modern Biodiversity Processing.
Strict Enforcement: NO SYNTHETIC FALLBACKS.
Fails loudly if official GBIF download file is missing.
"""

import sys
from pathlib import Path
import numpy as np
import pandas as pd

sys.path.append(str(Path(__file__).resolve().parent.parent))
import config

def process_phase1_gbif():
    print("=" * 80)
    print(f"PHASE 1: GBIF UNIQUE SPECIES RICHNESS PROCESSING")
    print("=" * 80)
    
    gbif_csv = config.GBIF_DIR / "gbif_download_raw.csv"
    gbif_tsv = config.GBIF_DIR / "gbif_download_raw.tsv"
    
    target_file = None
    if gbif_csv.exists() and gbif_csv.stat().st_size > 100:
        target_file = gbif_csv
    elif gbif_tsv.exists() and gbif_tsv.stat().st_size > 100:
        target_file = gbif_tsv
        
    if target_file is None:
        raise FileNotFoundError(
            f"[FATAL ERROR] Official GBIF occurrence download dataset missing in {config.GBIF_DIR}!\n"
            f"Expected file: gbif_download_raw.csv or gbif_download_raw.tsv\n"
            f"Pipeline cannot continue without a real verified GBIF occurrence dataset."
        )
        
    print(f"Reading official GBIF occurrence dataset: {target_file}")
    sep = '\t' if target_file.suffix == '.tsv' else ','
    df = pd.read_csv(target_file, sep=sep, low_memory=False)
    
    rows, cols = config.GRID_SHAPE
    species_richness_grid = np.zeros((rows, cols), dtype=np.int32)
    record_density_grid = np.zeros((rows, cols), dtype=np.int32)
    
    # Require latitude and longitude columns
    if 'decimalLatitude' not in df.columns or 'decimalLongitude' not in df.columns:
        raise ValueError("[FATAL ERROR] GBIF dataset missing decimalLatitude or decimalLongitude columns!")
        
    species_col = 'speciesKey' if 'speciesKey' in df.columns else ('species' if 'species' in df.columns else 'scientificName')
    
    # Filter valid coordinates
    df_valid = df.dropna(subset=['decimalLatitude', 'decimalLongitude', species_col])
    df_valid = df_valid[
        (df_valid['decimalLatitude'] >= -90) & (df_valid['decimalLatitude'] <= 90) &
        (df_valid['decimalLongitude'] >= -180) & (df_valid['decimalLongitude'] <= 180) &
        ~((df_valid['decimalLatitude'] == 0) & (df_valid['decimalLongitude'] == 0))
    ]
    
    row_indices = np.clip(((90.0 - df_valid['decimalLatitude']) / 180.0 * rows).astype(int), 0, rows - 1)
    col_indices = np.clip(((df_valid['decimalLongitude'] + 180.0) / 360.0 * cols).astype(int), 0, cols - 1)
    
    df_valid['row_idx'] = row_indices
    df_valid['col_idx'] = col_indices
    
    cell_species = df_valid.groupby(['row_idx', 'col_idx'])[species_col].nunique()
    cell_records = df_valid.groupby(['row_idx', 'col_idx'])[species_col].count()
    
    for (r, c), u_count in cell_species.items():
        species_richness_grid[r, c] = u_count
        record_density_grid[r, c] = cell_records.get((r, c), 0)
        
    np.save(config.PROCESSED_DATA_DIR / "gbif_species_richness.npy", species_richness_grid)
    np.save(config.PROCESSED_DATA_DIR / "gbif_record_density.npy", record_density_grid)
    
    print(f"[SUCCESS] Processed {len(df_valid)} valid occurrences across {len(cell_species)} grid cells.")

if __name__ == "__main__":
    process_phase1_gbif()
