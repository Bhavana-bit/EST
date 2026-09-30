"""
process_and_compute_results.py
------------------------------
Processes actual datasets found in EST/Palaeoclimate_Biodiversity_Project/:
Handles raster extent alignment between WorldClim 1.4 (900x2160) and PaleoClim.
Calculates real empirical temperature anomalies (DeltaT), GBIF unique species richness (S),
Spearman and Pearson correlations, and outputs statistics, maps, and tables.
"""

import os
import sys
import zipfile
from pathlib import Path
import numpy as np
import pandas as pd
import rasterio
from rasterio.enums import Resampling
from scipy import stats
import matplotlib.pyplot as plt

sys.path.append(str(Path(__file__).resolve().parent.parent))
import config

base_path = Path("EST/Palaeoclimate_Biodiversity_Project")

def run_analysis():
    print("=" * 80)
    print("RUNNING ANALYSIS ON EMPIRICAL REPOSITORY DATASETS")
    print("=" * 80)
    
    wc_dir = base_path / "data" / "worldclim"
    lgm_zip = wc_dir / "cclgmbi_10m.zip"
    mh_zip = wc_dir / "ccmidbi_10m.zip"
    
    lgm_extract = config.RAW_DATA_DIR / "worldclim" / "lgm"
    mh_extract = config.RAW_DATA_DIR / "worldclim" / "mid_holocene"
    
    if not (lgm_extract / "cclgmbi1.tif").exists():
        with zipfile.ZipFile(lgm_zip, 'r') as z:
            z.extractall(lgm_extract)
    if not (mh_extract / "ccmidbi1.tif").exists():
        with zipfile.ZipFile(mh_zip, 'r') as z:
            z.extractall(mh_extract)
            
    # Read WorldClim LGM BIO1 raster (cclgmbi1.tif)
    lgm_tif = lgm_extract / "cclgmbi1.tif"
    print(f"Reading LGM BIO1: {lgm_tif}")
    with rasterio.open(lgm_tif) as src_lgm:
        lgm_raw = src_lgm.read(1)
        meta_lgm = src_lgm.meta
        nodata_lgm = src_lgm.nodata

    # Read PaleoClim LGM & Current rasters
    pc_dir = base_path / "data" / "paleoclim"
    pc_curr_tif = pc_dir / "paleoclim_current_BIO1.tif"
    pc_lgm_tif = pc_dir / "paleoclim_LGM_BIO1.tif"
    
    print(f"Reading PaleoClim Current: {pc_curr_tif}")
    with rasterio.open(pc_curr_tif) as src_pc_cur:
        pc_cur_raw = src_pc_cur.read(
            1,
            out_shape=(meta_lgm['height'], meta_lgm['width']),
            resampling=Resampling.bilinear
        )
        
    print(f"Reading PaleoClim LGM: {pc_lgm_tif}")
    with rasterio.open(pc_lgm_tif) as src_pc_lgm:
        pc_lgm_raw = src_pc_lgm.read(
            1,
            out_shape=(meta_lgm['height'], meta_lgm['width']),
            resampling=Resampling.bilinear
        )
        
    # Scale temperatures (°C x 10 -> °C if needed)
    lgm_degc = np.where((lgm_raw == nodata_lgm) | (lgm_raw < -500), np.nan, lgm_raw / 10.0)
    
    pc_cur_degc = np.where((pc_cur_raw < -500) | (pc_cur_raw > 1000), np.nan, pc_cur_raw)
    if np.nanmax(pc_cur_degc) > 100:
        pc_cur_degc = pc_cur_degc / 10.0
        
    pc_lgm_degc = np.where((pc_lgm_raw < -500) | (pc_lgm_raw > 1000), np.nan, pc_lgm_raw)
    if np.nanmax(pc_lgm_degc) > 100:
        pc_lgm_degc = pc_lgm_degc / 10.0
        
    # Temperature difference DeltaT = |Present - LGM|
    delta_t_worldclim = np.abs(pc_cur_degc - lgm_degc)
    delta_t_paleoclim = np.abs(pc_cur_degc - pc_lgm_degc)
    
    # Process GBIF Occurrence Data
    gbif_csv = base_path / "data" / "gbif" / "gbif_all_continents.csv"
    print(f"\nProcessing GBIF dataset: {gbif_csv}")
    df_gbif = pd.read_csv(gbif_csv)
    print(f"Total GBIF occurrences loaded: {len(df_gbif)}")
    
    rows, cols = meta_lgm['height'], meta_lgm['width']
    species_richness_grid = np.zeros((rows, cols), dtype=np.int32)
    record_count_grid = np.zeros((rows, cols), dtype=np.int32)
    
    df_valid = df_gbif.dropna(subset=['decimalLatitude', 'decimalLongitude', 'speciesKey']).copy()
    
    # WorldClim LGM extent: 90N to 60S (rows = 900, cols = 2160)
    row_idx = np.clip(((90.0 - df_valid['decimalLatitude']) / 150.0 * rows).astype(int), 0, rows - 1)
    col_idx = np.clip(((df_valid['decimalLongitude'] + 180.0) / 360.0 * cols).astype(int), 0, cols - 1)
    
    df_valid.loc[:, 'row_idx'] = row_idx
    df_valid.loc[:, 'col_idx'] = col_idx
    
    cell_species = df_valid.groupby(['row_idx', 'col_idx'])['speciesKey'].nunique()
    cell_records = df_valid.groupby(['row_idx', 'col_idx'])['speciesKey'].count()
    
    for (r, c), count in cell_species.items():
        species_richness_grid[r, c] = count
        record_count_grid[r, c] = cell_records.get((r, c), 0)
        
    print(f"Mapped occurrences across {len(cell_species)} active grid cells.")
    print(f"Max unique species per cell: {species_richness_grid.max()}")
    
    # Join Climate DeltaT & Species Richness per cell
    valid_mask = (species_richness_grid > 0) & (~np.isnan(delta_t_worldclim))
    
    dt_vals = delta_t_worldclim[valid_mask]
    dt_pc_vals = delta_t_paleoclim[valid_mask]
    s_vals = species_richness_grid[valid_mask]
    log_s_vals = np.log1p(s_vals)
    n_cells = len(dt_vals)
    
    # Spearman & Pearson correlations
    rho_wc, p_rho_wc = stats.spearmanr(dt_vals, s_vals)
    r_wc, p_r_wc = stats.pearsonr(dt_vals, log_s_vals)
    
    rho_pc, p_rho_pc = stats.spearmanr(dt_pc_vals, s_vals)
    r_pc, p_r_pc = stats.pearsonr(dt_pc_vals, log_s_vals)
    
    print("\n" + "=" * 80)
    print("COMPUTED EMPIRICAL RESULTS")
    print("=" * 80)
    print(f"Sample Occupied Grid Cells (N): {n_cells}")
    print(f"WorldClim DeltaT vs. Species Richness S:")
    print(f"  Spearman rho: {rho_wc:.4f} (p = {p_rho_wc:.4e})")
    print(f"  Pearson r (log S+1): {r_wc:.4f} (p = {p_r_wc:.4e})")
    print(f"\nPaleoClim DeltaT vs. Species Richness S:")
    print(f"  Spearman rho: {rho_pc:.4f} (p = {p_rho_pc:.4e})")
    print(f"  Pearson r (log S+1): {r_pc:.4f} (p = {p_r_pc:.4e})")
    
    # Save statistical table
    res_df = pd.DataFrame([
        {
            'Climate Source': 'WorldClim 1.4 (Present vs LGM CCSM4)',
            'Sample Size (Cells N)': n_cells,
            'Spearman rho': round(rho_wc, 4),
            'Spearman p-value': f"{p_rho_wc:.4e}",
            'Pearson r (log S+1)': round(r_wc, 4),
            'Pearson p-value': f"{p_r_wc:.4e}"
        },
        {
            'Climate Source': 'PaleoClim v1.0 (Current vs LGM)',
            'Sample Size (Cells N)': n_cells,
            'Spearman rho': round(rho_pc, 4),
            'Spearman p-value': f"{p_rho_pc:.4e}",
            'Pearson r (log S+1)': round(r_pc, 4),
            'Pearson p-value': f"{p_r_pc:.4e}"
        }
    ])
    res_df.to_csv(config.TABLES_DIR / "empirical_statistical_results.csv", index=False)
    
    # Render Plots
    # MAP 1: DeltaT WorldClim
    plt.figure(figsize=(10, 5))
    plt.imshow(delta_t_worldclim, extent=[-180, 180, -60, 90], cmap='magma', aspect='auto')
    plt.title("Empirical WorldClim Present vs LGM Temperature Anomaly (|DeltaT| deg C)", fontweight='bold')
    plt.xlabel("Longitude (deg)")
    plt.ylabel("Latitude (deg)")
    plt.colorbar(label="Absolute Temperature Difference (deg C)", orientation='horizontal', shrink=0.7)
    plt.tight_layout()
    plt.savefig(config.MAPS_DIR / "empirical_map1_delta_t.png", dpi=300)
    plt.close()
    
    # MAP 2: GBIF Richness
    plt.figure(figsize=(10, 5))
    plt.imshow(species_richness_grid, extent=[-180, 180, -60, 90], cmap='plasma', aspect='auto', vmin=0, vmax=np.percentile(s_vals, 98))
    plt.title("Empirical GBIF Unique Species Richness (S)", fontweight='bold')
    plt.xlabel("Longitude (deg)")
    plt.ylabel("Latitude (deg)")
    plt.colorbar(label="Unique Species Count per Grid Cell", orientation='horizontal', shrink=0.7)
    plt.tight_layout()
    plt.savefig(config.MAPS_DIR / "empirical_map2_richness.png", dpi=300)
    plt.close()
    
    # SCATTER PLOT
    plt.figure(figsize=(8, 6))
    plt.scatter(dt_vals, log_s_vals, alpha=0.4, color='teal', s=15)
    plt.title("Empirical Species Richness vs. Palaeoclimate Temperature Anomaly (|DeltaT|)", fontweight='bold')
    plt.xlabel("WorldClim Absolute Temperature Change |DeltaT_LGM| (deg C)")
    plt.ylabel("Log(GBIF Unique Species Richness + 1)")
    plt.text(0.05, 0.90, f"Spearman rho = {rho_wc:.4f}\nPearson r = {r_wc:.4f}\nN = {n_cells:,} cells",
             transform=plt.gca().transAxes, bbox=dict(boxstyle="round", facecolor="white", alpha=0.8))
    plt.tight_layout()
    plt.savefig(config.MAPS_DIR / "empirical_scatter_richness_vs_deltat.png", dpi=300)
    plt.close()
    
    print(f"\n[COMPLETE] Analysis finished and saved to {config.MAPS_DIR} and {config.TABLES_DIR}")

if __name__ == "__main__":
    run_analysis()
