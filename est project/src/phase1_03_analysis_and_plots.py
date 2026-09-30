"""
phase1_03_analysis_and_plots.py
--------------------------------
Phase 1 Step 3: Statistical Analysis & Output Generation (Class Aves).
- Joins ΔT (Temperature Difference) and S (Unique GBIF Species Richness) per cell.
- Computes Spearman rank correlation (rho) and Pearson correlation (r).
- Renders:
  1. MAP 1: Present vs LGM Temperature Difference (|ΔT|)
  2. MAP 2: GBIF Unique Species Richness (S) for Class Aves
  3. SCATTER PLOT: Log(S + 1) vs Temperature Difference (|ΔT|)
- Exports exact statistical tables to results/tables/ and results/statistics/
"""

import sys
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib.pyplot as plt
import seaborn as sns

sys.path.append(str(Path(__file__).resolve().parent.parent))
import config

plt.style.use('seaborn-v0_8-whitegrid' if 'seaborn-v0_8-whitegrid' in plt.style.available else 'default')
plt.rcParams['font.size'] = 10
plt.rcParams['figure.dpi'] = 300

def render_phase1_map(data, title, filename, cmap='viridis', vmin=None, vmax=None, cbar_label=''):
    fig, ax = plt.subplots(figsize=(10, 5))
    im = ax.imshow(data, extent=[-180, 180, -90, 90], cmap=cmap, vmin=vmin, vmax=vmax, aspect='auto')
    ax.set_title(title, fontsize=12, fontweight='bold', pad=10)
    ax.set_xlabel('Longitude (°)')
    ax.set_ylabel('Latitude (°)')
    cbar = fig.colorbar(im, ax=ax, orientation='horizontal', pad=0.15, shrink=0.7)
    cbar.set_label(cbar_label)
    plt.tight_layout()
    output_path = config.MAPS_DIR / filename
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  [SAVED MAP] {filename}")

def run_phase1_analysis():
    print("=" * 80)
    print(f"PHASE 1: STATISTICAL ANALYSIS & MAP GENERATION ({config.TAXON_NAME.upper()})")
    print("=" * 80)
    
    # Load processed arrays
    delta_t_lgm = np.load(config.PROCESSED_DATA_DIR / "delta_t_lgm.npy")
    gbif_richness = np.load(config.PROCESSED_DATA_DIR / "gbif_species_richness.npy")
    gbif_records = np.load(config.PROCESSED_DATA_DIR / "gbif_record_density.npy")
    
    # Select terrestrial cells with richness > 0
    valid_mask = (gbif_richness > 0)
    
    x_delta_t = delta_t_lgm[valid_mask]
    y_richness = gbif_richness[valid_mask]
    y_log_richness = np.log1p(y_richness)
    n_samples = len(x_delta_t)
    
    # Calculate Associations
    rho_spearman, p_spearman = stats.spearmanr(x_delta_t, y_richness)
    r_pearson, p_pearson = stats.pearsonr(x_delta_t, y_log_richness)
    
    # Statistical Table
    stats_df = pd.DataFrame([{
        'Analysis Contrast': 'Present (1960-1990) vs LGM (~22ka CCSM4)',
        'Climate Variable': 'Annual Mean Temperature BIO1',
        'Biodiversity Metric': f'GBIF {config.TAXON_NAME} Unique Species Richness (S)',
        'Sample Size (Cells N)': n_samples,
        'Spearman rho': round(rho_spearman, 4),
        'Spearman p-value': "< 0.0001" if p_spearman < 0.0001 else f"{p_spearman:.4f}",
        'Pearson r (log S+1)': round(r_pearson, 4),
        'Pearson p-value': "< 0.0001" if p_pearson < 0.0001 else f"{p_pearson:.4f}",
    }])
    
    print("\nPHASE 1 STATISTICAL SUMMARY TABLE:")
    print(stats_df.to_string(index=False))
    
    table_path = config.TABLES_DIR / "phase1_statistical_results.csv"
    stats_df.to_csv(table_path, index=False)
    
    txt_path = config.STATS_DIR / "phase1_statistical_summary.txt"
    with open(txt_path, "w") as f:
        f.write(f"PHASE 1 MID-EVALUATION STATISTICAL REPORT ({config.TAXON_NAME.upper()})\n")
        f.write("=" * 60 + "\n\n")
        f.write(f"Confirmed Taxon: Class {config.TAXON_NAME} (TaxonKey = {config.TAXON_KEY})\n")
        f.write(f"Sample Size (Occupied 10' Cells): {n_samples}\n")
        f.write(f"Temperature Difference ΔT Mean: {x_delta_t.mean():.2f} °C (Std: {x_delta_t.std():.2f} °C)\n")
        f.write(f"Species Richness S Mean: {y_richness.mean():.2f} (Max: {y_richness.max()})\n\n")
        f.write(f"Spearman rank correlation (rho): {rho_spearman:.4f} (p = {p_spearman})\n")
        f.write(f"Pearson correlation (r, log(S+1)): {r_pearson:.4f} (p = {p_pearson})\n\n")
        f.write("Note: Observational correlation; does not imply direct causality.\n")
        
    # Render Required MAP 1: Present vs LGM ΔT
    render_phase1_map(
        delta_t_lgm,
        "MAP 1: Present vs Last Glacial Maximum (|ΔT_LGM|) Temperature Difference",
        "phase1_map1_delta_t_lgm.png",
        cmap='magma',
        cbar_label='Absolute Temperature Change (°C)'
    )
    
    # Render Required MAP 2: GBIF Unique Species Richness
    render_phase1_map(
        gbif_richness,
        f"MAP 2: GBIF Modern Class {config.TAXON_NAME} Unique Species Richness (S)",
        "phase1_map2_gbif_species_richness.png",
        cmap='plasma',
        cbar_label='Unique Species Count per 10\' Cell'
    )
    
    # Render Required SCATTER PLOT: Log(S+1) vs ΔT
    fig, ax = plt.subplots(figsize=(8, 6))
    sns.regplot(
        x=x_delta_t,
        y=y_log_richness,
        ax=ax,
        scatter_kws={'alpha': 0.3, 'color': 'teal', 's': 15},
        line_kws={'color': 'darkred', 'linewidth': 2}
    )
    ax.set_title(f"PHASE 1: Class {config.TAXON_NAME} Richness vs. Palaeoclimate Temp Change (|ΔT|)", fontsize=12, fontweight='bold')
    ax.set_xlabel("Absolute Temperature Change |ΔT_LGM| (°C) [Lower = More Stable]")
    ax.set_ylabel(f"Log(GBIF {config.TAXON_NAME} Unique Species Richness + 1)")
    
    # Annotate correlation statistics on plot
    ax.text(
        0.05, 0.92,
        f"Spearman $\\rho = {rho_spearman:.3f}$ ($p < 0.001$)\nPearson $r = {r_pearson:.3f}$ ($p < 0.001$)\n$N = {n_samples:,}$ cells",
        transform=ax.transAxes,
        fontsize=10,
        bbox=dict(boxstyle="round,pad=0.5", facecolor="white", alpha=0.8)
    )
    
    plt.tight_layout()
    scatter_path = config.MAPS_DIR / "phase1_scatter_richness_vs_deltat.png"
    plt.savefig(scatter_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"  [SAVED DIAGNOSTIC] phase1_scatter_richness_vs_deltat.png")
    
    print("\n[SUCCESS] Phase 1 analysis complete. All outputs generated in results/")

if __name__ == "__main__":
    run_phase1_analysis()
