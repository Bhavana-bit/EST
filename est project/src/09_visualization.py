"""
09_visualization.py
-------------------
Generates the 10 required scientific publication maps & diagnostic scatter plots.
All maps include titles, legends/colorbars, geographic axes, and correct colormaps.
"""

import sys
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

sys.path.append(str(Path(__file__).resolve().parent.parent))
import config

# Set scientific plotting aesthetics
plt.style.use('seaborn-v0_8-whitegrid' if 'seaborn-v0_8-whitegrid' in plt.style.available else 'default')
plt.rcParams['font.size'] = 10
plt.rcParams['figure.dpi'] = 300

def render_map(data, title, filename, cmap='viridis', vmin=None, vmax=None, cbar_label=''):
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

def generate_all_plots():
    print("=" * 80)
    print("GENERATING ALL 10 REQUIRED SCIENTIFIC MAPS & DIAGNOSTIC PLOTS")
    print("=" * 80)
    
    # Load processed arrays
    bio1_present = np.load(config.PROCESSED_DATA_DIR / "bio1_present.npy")
    bio1_lgm = np.load(config.PROCESSED_DATA_DIR / "bio1_lgm.npy")
    bio1_mh = np.load(config.PROCESSED_DATA_DIR / "bio1_mh.npy")
    
    delta_t_lgm = np.load(config.PROCESSED_DATA_DIR / "delta_t_lgm.npy")
    delta_t_mh = np.load(config.PROCESSED_DATA_DIR / "delta_t_mh.npy")
    delta_p_lgm = np.load(config.PROCESSED_DATA_DIR / "delta_p_lgm.npy")
    delta_p_mh = np.load(config.PROCESSED_DATA_DIR / "delta_p_mh.npy")
    
    stability_index = np.load(config.PROCESSED_DATA_DIR / "stability_index.npy")
    gbif_richness = np.load(config.PROCESSED_DATA_DIR / "gbif_species_richness.npy")
    
    # 1. Present BIO1
    render_map(bio1_present, "WorldClim Present Annual Mean Temperature (1960-1990)", "01_present_bio1.png", cmap='coolwarm', cbar_label='Temperature (°C)')
    
    # 2. LGM BIO1
    render_map(bio1_lgm, "WorldClim LGM (~22ka) Annual Mean Temperature (CCSM4)", "02_lgm_bio1.png", cmap='coolwarm', cbar_label='Temperature (°C)')
    
    # 3. Mid-Holocene BIO1
    render_map(bio1_mh, "WorldClim Mid-Holocene (~6ka) Annual Mean Temperature (CCSM4)", "03_mh_bio1.png", cmap='coolwarm', cbar_label='Temperature (°C)')
    
    # 4. Present vs LGM Temp Difference
    render_map(delta_t_lgm, "Temperature Difference: Present vs Last Glacial Maximum (|ΔT_LGM|)", "04_delta_t_lgm.png", cmap='magma', cbar_label='Absolute Temp Change (°C)')
    
    # 5. Present vs Mid-Holocene Temp Difference
    render_map(delta_t_mh, "Temperature Difference: Present vs Mid-Holocene (|ΔT_MH|)", "05_delta_t_mh.png", cmap='magma', cbar_label='Absolute Temp Change (°C)')
    
    # 6. Present vs LGM Precip Difference
    render_map(delta_p_lgm, "Precipitation Difference: Present vs Last Glacial Maximum (|ΔP_LGM|)", "06_delta_p_lgm.png", cmap='viridis', cbar_label='Absolute Precip Change (mm/yr)')
    
    # 7. Present vs Mid-Holocene Precip Difference
    render_map(delta_p_mh, "Precipitation Difference: Present vs Mid-Holocene (|ΔP_MH|)", "07_delta_p_mh.png", cmap='viridis', cbar_label='Absolute Precip Change (mm/yr)')
    
    # 8. Climate Stability Index
    render_map(stability_index, "Project Palaeoclimate Stability Index (0 = Instability, 1 = Stability)", "08_climate_stability_index.png", cmap='YlGnBu', cbar_label='Stability Index')
    
    # 9. GBIF Species Richness
    render_map(gbif_richness, "GBIF Mammal Species Richness (Unique Species per 10' Cell)", "09_gbif_species_richness.png", cmap='plasma', cbar_label='Species Richness')
    
    # 10. Climate Stability vs Species Richness Scatter Plot
    fig, ax = plt.subplots(figsize=(8, 6))
    valid = (gbif_richness > 0)
    x_vals = stability_index[valid]
    y_vals = np.log1p(gbif_richness[valid])
    
    sns.regplot(x=x_vals, y=y_vals, ax=ax, scatter_kws={'alpha': 0.3, 'color': 'teal', 's': 15}, line_kws={'color': 'darkred', 'linewidth': 2})
    ax.set_title("Modern Species Richness vs. Palaeoclimate Stability", fontsize=12, fontweight='bold')
    ax.set_xlabel("Palaeoclimate Stability Index (Higher = More Stable)")
    ax.set_ylabel("Log(GBIF Species Richness + 1)")
    plt.tight_layout()
    plt.savefig(config.MAPS_DIR / "10_stability_vs_richness_scatter.png", dpi=300)
    plt.close()
    print("  [SAVED DIAGNOSTIC] 10_stability_vs_richness_scatter.png")
    
    print(f"\n[COMPLETE] All 10 plots generated successfully in {config.MAPS_DIR}")

if __name__ == "__main__":
    generate_all_plots()
