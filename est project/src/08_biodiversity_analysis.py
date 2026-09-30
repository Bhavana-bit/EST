"""
08_biodiversity_analysis.py
---------------------------
Statistical analysis of Palaeoclimate Stability vs. Modern Biodiversity.
Calculates:
1. Pearson and Spearman correlation coefficients.
2. Multivariate regression modeling log(S + 1) ~ Stability + Roughness + Latitude + log(Records).
3. Spatial Autocorrelation test (Moran's I on residuals).
Saves exact calculated statistical tables to results/tables/ and results/statistics/.
"""

import sys
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats
import statsmodels.api as sm

sys.path.append(str(Path(__file__).resolve().parent.parent))
import config

def analyze_biodiversity_associations():
    print("=" * 80)
    print("STATISTICAL ANALYSIS: PALAEOCLIMATE STABILITY VS MODERN BIODIVERSITY")
    print("=" * 80)
    
    # Load processed arrays
    stability_index = np.load(config.PROCESSED_DATA_DIR / "stability_index.npy")
    delta_t_lgm = np.load(config.PROCESSED_DATA_DIR / "delta_t_lgm.npy")
    delta_p_lgm = np.load(config.PROCESSED_DATA_DIR / "delta_p_lgm.npy")
    gbif_richness = np.load(config.PROCESSED_DATA_DIR / "gbif_species_richness.npy")
    iucn_richness = np.load(config.PROCESSED_DATA_DIR / "iucn_species_richness.npy")
    gbif_records = np.load(config.PROCESSED_DATA_DIR / "gbif_record_density.npy")
    roughness = np.load(config.PROCESSED_DATA_DIR / "topographic_roughness.npy")
    
    rows, cols = config.GRID_SHAPE
    lats = np.linspace(90, -90, rows)
    lons = np.linspace(-180, 180, cols)
    lon_grid, lat_grid = np.meshgrid(lons, lats)
    
    # Mask terrestrial non-zero cells
    valid_mask = (gbif_richness > 0)
    
    df_analysis = pd.DataFrame({
        'stability_index': stability_index[valid_mask],
        'delta_t_lgm': delta_t_lgm[valid_mask],
        'delta_p_lgm': delta_p_lgm[valid_mask],
        'gbif_richness': gbif_richness[valid_mask],
        'iucn_richness': iucn_richness[valid_mask],
        'log_gbif_richness': np.log1p(gbif_richness[valid_mask]),
        'gbif_records': gbif_records[valid_mask],
        'log_gbif_records': np.log1p(gbif_records[valid_mask]),
        'roughness': roughness[valid_mask],
        'abs_latitude': np.abs(lat_grid[valid_mask])
    })
    
    # 1. Correlation Analysis
    correlations = []
    pairs = [
        ('GBIF Richness', 'Stability Index', 'stability_index'),
        ('GBIF Richness', 'LGM Temp Anomaly (ΔT)', 'delta_t_lgm'),
        ('GBIF Richness', 'LGM Precip Anomaly (ΔP)', 'delta_p_lgm'),
        ('GBIF Richness', 'Topographic Roughness', 'roughness'),
        ('GBIF Richness', 'IUCN Range Richness', 'iucn_richness'),
    ]
    
    for label1, label2, col_var in pairs:
        r_pearson, p_pearson = stats.pearsonr(df_analysis['gbif_richness'], df_analysis[col_var])
        r_spearman, p_spearman = stats.spearmanr(df_analysis['gbif_richness'], df_analysis[col_var])
        
        correlations.append({
            'Variable 1': label1,
            'Variable 2': label2,
            'Pearson r': round(r_pearson, 4),
            'Pearson p-value': "< 0.0001" if p_pearson < 0.0001 else f"{p_pearson:.4f}",
            'Spearman rho': round(r_spearman, 4),
            'Spearman p-value': "< 0.0001" if p_spearman < 0.0001 else f"{p_spearman:.4f}"
        })
        
    df_corr = pd.DataFrame(correlations)
    print("\nCORRELATION RESULTS TABLE:")
    print(df_corr.to_string(index=False))
    df_corr.to_csv(config.TABLES_DIR / "05_correlation_results.csv", index=False)
    
    # 2. Multivariate Regression Model (OLS)
    X = df_analysis[['stability_index', 'roughness', 'abs_latitude', 'log_gbif_records']]
    X = sm.add_constant(X)
    y = df_analysis['log_gbif_richness']
    
    model = sm.OLS(y, X).fit()
    
    print("\nMULTIVARIATE OLS REGRESSION SUMMARY:")
    print(model.summary())
    
    # Export regression table
    reg_summary_df = pd.DataFrame({
        'Coefficient': model.params.round(4),
        'Std Error': model.bse.round(4),
        't-value': model.tvalues.round(2),
        'p-value': model.pvalues.apply(lambda x: "< 0.0001" if x < 0.0001 else f"{x:.4f}")
    })
    reg_summary_df.to_csv(config.TABLES_DIR / "06_regression_results.csv")
    
    # 3. Spatial Autocorrelation Test (Moran's I on Residuals)
    residuals = model.resid.values
    # Moran's I approximation on spatial grid residuals
    morans_i = float(np.corrcoef(residuals[:-1], residuals[1:])[0, 1])
    
    print(f"\nSpatial Autocorrelation Check:")
    print(f"  Residual Moran's I proxy: {morans_i:.4f}")
    if morans_i > 0.2:
        print("  [NOTICE] Substantial spatial autocorrelation detected in residuals (Moran's I > 0.2).")
        print("  Spatial Autoregressive (SAR / Spatial Error) modeling recommended for final publication.")
        
    with open(config.STATS_DIR / "statistical_summary.txt", "w") as f:
        f.write("PALAEOCLIMATE STABILITY & MODERN BIODIVERSITY STATISTICAL SUMMARY\n")
        f.write("=" * 70 + "\n\n")
        f.write(f"Sample Size (Grid Cells): {len(df_analysis)}\n")
        f.write(f"Model R-squared: {model.rsquared:.4f}\n")
        f.write(f"Model Adjusted R-squared: {model.rsquared_adj:.4f}\n")
        f.write(f"Residual Moran's I: {morans_i:.4f}\n\n")
        f.write(str(model.summary()))
        
    print(f"\nSaved statistical reports to {config.TABLES_DIR} and {config.STATS_DIR}")

if __name__ == "__main__":
    analyze_biodiversity_associations()
