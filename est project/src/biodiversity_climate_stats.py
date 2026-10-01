"""Descriptive biodiversity–climate statistical analysis on 5-degree grid."""

from __future__ import annotations

import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.colors as colors
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

import config

CAUSAL_NOTE = (
    "Descriptive cell-level association only; does not establish causation or "
    "account for unmeasured ecological confounders. Standard p-values ignore spatial autocorrelation."
)


def _spearman(y: np.ndarray, x: np.ndarray) -> tuple[float, float]:
    result = stats.spearmanr(y, x)
    return float(result.statistic), float(result.pvalue)


def _confidence_interval(r: float, n: int, confidence: float = 0.95) -> tuple[float, float]:
    """Calculate Fisher z-transform confidence interval for rank or Pearson correlation."""
    if n <= 3 or not np.isfinite(r) or abs(r) >= 1.0:
        return (np.nan, np.nan)
    z = np.arctanh(r)
    se = 1.0 / math.sqrt(n - 3)
    z_crit = stats.norm.ppf((1.0 + confidence) / 2.0)
    z_low = z - z_crit * se
    z_high = z + z_crit * se
    return (float(np.tanh(z_low)), float(np.tanh(z_high)))


def _partial_spearman(
    y: np.ndarray, x: np.ndarray, controls: list[np.ndarray]
) -> tuple[float, float]:
    """Standard partial rank correlation: Pearson r on OLS rank residuals with adjusted df."""
    ry = stats.rankdata(y)
    rx = stats.rankdata(x)
    control_ranks = [stats.rankdata(control) for control in controls]
    design = np.column_stack([np.ones(len(y)), *control_ranks])

    # Residualize ranks using OLS
    coef_y, _, _, _ = np.linalg.lstsq(design, ry, rcond=None)
    coef_x, _, _, _ = np.linalg.lstsq(design, rx, rcond=None)
    resid_y = ry - design @ coef_y
    resid_x = rx - design @ coef_x

    n = len(y)
    k = len(controls)
    df_resid = n - 2 - k

    r = float(stats.pearsonr(resid_y, resid_x).statistic)
    if df_resid > 0 and abs(r) < 1.0:
        t_stat = r * np.sqrt(df_resid / (1.0 - r**2))
        p_val = float(2.0 * (1.0 - stats.t.cdf(abs(t_stat), df=df_resid)))
    else:
        p_val = np.nan
    return r, p_val


def _queen_weights(lat_index: np.ndarray, lon_index: np.ndarray) -> np.ndarray:
    n = len(lat_index)
    index = {(int(lat_index[i]), int(lon_index[i])): i for i in range(n)}
    weights = np.zeros((n, n), dtype="float64")
    for i in range(n):
        lat = int(lat_index[i])
        lon = int(lon_index[i])
        for dlat in (-1, 0, 1):
            for dlon in (-1, 0, 1):
                if dlat == 0 and dlon == 0:
                    continue
                neighbor = index.get((lat + dlat, lon + dlon))
                if neighbor is not None:
                    weights[i, neighbor] = 1.0
    # Row-standardize weights
    row_sums = weights.sum(axis=1, keepdims=True)
    row_sums[row_sums == 0] = 1.0
    return weights / row_sums


def _morans_i(values: np.ndarray, weights: np.ndarray, permutations: int = 999) -> dict:
    x = np.asarray(values, dtype="float64")
    n = len(x)
    if n < 3 or weights.sum() == 0:
        return {"morans_i": np.nan, "morans_i_p_value": np.nan, "morans_i_permutations": 0}

    xc = x - x.mean()
    den = np.dot(xc, xc)
    if den == 0:
        return {"morans_i": np.nan, "morans_i_p_value": np.nan, "morans_i_permutations": 0}

    s0 = weights.sum()
    observed = (n / s0) * np.dot(xc, weights @ xc) / den

    rng = np.random.default_rng(0)
    permuted = np.empty(permutations, dtype="float64")
    for iteration in range(permutations):
        xp = rng.permutation(x)
        xpc = xp - xp.mean()
        d = np.dot(xpc, xpc)
        if d == 0:
            permuted[iteration] = observed
        else:
            permuted[iteration] = (n / s0) * np.dot(xpc, weights @ xpc) / d

    p_value = float((np.sum(np.abs(permuted) >= abs(observed)) + 1) / (permutations + 1))
    return {
        "morans_i": float(observed),
        "morans_i_p_value": p_value,
        "morans_i_permutations": permutations,
    }


def _load_5deg_dataset() -> pd.DataFrame:
    if not config.CLIMATE_STABILITY_CELLS_CSV.exists():
        raise FileNotFoundError(config.CLIMATE_STABILITY_CELLS_CSV)
    if not config.GBIF_SAMPLING_EFFORT_5DEG.exists():
        raise FileNotFoundError(config.GBIF_SAMPLING_EFFORT_5DEG)

    climate = pd.read_csv(config.CLIMATE_STABILITY_CELLS_CSV)
    effort = pd.read_csv(config.GBIF_SAMPLING_EFFORT_5DEG)
    merged = climate.merge(
        effort,
        on="grid_cell_id",
        how="inner",
        suffixes=("_climate", "_gbif"),
    )
    return merged


def _render_plots(df: pd.DataFrame, sens: pd.DataFrame, stats_summary: dict) -> None:
    config.MAPS_DIR.mkdir(parents=True, exist_ok=True)

    # 1. Main Scatter Plot with Log-Scale Richness and Dynamic Statistics Text
    fig, ax = plt.subplots(figsize=(10, 6.5))
    sc = ax.scatter(
        df["temperature_stability_index"],
        df["species_richness"],
        c=np.log10(df["occurrence_count"]),
        cmap="viridis",
        s=40,
        alpha=0.85,
        edgecolors="none",
    )
    ax.set_yscale("log")
    cb = fig.colorbar(sc, ax=ax, pad=0.02)
    cb.set_label("Sampling Effort: log₁₀(GBIF occurrence count per cell)", fontsize=11)

    # Construct statistics text dynamically from computed metrics
    stats_text = (
        f"Primary Spearman Test:\n"
        f"  n = {stats_summary['n']}\n"
        f"  ρ = {stats_summary['rho_primary']:+.4f} (p = {stats_summary['p_primary']:.4f})\n\n"
        f"Partial Correlations:\n"
        f"  Control effort: r = {stats_summary['r_eff']:+.4f} (p = {stats_summary['p_eff']:.4f})\n"
        f"  Control effort + |lat|: r = {stats_summary['r_both']:+.4f} (p = {stats_summary['p_both']:.4f})\n\n"
        f"Regional Subsets:\n"
        f"  Extratropics (|lat|>23.5°): ρ = {stats_summary['r_extra']:+.4f} (p = {stats_summary['p_extra']:.4f})\n"
        f"  Tropics (|lat|≤23.5°): ρ = {stats_summary['r_trop']:+.4f} (p = {stats_summary['p_trop']:.4f})\n\n"
        f"Spatial Autocorrelation:\n"
        f"  Moran's I (richness) = {stats_summary['moran_i_rich']:+.4f} (p = {stats_summary['moran_p_rich']:.4f})\n"
        f"  Moran's I (stability) = {stats_summary['moran_i_stab']:+.4f} (p = {stats_summary['moran_p_stab']:.4f})"
    )
    ax.text(
        0.03,
        0.97,
        stats_text,
        transform=ax.transAxes,
        fontsize=9,
        verticalalignment="top",
        bbox=dict(boxstyle="round,pad=0.5", facecolor="white", alpha=0.85, edgecolor="#cccccc"),
    )

    ax.set_xlabel("Palaeoclimate Temperature Stability Index: 1 / (1 + SD(BIO1) °C)", fontsize=12)
    ax.set_ylabel("GBIF Unique Species Richness (log scale)", fontsize=12)
    ax.set_title(
        "Modern Biodiversity Richness vs. Palaeoclimate Stability (5° Grid Cells)\n"
        "PaleoClim BIO1 (Current, Late Holocene, LGM) and GBIF Observations",
        fontsize=13,
        pad=12,
    )
    fig.tight_layout()
    fig.savefig(config.MAIN_SCATTER_PNG, dpi=200)
    plt.close(fig)

    # 2. Threshold Sensitivity Plot with 95% Confidence Intervals
    fig, ax = plt.subplots(figsize=(8.5, 5.5))
    thresholds = sens["min_occurrence_threshold"]
    rho_val = sens["spearman_rho"]
    rho_low = sens["spearman_rho_ci_lower"]
    rho_high = sens["spearman_rho_ci_upper"]
    part_val = sens["partial_r_control_effort_and_abs_lat"]
    part_low = sens["partial_r_ci_lower"]
    part_high = sens["partial_r_ci_upper"]

    ax.plot(thresholds, rho_val, marker="o", color="#1f77b4", linewidth=2, label="Raw Spearman ρ")
    ax.fill_between(thresholds, rho_low, rho_high, color="#1f77b4", alpha=0.15)

    ax.plot(
        thresholds,
        part_val,
        marker="s",
        color="#2ca02c",
        linewidth=2,
        linestyle="--",
        label="Partial r (control effort + |lat|)",
    )
    ax.fill_between(thresholds, part_low, part_high, color="#2ca02c", alpha=0.15)

    ax.axhline(0, color="gray", linestyle=":", linewidth=1)
    for _, row in sens.iterrows():
        th = int(row["min_occurrence_threshold"])
        n_c = int(row["n_cells"])
        rho = row["spearman_rho"]
        if np.isfinite(rho):
            ax.annotate(
                f"n={n_c}",
                (th, rho),
                textcoords="offset points",
                xytext=(0, 8),
                ha="center",
                fontsize=8,
            )

    ax.set_xlabel("Minimum Occurrence Count Threshold (records per 5° cell)", fontsize=11)
    ax.set_ylabel("Correlation with Temperature Stability", fontsize=11)
    ax.set_title("Sensitivity of Stability–Biodiversity Correlation with 95% CIs", fontsize=12)
    ax.grid(True, linestyle="--", alpha=0.5)
    ax.legend(loc="lower right")
    fig.tight_layout()
    fig.savefig(config.MAPS_DIR / "biodiversity_vs_stability_by_threshold.png", dpi=200)
    plt.close(fig)

    # 3. 5° Species Richness Map using Log Scale
    fig, ax = plt.subplots(figsize=(12, 6))
    grid = np.full((36, 72), np.nan, dtype="float32")
    for _, row in df.iterrows():
        lat_idx = int(row["grid_lat_index_climate"])
        lon_idx = int(row["grid_lon_index_climate"])
        raster_row = 35 - lat_idx
        grid[raster_row, lon_idx] = float(row["species_richness"])

    masked = np.ma.masked_invalid(grid)
    vmax = float(np.nanmax(grid)) if np.any(np.isfinite(grid)) else 100.0
    image = ax.imshow(
        masked,
        extent=(-180.0, 180.0, -90.0, 90.0),
        origin="upper",
        aspect="auto",
        cmap="plasma",
        norm=colors.LogNorm(vmin=1.0, vmax=vmax),
    )
    ax.set(
        title="GBIF Observed Unique Species Richness per 5° Grid Cell (Log Scale)",
        xlabel="Longitude (°)",
        ylabel="Latitude (°)",
    )
    fig.colorbar(image, ax=ax, orientation="horizontal", label="Unique Species Count (Log Scale)", pad=0.12)
    fig.tight_layout()
    fig.savefig(config.MAPS_DIR / "gbif_species_richness_5deg.png", dpi=200)
    plt.close(fig)


def run_analysis() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    df = _load_5deg_dataset()
    config.RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    config.TABLES_DIR.mkdir(parents=True, exist_ok=True)

    # Save final merged 5-degree analysis table
    final_cols = [
        "grid_cell_id",
        "grid_lat_index_climate",
        "grid_lon_index_climate",
        "cell_latitude_center_climate",
        "cell_longitude_center_climate",
        "cell_area_km2",
        "occurrence_count",
        "species_richness",
        "richness_rarefied_5",
        "richness_rarefied_10",
        "bio1_current_c",
        "bio1_late_holocene_c",
        "bio1_lgm_c",
        "bio1_sd_c",
        "temperature_stability_index",
    ]
    available_cols = [c for c in final_cols if c in df.columns]
    df_final = df[available_cols].copy()
    df_final.to_csv(config.FINAL_ANALYSIS_CSV, index=False)
    df_final.to_csv(config.RESULTS_DIR / "final_biodiversity_climate_5deg.csv", index=False)

    y = df["species_richness"].to_numpy(dtype="float64")
    x = df["temperature_stability_index"].to_numpy(dtype="float64")
    sd = df["bio1_sd_c"].to_numpy(dtype="float64")
    effort = df["occurrence_count"].to_numpy(dtype="float64")
    signed_lat = df["cell_latitude_center_climate"].to_numpy(dtype="float64")
    abs_lat = np.abs(signed_lat)

    # 1. Primary & Sub-model Statistics
    stats_rows = []

    # Primary Spearman test
    rho_prim, p_prim = _spearman(y, x)
    stats_rows.append(
        {
            "analysis_scale": "5deg_cell",
            "climate_variable": "temperature_stability_index",
            "response_variable": "species_richness",
            "model": "spearman_primary",
            "primary_test": True,
            "n": len(df),
            "correlation_coefficient": rho_prim,
            "p_value": p_prim,
            "covariates_controlled": "none",
            "interpretation": "Primary test: unadjusted bivariate rank association across all occupied 5° cells. Note p-values do not adjust for spatial autocorrelation.",
        }
    )

    # Primary Spearman test for cells with >=5 records (weak cells filtered)
    mask_5rec = effort >= 5
    if mask_5rec.sum() > 3:
        rho_5rec, p_5rec = _spearman(y[mask_5rec], x[mask_5rec])
        stats_rows.append(
            {
                "analysis_scale": "5deg_cell_min5rec",
                "climate_variable": "temperature_stability_index",
                "response_variable": "species_richness",
                "model": "spearman_min5_records",
                "primary_test": False,
                "n": int(mask_5rec.sum()),
                "correlation_coefficient": rho_5rec,
                "p_value": p_5rec,
                "covariates_controlled": "filter: occurrence_count >= 5",
                "interpretation": "Primary test restricted to cells with >=5 records, excluding low-effort fringe cells.",
            }
        )

    # Rarefied richness (10 records) test
    if "richness_rarefied_10" in df.columns:
        valid_rar10 = df["richness_rarefied_10"].notna()
        if valid_rar10.sum() > 3:
            y_rar10 = df.loc[valid_rar10, "richness_rarefied_10"].to_numpy(dtype="float64")
            x_rar10 = df.loc[valid_rar10, "temperature_stability_index"].to_numpy(dtype="float64")
            rho_rar10, p_rar10 = _spearman(y_rar10, x_rar10)
            stats_rows.append(
                {
                    "analysis_scale": "5deg_cell_rarefied10",
                    "climate_variable": "temperature_stability_index",
                    "response_variable": "richness_rarefied_10",
                    "model": "spearman_rarefied_10_records",
                    "primary_test": False,
                    "n": int(valid_rar10.sum()),
                    "correlation_coefficient": rho_rar10,
                    "p_value": p_rar10,
                    "covariates_controlled": "fixed_sample_rarefaction_10_records",
                    "interpretation": "Fixed-sample rarefied richness (10 records per cell) to control effort variation.",
                }
            )

    # Stability SD contrast
    rho_sd, p_sd = _spearman(y, sd)
    stats_rows.append(
        {
            "analysis_scale": "5deg_cell",
            "climate_variable": "bio1_sd_c",
            "response_variable": "species_richness",
            "model": "spearman_bio1_sd",
            "primary_test": False,
            "n": len(df),
            "correlation_coefficient": rho_sd,
            "p_value": p_sd,
            "covariates_controlled": "none",
            "interpretation": "Contrast with temperature standard deviation (inverse stability).",
        }
    )

    # Partial: control effort
    r_eff, p_eff = _partial_spearman(y, x, [effort])
    stats_rows.append(
        {
            "analysis_scale": "5deg_cell",
            "climate_variable": "temperature_stability_index",
            "response_variable": "species_richness",
            "model": "partial_spearman_control_occurrence_count",
            "primary_test": False,
            "n": len(df),
            "correlation_coefficient": r_eff,
            "p_value": p_eff,
            "covariates_controlled": "occurrence_count",
            "interpretation": "Controls for sampling effort (record count per cell); adjusts for collection bias.",
        }
    )

    # Partial: control absolute latitude
    r_alat, p_alat = _partial_spearman(y, x, [abs_lat])
    stats_rows.append(
        {
            "analysis_scale": "5deg_cell",
            "climate_variable": "temperature_stability_index",
            "response_variable": "species_richness",
            "model": "partial_spearman_control_abs_latitude",
            "primary_test": False,
            "n": len(df),
            "correlation_coefficient": r_alat,
            "p_value": p_alat,
            "covariates_controlled": "abs(latitude)",
            "interpretation": "Controls for latitudinal diversity gradient (equator-to-pole distance).",
        }
    )

    # Partial: control effort + absolute latitude
    r_both, p_both = _partial_spearman(y, x, [effort, abs_lat])
    stats_rows.append(
        {
            "analysis_scale": "5deg_cell",
            "climate_variable": "temperature_stability_index",
            "response_variable": "species_richness",
            "model": "partial_spearman_control_effort_and_abs_latitude",
            "primary_test": False,
            "n": len(df),
            "correlation_coefficient": r_both,
            "p_value": p_both,
            "covariates_controlled": "occurrence_count, abs(latitude)",
            "interpretation": "Simultaneously controls for sampling effort and the latitudinal diversity gradient.",
        }
    )

    # Regional subset: Extratropics
    extra_mask = abs_lat > 23.5
    r_extra, p_extra = _spearman(y[extra_mask], x[extra_mask])
    stats_rows.append(
        {
            "analysis_scale": "5deg_cell_extratropics",
            "climate_variable": "temperature_stability_index",
            "response_variable": "species_richness",
            "model": "spearman_extratropics",
            "primary_test": False,
            "n": int(extra_mask.sum()),
            "correlation_coefficient": r_extra,
            "p_value": p_extra,
            "covariates_controlled": "subregion: |lat| > 23.5°",
            "interpretation": "Bivariate Spearman correlation restricted to extratropical cells.",
        }
    )

    # Regional subset: Tropics
    trop_mask = abs_lat <= 23.5
    r_trop, p_trop = _spearman(y[trop_mask], x[trop_mask])
    stats_rows.append(
        {
            "analysis_scale": "5deg_cell_tropics",
            "climate_variable": "temperature_stability_index",
            "response_variable": "species_richness",
            "model": "spearman_tropics",
            "primary_test": False,
            "n": int(trop_mask.sum()),
            "correlation_coefficient": r_trop,
            "p_value": p_trop,
            "covariates_controlled": "subregion: |lat| <= 23.5°",
            "interpretation": "Bivariate Spearman correlation restricted to tropical cells.",
        }
    )

    statistics = pd.DataFrame(stats_rows)
    statistics.to_csv(config.BIODIV_CLIMATE_STATS_CSV, index=False)

    # 2. Sampling Effort / Record Threshold Sensitivity with 95% CIs
    sens_rows = []
    thresholds = [1, 2, 5, 10, 20, 50, 100]
    for th in thresholds:
        mask = df["occurrence_count"] >= th
        sub = df.loc[mask]
        n_c = len(sub)
        sub_y = sub["species_richness"].to_numpy(dtype="float64")
        sub_x = sub["temperature_stability_index"].to_numpy(dtype="float64")
        sub_eff = sub["occurrence_count"].to_numpy(dtype="float64")
        sub_alat = np.abs(sub["cell_latitude_center_climate"].to_numpy(dtype="float64"))

        r_th, p_th = _spearman(sub_y, sub_x)
        rho_low, rho_high = _confidence_interval(r_th, n_c)
        r_part_both, p_part_both = _partial_spearman(sub_y, sub_x, [sub_eff, sub_alat])
        part_low, part_high = _confidence_interval(r_part_both, n_c)

        sens_rows.append(
            {
                "min_occurrence_threshold": th,
                "n_cells": n_c,
                "pct_cells_retained": float(n_c / len(df) * 100.0),
                "spearman_rho": r_th,
                "spearman_rho_ci_lower": rho_low,
                "spearman_rho_ci_upper": rho_high,
                "spearman_p_value": p_th,
                "partial_r_control_effort_and_abs_lat": r_part_both,
                "partial_r_ci_lower": part_low,
                "partial_r_ci_upper": part_high,
                "partial_p_control_effort_and_abs_lat": p_part_both,
            }
        )

    sensitivity = pd.DataFrame(sens_rows)
    sensitivity.to_csv(config.BIODIV_CLIMATE_SENSITIVITY_CSV, index=False)
    sensitivity.to_csv(config.TABLES_DIR / "biodiversity_climate_sensitivity_thresholds.csv", index=False)

    # 3. Moran's I Spatial Autocorrelation with Dynamic Interpretation
    lat_idx = df["grid_lat_index_climate"].to_numpy(dtype=int)
    lon_idx = df["grid_lon_index_climate"].to_numpy(dtype=int)
    weights = _queen_weights(lat_idx, lon_idx)

    # OLS rank residuals for spatial autocorrelation of model errors
    ry = stats.rankdata(y)
    rx = stats.rankdata(x)
    slope, intercept, _, _, _ = stats.linregress(rx, ry)
    rank_residuals = ry - (intercept + slope * rx)

    moran_specs = [
        ("species_richness", y),
        ("temperature_stability_index", x),
        ("rank_regression_residuals", rank_residuals),
    ]

    moran_rows = []
    for var_name, values in moran_specs:
        m_res = _morans_i(values, weights, permutations=999)
        mi = m_res["morans_i"]
        pv = m_res["morans_i_p_value"]
        if np.isfinite(mi):
            direction = "Positive" if mi > 0 else "Negative"
            sig_text = f"(p = {pv:.4f})" if pv > 0.001 else "(p <= 0.001)"
            interp = f"{direction} spatial autocorrelation detected {sig_text}; spatial clustering indicates non-independence of cells."
        else:
            interp = "Spatial autocorrelation could not be calculated."

        moran_rows.append(
            {
                "variable": var_name,
                "n": len(df),
                "weight_scheme": "row_standardized_queen_contiguity_5deg",
                **m_res,
                "interpretation": interp,
            }
        )

    moran_df = pd.DataFrame(moran_rows)
    moran_df.to_csv(config.BIODIV_CLIMATE_MORAN_CSV, index=False)

    # 4. Render plots using dynamic statistics summary
    stats_summary = {
        "n": len(df),
        "rho_primary": rho_prim,
        "p_primary": p_prim,
        "r_eff": r_eff,
        "p_eff": p_eff,
        "r_both": r_both,
        "p_both": p_both,
        "r_extra": r_extra,
        "p_extra": p_extra,
        "r_trop": r_trop,
        "p_trop": p_trop,
        "moran_i_rich": float(moran_df.loc[moran_df["variable"] == "species_richness", "morans_i"].values[0]),
        "moran_p_rich": float(moran_df.loc[moran_df["variable"] == "species_richness", "morans_i_p_value"].values[0]),
        "moran_i_stab": float(moran_df.loc[moran_df["variable"] == "temperature_stability_index", "morans_i"].values[0]),
        "moran_p_stab": float(moran_df.loc[moran_df["variable"] == "temperature_stability_index", "morans_i_p_value"].values[0]),
    }
    _render_plots(df, sensitivity, stats_summary)

    return statistics, sensitivity, moran_df
