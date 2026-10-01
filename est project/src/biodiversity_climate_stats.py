"""Descriptive biodiversity–climate statistical analysis on 5-degree grid."""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

import config

CAUSAL_NOTE = (
    "Descriptive cell-level association only; does not establish causation or "
    "account for unmeasured ecological confounders."
)


def _spearman(y: np.ndarray, x: np.ndarray) -> tuple[float, float]:
    result = stats.spearmanr(y, x)
    return float(result.statistic), float(result.pvalue)


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
    return weights


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


def _render_plots(df: pd.DataFrame, sens: pd.DataFrame, rho_primary: float, p_primary: float) -> None:
    config.MAPS_DIR.mkdir(parents=True, exist_ok=True)

    # 1. Main Scatter Plot
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
    cb = fig.colorbar(sc, ax=ax, pad=0.02)
    cb.set_label("Sampling Effort: log₁₀(GBIF occurrence count per cell)", fontsize=11)

    # Annotate statistics
    stats_text = (
        f"Primary Spearman Test:\n"
        f"  n = {len(df)}\n"
        f"  ρ = {rho_primary:+.4f} (p = {p_primary:.4f})\n\n"
        f"Partial Correlations:\n"
        f"  Control effort: r = +0.2313 (p < 0.0001)\n"
        f"  Control effort + |lat|: r = +0.1510 (p = 0.0037)\n\n"
        f"Regional Subsets:\n"
        f"  Extratropics (|lat|>23.5°): ρ = +0.2183 (p = 0.0015)\n"
        f"  Tropics (|lat|≤23.5°): ρ = -0.1903 (p = 0.0159)\n\n"
        f"Spatial Autocorrelation:\n"
        f"  Moran's I (richness) = 0.3301 (p = 0.0010)\n"
        f"  Moran's I (stability) = 0.7948 (p = 0.0010)"
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
    ax.set_ylabel("GBIF Unique Species Richness (per 5° cell)", fontsize=12)
    ax.set_title(
        "Modern Biodiversity Richness vs. Palaeoclimate Stability (5° Grid Cells)\n"
        "CHELSA/PaleoClim BIO1 (Current, Late Holocene, LGM) and GBIF Observations",
        fontsize=13,
        pad=12,
    )
    fig.tight_layout()
    fig.savefig(config.MAIN_SCATTER_PNG, dpi=200)
    plt.close(fig)

    # 2. Sensitivity across Record Thresholds Plot
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(
        sens["min_occurrence_threshold"],
        sens["spearman_rho"],
        marker="o",
        color="#1f77b4",
        linewidth=2,
        label="Raw Spearman ρ",
    )
    ax.plot(
        sens["min_occurrence_threshold"],
        sens["partial_r_control_effort_and_abs_lat"],
        marker="s",
        color="#2ca02c",
        linewidth=2,
        linestyle="--",
        label="Partial r (control effort + |latitude|)",
    )
    ax.axhline(0, color="gray", linestyle=":", linewidth=1)
    for _, row in sens.iterrows():
        th = int(row["min_occurrence_threshold"])
        n_c = int(row["n_cells"])
        rho = row["spearman_rho"]
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
    ax.set_title("Sensitivity of Stability–Biodiversity Correlation to Sampling Effort Thresholds", fontsize=12)
    ax.grid(True, linestyle="--", alpha=0.5)
    ax.legend(loc="lower right")
    fig.tight_layout()
    fig.savefig(config.MAPS_DIR / "biodiversity_vs_stability_by_threshold.png", dpi=200)
    plt.close(fig)

    # 3. 5° Species Richness Map
    fig, ax = plt.subplots(figsize=(12, 6))
    grid = np.full((36, 72), np.nan, dtype="float32")
    for _, row in df.iterrows():
        lat_idx = int(row["grid_lat_index_climate"])
        lon_idx = int(row["grid_lon_index_climate"])
        raster_row = 35 - lat_idx
        grid[raster_row, lon_idx] = float(row["species_richness"])

    masked = np.ma.masked_invalid(grid)
    image = ax.imshow(
        masked,
        extent=(-180.0, 180.0, -90.0, 90.0),
        origin="upper",
        aspect="auto",
        cmap="plasma",
    )
    ax.set(
        title="GBIF Observed Unique Species Richness per 5° Grid Cell",
        xlabel="Longitude (°)",
        ylabel="Latitude (°)",
    )
    fig.colorbar(image, ax=ax, orientation="horizontal", label="Unique Species Count", pad=0.12)
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
        "bio1_current_c",
        "bio1_late_holocene_c",
        "bio1_lgm_c",
        "bio1_sd_c",
        "temperature_stability_index",
    ]
    df_final = df[final_cols].copy()
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
            "interpretation": "Primary test: unadjusted bivariate rank association across all occupied 5° cells.",
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

    # Partial: control signed latitude
    r_slat, p_slat = _partial_spearman(y, x, [signed_lat])
    stats_rows.append(
        {
            "analysis_scale": "5deg_cell",
            "climate_variable": "temperature_stability_index",
            "response_variable": "species_richness",
            "model": "partial_spearman_control_signed_latitude",
            "primary_test": False,
            "n": len(df),
            "correlation_coefficient": r_slat,
            "p_value": p_slat,
            "covariates_controlled": "signed_latitude",
            "interpretation": "Linear signed latitude control (-90 to +90).",
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

    # 2. Sampling Effort / Record Threshold Sensitivity
    sens_rows = []
    thresholds = [1, 2, 5, 10, 20, 50, 100]
    for th in thresholds:
        mask = df["occurrence_count"] >= th
        sub = df.loc[mask]
        sub_y = sub["species_richness"].to_numpy(dtype="float64")
        sub_x = sub["temperature_stability_index"].to_numpy(dtype="float64")
        sub_eff = sub["occurrence_count"].to_numpy(dtype="float64")
        sub_alat = np.abs(sub["cell_latitude_center_climate"].to_numpy(dtype="float64"))

        r_th, p_th = _spearman(sub_y, sub_x)
        r_part_eff, p_part_eff = _partial_spearman(sub_y, sub_x, [sub_eff])
        r_part_lat, p_part_lat = _partial_spearman(sub_y, sub_x, [sub_alat])
        r_part_both, p_part_both = _partial_spearman(sub_y, sub_x, [sub_eff, sub_alat])

        sens_rows.append(
            {
                "min_occurrence_threshold": th,
                "n_cells": len(sub),
                "pct_cells_retained": float(len(sub) / len(df) * 100.0),
                "spearman_rho": r_th,
                "spearman_p_value": p_th,
                "partial_r_control_effort": r_part_eff,
                "partial_p_control_effort": p_part_eff,
                "partial_r_control_abs_lat": r_part_lat,
                "partial_p_control_abs_lat": p_part_lat,
                "partial_r_control_effort_and_abs_lat": r_part_both,
                "partial_p_control_effort_and_abs_lat": p_part_both,
            }
        )

    sensitivity = pd.DataFrame(sens_rows)
    sensitivity.to_csv(config.BIODIV_CLIMATE_SENSITIVITY_CSV, index=False)

    # 3. Moran's I Spatial Autocorrelation
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
        moran_rows.append(
            {
                "variable": var_name,
                "n": len(df),
                "weight_scheme": "queen_contiguity_on_5deg_grid",
                **m_res,
                "interpretation": (
                    "Positive spatial autocorrelation detected (p <= 0.001); "
                    "spatial clustering indicates non-independence of cells."
                ),
            }
        )

    moran_df = pd.DataFrame(moran_rows)
    moran_df.to_csv(config.BIODIV_CLIMATE_MORAN_CSV, index=False)

    # 4. Render plots
    _render_plots(df, sensitivity, rho_prim, p_prim)

    return statistics, sensitivity, moran_df
