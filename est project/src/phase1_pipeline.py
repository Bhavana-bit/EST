"""Real-data Phase 1 workflow for WorldClim and GBIF."""

from __future__ import annotations

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import rasterio
from rasterio.crs import CRS
from rasterio.enums import Resampling
from rasterio.transform import xy
from rasterio.warp import reproject, transform as warp_transform
from scipy import stats

import config

BIO1_SCALE = 10.0
GBIF_REQUIRED_COLUMNS = {
    "key",
    "speciesKey",
    "taxonRank",
    "decimalLatitude",
    "decimalLongitude",
}


def _check_file(path, label):
    if not path.exists():
        raise FileNotFoundError(f"Missing {label}: {path}")


def _read_grid(path, label):
    with rasterio.open(path) as source:
        if source.crs is None:
            raise ValueError(f"{label} has no CRS: {path}")
        if source.nodata is None:
            raise ValueError(f"{label} has no declared NoData value: {path}")
        if source.count < 1:
            raise ValueError(f"{label} contains no raster bands: {path}")

        raw = source.read(1, masked=True).astype("float32")
        values = np.asarray(raw.filled(np.nan), dtype="float32") / BIO1_SCALE
        metadata = source.profile.copy()
        metadata.update(transform=source.transform, crs=source.crs)
        return values, metadata


def _align_to_target(source_values, source_meta, target_meta):
    aligned = np.full(
        (target_meta["height"], target_meta["width"]), np.nan, dtype="float32"
    )
    reproject(
        source=source_values,
        destination=aligned,
        src_transform=source_meta["transform"],
        src_crs=source_meta["crs"],
        src_nodata=np.nan,
        dst_transform=target_meta["transform"],
        dst_crs=target_meta["crs"],
        dst_nodata=np.nan,
        resampling=Resampling.bilinear,
        init_dest_nodata=True,
    )
    return aligned


def _validate_gbif():
    _check_file(config.GBIF_CSV, "clean GBIF CSV")
    data = pd.read_csv(config.GBIF_CSV, dtype={"speciesKey": "string"}, low_memory=False)
    missing = GBIF_REQUIRED_COLUMNS.difference(data.columns)
    if missing:
        raise ValueError(f"GBIF CSV is missing required columns: {sorted(missing)}")
    if not data["taxonRank"].eq("SPECIES").all():
        raise ValueError("The selected clean GBIF table contains non-SPECIES ranks.")
    if data["key"].isna().any():
        raise ValueError("The selected GBIF table contains records without occurrence keys.")

    duplicate_occurrences = data["key"].duplicated(keep="first")
    duplicate_species_coordinates = data.duplicated(
        subset=["speciesKey", "decimalLatitude", "decimalLongitude"],
        keep="first",
    )
    data.attrs["duplicate_occurrence_keys_removed"] = int(duplicate_occurrences.sum())
    data.attrs["duplicate_species_coordinate_records_removed"] = int(
        duplicate_species_coordinates.sum()
    )
    if duplicate_occurrences.any() or duplicate_species_coordinates.any():
        raise ValueError(
            "The selected GBIF table still contains duplicate occurrence keys or "
            "species-coordinate duplicates. Regenerate it with run_gbif_pipeline.py."
        )

    data["decimalLatitude"] = pd.to_numeric(data["decimalLatitude"], errors="coerce")
    data["decimalLongitude"] = pd.to_numeric(data["decimalLongitude"], errors="coerce")
    return data


def validate_inputs():
    """Check required local files and source metadata without writing outputs."""
    for path, label in (
        (config.WORLDCLIM_PRESENT_BIO1, "WorldClim present BIO1 grid"),
        (config.WORLDCLIM_LGM_BIO1, "WorldClim LGM BIO1 raster"),
        (config.WORLDCLIM_MID_HOLOCENE_BIO1, "WorldClim mid-Holocene BIO1 raster"),
        (config.CHELSA_CURRENT_BIO1, "CHELSA current BIO1 raster"),
        (config.CHELSA_LGM_BIO1, "CHELSA LGM BIO1 raster"),
        (config.GBIF_CSV, "clean GBIF CSV"),
    ):
        _check_file(path, label)

    grids = []
    for path, label in (
        (config.WORLDCLIM_PRESENT_BIO1, "WorldClim present BIO1"),
        (config.WORLDCLIM_LGM_BIO1, "WorldClim LGM BIO1"),
        (config.WORLDCLIM_MID_HOLOCENE_BIO1, "WorldClim mid-Holocene BIO1"),
        (config.CHELSA_CURRENT_BIO1, "CHELSA current BIO1"),
        (config.CHELSA_LGM_BIO1, "CHELSA LGM BIO1"),
    ):
        _, metadata = _read_grid(path, label)
        grids.append((label, metadata))

    target = grids[0][1]
    if CRS.from_user_input(target["crs"]) != CRS.from_epsg(4326):
        raise ValueError(f"Expected GBIF decimal coordinates in EPSG:4326; target is {target['crs']}.")
    print("Input raster checks:")
    for label, metadata in grids:
        print(
            f"  {label}: {metadata['width']} x {metadata['height']}, "
            f"CRS={metadata['crs']}, transform={metadata['transform']}"
        )

    data = _validate_gbif()
    valid_coordinates = (
        np.isfinite(data["decimalLatitude"])
        & np.isfinite(data["decimalLongitude"])
        & data["decimalLatitude"].between(-90, 90)
        & data["decimalLongitude"].between(-180, 180)
    )
    valid_keys = data["speciesKey"].notna() & data["speciesKey"].str.strip().ne("")
    print(
        f"GBIF checks: {len(data)} rows, {int((valid_coordinates & valid_keys).sum())} "
        "with valid WGS84 coordinates and speciesKey."
    )
    print("No output files created during validation.")


def _save_delta_raster(values, target_meta, path):
    profile = target_meta.copy()
    profile.pop("blockxsize", None)
    profile.pop("blockysize", None)
    profile.update(
        driver="GTiff",
        count=1,
        dtype="float32",
        nodata=np.nan,
        compress="deflate",
        predictor=3,
        tiled=False,
    )
    with rasterio.open(path, "w", **profile) as destination:
        destination.write(values.astype("float32"), 1)
        destination.set_band_description(1, "Absolute BIO1 temperature change (degrees C)")


def _mapped_occurrences(data, target_meta, climate_valid):
    coordinates_valid = (
        np.isfinite(data["decimalLatitude"])
        & np.isfinite(data["decimalLongitude"])
        & data["decimalLatitude"].between(-90, 90)
        & data["decimalLongitude"].between(-180, 180)
    )
    key_valid = data["speciesKey"].notna() & data["speciesKey"].str.strip().ne("")
    usable = data.loc[coordinates_valid & key_valid].copy()
    rejected = {
        "input_rows": len(data),
        "duplicate_occurrence_keys_removed": int(
            data.attrs.get("duplicate_occurrence_keys_removed", 0)
        ),
        "rejected_invalid_coordinates": int((~coordinates_valid).sum()),
        "rejected_missing_species_key": int((~key_valid).sum()),
    }

    source_crs = CRS.from_epsg(4326)
    target_crs = CRS.from_user_input(target_meta["crs"])
    if source_crs != target_crs:
        xs, ys = warp_transform(
            source_crs,
            target_crs,
            usable["decimalLongitude"].to_numpy().tolist(),
            usable["decimalLatitude"].to_numpy().tolist(),
        )
        xs = np.asarray(xs)
        ys = np.asarray(ys)
    else:
        xs = usable["decimalLongitude"].to_numpy()
        ys = usable["decimalLatitude"].to_numpy()

    rows, cols = rasterio.transform.rowcol(target_meta["transform"], xs, ys)
    usable["row"] = np.asarray(rows, dtype=np.int64)
    usable["col"] = np.asarray(cols, dtype=np.int64)
    inside = (
        usable["row"].between(0, target_meta["height"] - 1)
        & usable["col"].between(0, target_meta["width"] - 1)
    )
    rejected["rejected_outside_target_grid"] = int((~inside).sum())
    usable = usable.loc[inside].copy()

    climate_at_points = climate_valid[
        usable["row"].to_numpy(), usable["col"].to_numpy()
    ]
    rejected["rejected_climate_nodata_cells"] = int((~climate_at_points).sum())
    usable = usable.loc[climate_at_points].copy()
    rejected["records_used"] = len(usable)
    return usable, rejected


def _render_map(values, target_meta, title, path, label, use_log=False):
    bounds = rasterio.transform.array_bounds(
        target_meta["height"], target_meta["width"], target_meta["transform"]
    )
    masked = np.ma.masked_invalid(values)
    fig, ax = plt.subplots(figsize=(12, 6))
    norm = None
    if use_log:
        vmax = float(np.nanmax(values)) if np.any(np.isfinite(values)) else 100.0
        norm = matplotlib.colors.LogNorm(vmin=1.0, vmax=vmax)
    image = ax.imshow(
        masked,
        extent=(bounds[0], bounds[2], bounds[1], bounds[3]),
        origin="upper",
        aspect="auto",
        cmap="plasma" if use_log else "magma",
        norm=norm,
    )
    ax.set(title=title, xlabel="Longitude", ylabel="Latitude")
    fig.colorbar(image, ax=ax, orientation="horizontal", label=label, pad=0.12)
    fig.tight_layout()
    fig.savefig(path, dpi=200)
    plt.close(fig)


def _statistics(analysis, records_used, target_cells):
    result = {
        "gbif_input_records_used": records_used,
        "valid_analysis_cells": len(analysis),
        "occupied_target_cells": target_cells,
    }
    for name in (
        "delta_t_lgm_c",
        "delta_t_mh_c",
        "observed_species_richness",
        "gbif_occurrence_count",
    ):
        values = analysis[name].to_numpy(dtype=float)
        result[f"{name}_mean"] = float(np.mean(values))
        result[f"{name}_median"] = float(np.median(values))
        result[f"{name}_minimum"] = float(np.min(values))
        result[f"{name}_maximum"] = float(np.max(values))

    for climate_name, climate_column in (
        ("lgm", "delta_t_lgm_c"),
        ("mid_holocene", "delta_t_mh_c"),
    ):
        pearson = stats.pearsonr(analysis[climate_column], analysis["observed_species_richness"])
        spearman = stats.spearmanr(analysis[climate_column], analysis["observed_species_richness"])
        result[f"pearson_r_{climate_name}"] = float(pearson.statistic)
        result[f"pearson_p_{climate_name}"] = float(pearson.pvalue)
        result[f"spearman_rho_{climate_name}"] = float(spearman.statistic)
        result[f"spearman_p_{climate_name}"] = float(spearman.pvalue)
    return pd.DataFrame([result])


def run_pipeline():
    validate_inputs()
    config.PROCESSED_DATA_DIR.mkdir(parents=True, exist_ok=True)
    config.MAPS_DIR.mkdir(parents=True, exist_ok=True)
    config.TABLES_DIR.mkdir(parents=True, exist_ok=True)
    config.RESULTS_DIR.mkdir(parents=True, exist_ok=True)

    present, target_meta = _read_grid(config.WORLDCLIM_PRESENT_BIO1, "WorldClim present BIO1")
    lgm_raw, lgm_meta = _read_grid(config.WORLDCLIM_LGM_BIO1, "WorldClim LGM BIO1")
    mh_raw, mh_meta = _read_grid(config.WORLDCLIM_MID_HOLOCENE_BIO1, "WorldClim mid-Holocene BIO1")
    lgm = _align_to_target(lgm_raw, lgm_meta, target_meta)
    mid_holocene = _align_to_target(mh_raw, mh_meta, target_meta)

    paleoclim_current_raw, paleoclim_current_meta = _read_grid(
        config.PALEOCLIM_CURRENT_BIO1, "PaleoClim current BIO1"
    )
    paleoclim_lgm_raw, paleoclim_lgm_meta = _read_grid(
        config.PALEOCLIM_LGM_BIO1, "PaleoClim LGM BIO1"
    )
    paleoclim_current = _align_to_target(paleoclim_current_raw, paleoclim_current_meta, target_meta)
    paleoclim_lgm = _align_to_target(paleoclim_lgm_raw, paleoclim_lgm_meta, target_meta)
    delta_paleoclim_lgm = np.abs(paleoclim_current - paleoclim_lgm)
    paleoclim_valid = np.isfinite(paleoclim_current) & np.isfinite(paleoclim_lgm)
    delta_paleoclim_lgm[~paleoclim_valid] = np.nan

    delta_lgm = np.abs(present - lgm)
    delta_mh = np.abs(present - mid_holocene)
    climate_valid = np.isfinite(present) & np.isfinite(delta_lgm) & np.isfinite(delta_mh)
    delta_lgm[~climate_valid] = np.nan
    delta_mh[~climate_valid] = np.nan

    _save_delta_raster(delta_lgm, target_meta, config.PROCESSED_DATA_DIR / "delta_t_lgm.tif")
    _save_delta_raster(delta_mh, target_meta, config.PROCESSED_DATA_DIR / "delta_t_mid_holocene.tif")
    _save_delta_raster(
        delta_paleoclim_lgm,
        target_meta,
        config.PROCESSED_DATA_DIR / "delta_t_paleoclim_current_lgm.tif",
    )

    occurrences, filtering = _mapped_occurrences(_validate_gbif(), target_meta, climate_valid)
    grouped = (
        occurrences.groupby(["row", "col"], sort=True)
        .agg(
            observed_species_richness=("speciesKey", "nunique"),
            gbif_occurrence_count=("speciesKey", "size"),
        )
        .reset_index()
    )
    rows = grouped["row"].to_numpy(dtype=int)
    cols = grouped["col"].to_numpy(dtype=int)
    center_x, center_y = xy(target_meta["transform"], rows, cols, offset="center")
    grouped.insert(0, "grid_cell_id", [f"{row}_{col}" for row, col in zip(rows, cols)])
    grouped["longitude"] = center_x
    grouped["latitude"] = center_y
    grouped["delta_t_lgm_c"] = delta_lgm[rows, cols]
    grouped["delta_t_mh_c"] = delta_mh[rows, cols]

    final_columns = [
        "grid_cell_id", "row", "col", "longitude", "latitude",
        "delta_t_lgm_c", "delta_t_mh_c", "observed_species_richness",
        "gbif_occurrence_count",
    ]
    analysis = grouped[final_columns].copy()
    analysis[["grid_cell_id", "row", "col", "longitude", "latitude", "observed_species_richness"]].to_csv(
        config.TABLES_DIR / "gbif_species_richness_per_cell.csv", index=False
    )
    analysis[["grid_cell_id", "row", "col", "longitude", "latitude", "gbif_occurrence_count"]].to_csv(
        config.TABLES_DIR / "gbif_occurrence_count_per_cell.csv", index=False
    )
    analysis.to_csv(config.TABLES_DIR / "final_analysis.csv", index=False)
    analysis.to_csv(config.RESULTS_DIR / "final_analysis.csv", index=False)

    summary = _statistics(analysis, filtering["records_used"], len(grouped))
    for name in (
        "duplicate_occurrence_keys_removed",
        "rejected_invalid_coordinates",
        "rejected_missing_species_key",
        "rejected_outside_target_grid",
        "rejected_climate_nodata_cells",
    ):
        summary[name] = filtering[name]
    summary.to_csv(config.TABLES_DIR / "summary_statistics.csv", index=False)
    summary.to_csv(config.RESULTS_DIR / "summary_statistics.csv", index=False)

    paleoclim_mask = paleoclim_valid[rows, cols]
    paleoclim_cells = grouped.loc[paleoclim_mask, [
        "grid_cell_id", "row", "col", "longitude", "latitude",
        "observed_species_richness", "gbif_occurrence_count",
    ]].copy()
    paleoclim_rows = paleoclim_cells["row"].to_numpy(dtype=int)
    paleoclim_cols = paleoclim_cells["col"].to_numpy(dtype=int)
    paleoclim_cells["delta_t_paleoclim_lgm_c"] = delta_paleoclim_lgm[paleoclim_rows, paleoclim_cols]
    paleoclim_cells.to_csv(config.TABLES_DIR / "paleoclim_lgm_cell_analysis.csv", index=False)
    paleoclim_cells.to_csv(config.TABLES_DIR / "chelsa_lgm_cell_analysis.csv", index=False)

    paleoclim_summary = {
        "product": "PaleoClim v1.2B current versus LGM BIO1",
        "gbif_records_used": int(paleoclim_cells["gbif_occurrence_count"].sum()),
        "valid_occupied_cells": len(paleoclim_cells),
        "delta_t_paleoclim_lgm_c_mean": float(paleoclim_cells["delta_t_paleoclim_lgm_c"].mean()),
        "delta_t_paleoclim_lgm_c_median": float(paleoclim_cells["delta_t_paleoclim_lgm_c"].median()),
        "delta_t_paleoclim_lgm_c_minimum": float(paleoclim_cells["delta_t_paleoclim_lgm_c"].min()),
        "delta_t_paleoclim_lgm_c_maximum": float(paleoclim_cells["delta_t_paleoclim_lgm_c"].max()),
        "observed_species_richness_mean": float(paleoclim_cells["observed_species_richness"].mean()),
        "gbif_occurrence_count_mean": float(paleoclim_cells["gbif_occurrence_count"].mean()),
    }
    paleoclim_pearson = stats.pearsonr(
        paleoclim_cells["delta_t_paleoclim_lgm_c"], paleoclim_cells["observed_species_richness"]
    )
    paleoclim_spearman = stats.spearmanr(
        paleoclim_cells["delta_t_paleoclim_lgm_c"], paleoclim_cells["observed_species_richness"]
    )
    paleoclim_summary["pearson_r"] = float(paleoclim_pearson.statistic)
    paleoclim_summary["pearson_p"] = float(paleoclim_pearson.pvalue)
    paleoclim_summary["spearman_rho"] = float(paleoclim_spearman.statistic)
    paleoclim_summary["spearman_p"] = float(paleoclim_spearman.pvalue)
    pd.DataFrame([paleoclim_summary]).to_csv(
        config.TABLES_DIR / "paleoclim_lgm_summary.csv", index=False
    )
    pd.DataFrame([paleoclim_summary]).to_csv(
        config.TABLES_DIR / "chelsa_lgm_summary.csv", index=False
    )

    richness_grid = np.full(delta_lgm.shape, np.nan, dtype="float32")
    richness_grid[rows, cols] = grouped["observed_species_richness"].to_numpy(dtype="float32")
    _render_map(
        delta_lgm, target_meta, "WorldClim Present–LGM Thermal Change Proxy",
        config.MAPS_DIR / "phase1_delta_t_lgm.png", "Absolute BIO1 change (°C)"
    )
    _render_map(
        delta_mh, target_meta, "WorldClim Present–Mid-Holocene Thermal Change",
        config.MAPS_DIR / "phase1_delta_t_mid_holocene.png", "Absolute BIO1 change (°C)"
    )
    _render_map(
        delta_paleoclim_lgm,
        target_meta,
        "PaleoClim v1.2B Current–LGM Thermal Change Proxy",
        config.MAPS_DIR / "phase1_paleoclim_current_lgm_delta_t.png",
        "Absolute BIO1 change (°C)",
    )
    _render_map(
        richness_grid, target_meta, "GBIF-Observed Species Richness Across Recorded Taxa (10' Grid)",
        config.MAPS_DIR / "phase1_gbif_observed_richness.png", "Unique speciesKey values per cell (log scale)",
        use_log=True,
    )

    for contrast, climate_column in (
        ("lgm", "delta_t_lgm_c"),
        ("mid_holocene", "delta_t_mh_c"),
        ("paleoclim_current_lgm", "delta_t_paleoclim_lgm_c"),
    ):
        plot_data = paleoclim_cells if contrast == "paleoclim_current_lgm" else analysis
        fig, ax = plt.subplots(figsize=(8, 6))
        ax.scatter(plot_data[climate_column], plot_data["observed_species_richness"], alpha=0.35, s=12)
        ax.set_yscale("log")
        ax.set(
            title=(
                "GBIF-observed richness vs PaleoClim current-LGM thermal change"
                if contrast == "paleoclim_current_lgm"
                else f"GBIF-observed richness vs present-{contrast.replace('_', ' ')} thermal change"
            ),
            xlabel="Absolute BIO1 change (°C)",
            ylabel="GBIF-observed species richness across recorded taxa (log scale)",
        )
        fig.tight_layout()
        output_name = (
            "phase1_richness_vs_paleoclim_current_lgm.png"
            if contrast == "paleoclim_current_lgm"
            else f"phase1_richness_vs_delta_t_{contrast}.png"
        )
        fig.savefig(config.MAPS_DIR / output_name, dpi=200)
        if contrast == "paleoclim_current_lgm":
            fig.savefig(config.MAPS_DIR / "phase1_richness_vs_chelsa_current_lgm.png", dpi=200)
        plt.close(fig)

    print(f"GBIF records used in climate-valid cells: {filtering['records_used']}")
    print(f"Occupied cells with both climate contrasts valid: {len(analysis)}")
    print(
        f"Pearson r: LGM={summary.loc[0, 'pearson_r_lgm']:.6f}; "
        f"mid-Holocene={summary.loc[0, 'pearson_r_mid_holocene']:.6f}"
    )
    print(
        f"Spearman rho: LGM={summary.loc[0, 'spearman_rho_lgm']:.6f}; "
        f"mid-Holocene={summary.loc[0, 'spearman_rho_mid_holocene']:.6f}"
    )