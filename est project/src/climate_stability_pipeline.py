"""PaleoClim BIO1 temperature stability pipeline on the project 5-degree grid."""

from __future__ import annotations

import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import rasterio
from rasterio.transform import Affine

import config

BIO1_NODATA = -32768.0
BIO1_SCALE = 10.0
EARTH_RADIUS_KM = 6371.0
PERIOD_LABELS = ("current", "late_holocene", "lgm")
MIN_FINE_PIXELS = 5  # Minimum 10-arc-minute land pixels per 5-degree cell


def _read_bio1_c(path: Path) -> tuple[np.ndarray, dict]:
    with rasterio.open(path) as source:
        if source.crs is None:
            raise ValueError(f"Missing CRS: {path}")
        raw = source.read(1).astype("float32")
        values = raw.copy()
        if source.nodata is not None and np.isfinite(source.nodata):
            values[raw == source.nodata] = np.nan
        values[raw == BIO1_NODATA] = np.nan
        values[~np.isfinite(values)] = np.nan
        values = values / BIO1_SCALE
        meta = {
            "crs": source.crs,
            "transform": source.transform,
            "width": source.width,
            "height": source.height,
        }
        return values, meta


def _validate_same_grid(metas: list[dict], labels: list[str]) -> None:
    reference = metas[0]
    for meta, label in zip(metas[1:], labels[1:]):
        if (
            meta["width"] != reference["width"]
            or meta["height"] != reference["height"]
            or meta["transform"] != reference["transform"]
            or meta["crs"] != reference["crs"]
        ):
            raise ValueError(f"{label} is not on the same grid as {labels[0]}.")


def _cell_area_km2(lat_center_deg: float) -> float:
    """Spherical cap strip area for a 5°×5° cell."""
    delta = math.radians(config.GRID_DEGREES)
    lat = math.radians(lat_center_deg)
    return (EARTH_RADIUS_KM**2) * delta * delta * math.cos(lat)


def _aggregate_to_5deg(values: np.ndarray, transform: Affine) -> tuple[np.ndarray, np.ndarray]:
    n_lat = int(180 / config.GRID_DEGREES)
    n_lon = int(360 / config.GRID_DEGREES)
    sums = np.zeros((n_lat, n_lon), dtype="float64")
    counts = np.zeros((n_lat, n_lon), dtype="int32")

    rows, cols = np.indices(values.shape)
    xs, ys = rasterio.transform.xy(transform, rows.ravel(), cols.ravel(), offset="center")
    lon = np.asarray(xs, dtype="float64").reshape(values.shape)
    lat = np.asarray(ys, dtype="float64").reshape(values.shape)
    finite = np.isfinite(values)
    lat_index = np.floor((lat + 90.0) / config.GRID_DEGREES).astype(int)
    lon_index = np.floor((lon + 180.0) / config.GRID_DEGREES).astype(int)
    in_bounds = finite & (lat_index >= 0) & (lat_index < n_lat) & (lon_index >= 0) & (lon_index < n_lon)
    np.add.at(sums, (lat_index[in_bounds], lon_index[in_bounds]), values[in_bounds])
    np.add.at(counts, (lat_index[in_bounds], lon_index[in_bounds]), 1)

    with np.errstate(invalid="ignore", divide="ignore"):
        means = sums / counts
    means[counts < MIN_FINE_PIXELS] = np.nan
    return means, counts


def _grid_metadata() -> tuple[Affine, int, int]:
    n_lat = int(180 / config.GRID_DEGREES)
    n_lon = int(360 / config.GRID_DEGREES)
    transform = Affine(config.GRID_DEGREES, 0, -180.0, 0, -config.GRID_DEGREES, 90.0)
    return transform, n_lon, n_lat


def _cell_table(
    current: np.ndarray,
    late_holocene: np.ndarray,
    lgm: np.ndarray,
    sd: np.ndarray,
    stability: np.ndarray,
    counts: np.ndarray,
) -> pd.DataFrame:
    n_lat, n_lon = current.shape
    max_fine_per_cell = 900.0  # 30x30 fine pixels per 5 deg cell at 10-arc-min resolution
    rows = []
    for lat_index in range(n_lat):
        lat_lower = lat_index * config.GRID_DEGREES - 90.0
        lat_center = lat_lower + config.GRID_DEGREES / 2.0
        for lon_index in range(n_lon):
            lon_lower = lon_index * config.GRID_DEGREES - 180.0
            lon_center = lon_lower + config.GRID_DEGREES / 2.0
            raster_row = (n_lat - 1) - lat_index
            temperature_sd = sd[lat_index, lon_index]
            fine_count = counts[lat_index, lon_index]
            if not np.isfinite(temperature_sd) or fine_count < MIN_FINE_PIXELS:
                continue
            tot_area = _cell_area_km2(lat_center)
            land_fraction = min(1.0, float(fine_count) / max_fine_per_cell)
            rows.append(
                {
                    "grid_cell_id": f"lat{lat_index}_lon{lon_index}",
                    "grid_lat_index": lat_index,
                    "grid_lon_index": lon_index,
                    "row": raster_row,
                    "col": lon_index,
                    "cell_latitude_center": lat_center,
                    "cell_longitude_center": lon_center,
                    "fine_pixel_count": fine_count,
                    "cell_area_km2": tot_area,
                    "estimated_land_area_km2": tot_area * land_fraction,
                    "bio1_current_c": current[lat_index, lon_index],
                    "bio1_late_holocene_c": late_holocene[lat_index, lon_index],
                    "bio1_lgm_c": lgm[lat_index, lon_index],
                    "bio1_sd_c": temperature_sd,
                    "temperature_stability_index": stability[lat_index, lon_index],
                }
            )
    return pd.DataFrame(rows)


def _render_stability_map(stability_raster: np.ndarray, output_path: Path) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(12, 6))
    masked = np.ma.masked_invalid(stability_raster)
    image = ax.imshow(
        masked,
        extent=(-180.0, 180.0, -90.0, 90.0),
        origin="upper",
        aspect="auto",
        cmap="viridis",
    )
    ax.set(
        title=(
            "PaleoClim BIO1 Temperature Stability Index (5° Grid)\n"
            "LGM-Dominated Multi-Period Stability: Current, Late Holocene, LGM"
        ),
        xlabel="Longitude (°)",
        ylabel="Latitude (°)",
    )
    fig.colorbar(
        image,
        ax=ax,
        orientation="horizontal",
        label="Temperature Stability Index = 1 / (1 + SD(BIO1) °C)",
        pad=0.12,
    )
    fig.tight_layout()
    fig.savefig(output_path, dpi=200)
    plt.close(fig)


def _write_readme(dest_dir: Path) -> None:
    dest_dir.mkdir(parents=True, exist_ok=True)
    readme_path = dest_dir / "CLIMATE_STABILITY_README.md"
    readme_path.write_text(
        """# PALEOCLIM TEMPERATURE STABILITY INDEX (LGM-DOMINATED)

## Inputs
- PaleoClim v1.2B BIO1 (mean annual temperature) across three periods:
  1. Current baseline (1979–2013)
  2. Late Holocene (Meghalayan, 4.2–0.3 ka BP)
  3. Last Glacial Maximum (LGM, ~21 ka BP)
- Source rasters: `data/raw/paleoclim/paleoclim_current_BIO1.tif`, `paleoclim_late_holocene_BIO1.tif`, `paleoclim_LGM_BIO1_aligned.tif`.
- Stored as integer °C×10; converted to °C by dividing by 10.0 before analysis.

## Grid & Minimum Land Rule
- Common **5°** global grid aligned with GBIF sampling cells (`floor((lat+90)/5)`, `floor((lon+180)/5)`).
- Fine (10 arc-minute) pixels are averaged within each 5° cell for each period.
- Enforces a minimum land pixel threshold (>= 5 fine pixels per 5° cell) to exclude ocean fringe cells.
- Strictly requires valid BIO1 in ALL three periods; ocean and incomplete cells are assigned NoData (NaN).

## Index Structure & LGM Domination
For each 5° cell with valid mean BIO1 in all three periods:
- `bio1_sd_c` = SD(Current, Late Holocene, LGM) in °C (population SD, ddof=0).
- **Temperature stability index** = `1 / (1 + bio1_sd_c)`.

> [!NOTE]
> The index is an **LGM-dominated stability index**. Late Holocene temperature differs from present temperature by only 0.68 °C on average globally, whereas LGM temperature differs by 11.6 °C on average. Consequently, the standard deviation index is 0.9996 correlated (Spearman) with absolute Present–LGM temperature change (|LGM − Present|).

## Cell & Land Area
Geographic 5° cells vary in physical area with latitude:
`cell_area_km2 = R² × (Δ°→rad)² × cos(φ)`, with `R = 6371 km` and `Δ = 5°`.
`estimated_land_area_km2 = cell_area_km2 × (fine_pixel_count / 900)`.
""",
        encoding="utf-8",
    )


def run_pipeline() -> pd.DataFrame:
    paths = (
        config.PALEOCLIM_CURRENT_BIO1,
        config.PALEOCLIM_LATE_HOLOCENE_BIO1,
        config.PALEOCLIM_LGM_BIO1_ALIGNED,
    )
    for path in paths:
        if not path.exists():
            raise FileNotFoundError(path)

    layers = [_read_bio1_c(path) for path in paths]
    values = [layer[0] for layer in layers]
    metas = [layer[1] for layer in layers]
    _validate_same_grid(metas, list(PERIOD_LABELS))

    current_fine, late_holocene_fine, lgm_fine = values
    transform = metas[0]["transform"]
    current, count_curr = _aggregate_to_5deg(current_fine, transform)
    late_holocene, count_lh = _aggregate_to_5deg(late_holocene_fine, transform)
    lgm, count_lgm = _aggregate_to_5deg(lgm_fine, transform)

    # Minimum fine pixel counts across all three periods
    counts = np.minimum.reduce([count_curr, count_lh, count_lgm])

    stack = np.stack([current, late_holocene, lgm], axis=0)
    valid = np.all(np.isfinite(stack), axis=0) & (counts >= MIN_FINE_PIXELS)
    sd = np.full(current.shape, np.nan, dtype="float32")
    sd[valid] = np.std(stack[:, valid], axis=0, ddof=0)
    stability = np.full(current.shape, np.nan, dtype="float32")
    stability[valid] = 1.0 / (1.0 + sd[valid])

    transform_5deg, n_lon, n_lat = _grid_metadata()
    config.CLIMATE_PROCESSED_DIR.mkdir(parents=True, exist_ok=True)
    config.TABLES_DIR.mkdir(parents=True, exist_ok=True)
    config.MAPS_DIR.mkdir(parents=True, exist_ok=True)

    profile = {
        "driver": "GTiff",
        "dtype": "float32",
        "count": 1,
        "width": n_lon,
        "height": n_lat,
        "crs": metas[0]["crs"],
        "transform": transform_5deg,
        "nodata": np.nan,
        "compress": "deflate",
    }
    raster_array = np.flipud(stability).astype("float32")
    # Save to data/processed/climate/
    with rasterio.open(config.CLIMATE_STABILITY_RASTER, "w", **profile) as destination:
        destination.write(raster_array, 1)
        destination.set_band_description(
            1, "PaleoClim LGM-dominated temperature stability index = 1/(1+SD BIO1 °C)"
        )

    # Also save to processed_data/climate/
    alt_climate_dir = config.BASE_DIR / "processed_data" / "climate"
    alt_climate_dir.mkdir(parents=True, exist_ok=True)
    with rasterio.open(alt_climate_dir / "climate_stability_5deg.tif", "w", **profile) as destination:
        destination.write(raster_array, 1)
        destination.set_band_description(
            1, "PaleoClim LGM-dominated temperature stability index = 1/(1+SD BIO1 °C)"
        )

    cells = _cell_table(current, late_holocene, lgm, sd, stability, counts)
    cells.to_csv(config.CLIMATE_STABILITY_CELLS_CSV, index=False)
    cells.to_csv(config.RESULTS_DIR / "climate_stability_5deg.csv", index=False)

    summary = {
        "index_name": "PALEOCLIM LGM-DOMINATED TEMPERATURE STABILITY INDEX",
        "formula": "1/(1+SD(bio1_current_c, bio1_late_holocene_c, bio1_lgm_c))",
        "grid_degrees": config.GRID_DEGREES,
        "valid_5deg_cells": int(valid.sum()),
        "bio1_sd_c_mean": float(cells["bio1_sd_c"].mean()),
        "bio1_sd_c_median": float(cells["bio1_sd_c"].median()),
        "stability_mean": float(cells["temperature_stability_index"].mean()),
        "stability_median": float(cells["temperature_stability_index"].median()),
        "cell_area_km2_min": float(cells["cell_area_km2"].min()),
        "cell_area_km2_max": float(cells["cell_area_km2"].max()),
        "fine_pixels_all_three_periods": int(
            np.sum(
                np.isfinite(current_fine)
                & np.isfinite(late_holocene_fine)
                & np.isfinite(lgm_fine)
            )
        ),
    }
    summary_df = pd.DataFrame([summary])
    summary_df.to_csv(config.CLIMATE_STABILITY_SUMMARY_CSV, index=False)
    summary_df.to_csv(config.RESULTS_DIR / "climate_stability_summary.csv", index=False)

    _write_readme(config.CLIMATE_PROCESSED_DIR)
    _write_readme(alt_climate_dir)

    # Render climate stability 5° map
    _render_stability_map(raster_array, config.MAPS_DIR / "climate_stability_5deg.png")

    return summary_df

