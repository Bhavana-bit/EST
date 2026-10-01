"""Reproducible GBIF occurrence cleaning and 5-degree sampling-effort aggregation."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
import rasterio
from rasterio.transform import Affine

import config

GBIF_RAW_COLUMNS = {
    "key",
    "scientificName",
    "taxonRank",
    "decimalLatitude",
    "decimalLongitude",
    "speciesKey",
}
REGIONAL_RAW_FILES = (
    "gbif_africa_5000.csv",
    "gbif_asia_5000.csv",
    "gbif_europe_5000.csv",
    "gbif_north_america_1000.csv",
    "gbif_south_america_1000.csv",
    "gbif_oceania_1000.csv",
)
FILENAME_TO_CONTINENT = {
    "gbif_africa_5000.csv": "Africa",
    "gbif_asia_5000.csv": "Asia",
    "gbif_europe_5000.csv": "Europe",
    "gbif_north_america_1000.csv": "North America",
    "gbif_south_america_1000.csv": "South America",
    "gbif_oceania_1000.csv": "Oceania",
}
GRID_DEGREES = 5.0


def _verify_gbif_metadata_provenance() -> None:
    metadata_path = config.GBIF_RAW_DIR / "GBIF_METADATA.md"
    if not metadata_path.exists():
        raise FileNotFoundError(
            f"GBIF cleaning stopped: Missing metadata file {metadata_path}. "
            "Follow MANUAL_DATA_ACQUISITION.md to record official GBIF download provenance."
        )
    content = metadata_path.read_text(encoding="utf-8")
    content_lower = content.lower()

    if "mixed-taxonomic" in content_lower:
        raise ValueError(
            "GBIF cleaning stopped: GBIF_METADATA.md contains 'mixed-taxonomic'. "
            "An official GBIF download restricted to Aves (taxonKey 212) is required per MANUAL_DATA_ACQUISITION.md."
        )

    required_fields = {
        "Download Key:": r"Download\s+Key\s*:\s*(.*)",
        "DOI:": r"DOI\s*:\s*(.*)",
        "Download Date:": r"Download\s+Date\s*:\s*(.*)",
        "Query Filters:": r"Query\s+Filters\s*:\s*(.*)",
        "Record Count:": r"Record\s+Count\s*:\s*(.*)",
    }

    placeholders = {"none", "tbd", "xxx", "empty", "null", "no download key"}
    missing_or_invalid = []
    parsed_values = {}

    import re

    for label, pattern in required_fields.items():
        match = re.search(pattern, content, re.IGNORECASE)
        if not match:
            missing_or_invalid.append(f"Missing label '{label}'")
            continue
        val = match.group(1).strip().strip("*`_# ").strip()
        val_lower = val.lower()

        if not val or val_lower in placeholders or any(val_lower.startswith(p) for p in placeholders):
            missing_or_invalid.append(f"Placeholder/empty value for '{label}' (got '{val}')")
        else:
            parsed_values[label] = val

    if "Query Filters:" in parsed_values:
        qf_clean = parsed_values["Query Filters:"].lower().replace(" ", "").replace(":", "=")
        if "taxonkey=212" not in qf_clean:
            missing_or_invalid.append(
                f"Query Filters: must include 'taxonKey=212' (got '{parsed_values['Query Filters:']}')"
            )

    if "Record Count:" in parsed_values:
        rc_clean = parsed_values["Record Count:"].replace(",", "").replace(".", "").strip()
        if not rc_clean.isdigit():
            missing_or_invalid.append(
                f"Record Count: must be an integer (got '{parsed_values['Record Count:']}')"
            )

    if missing_or_invalid:
        raise ValueError(
            "GBIF cleaning stopped: GBIF_METADATA.md has invalid or missing provenance fields:\n - "
            + "\n - ".join(missing_or_invalid)
            + "\nFollow MANUAL_DATA_ACQUISITION.md to record official GBIF download metadata for Aves (taxonKey 212)."
        )


def _load_raw_gbif() -> pd.DataFrame:
    # Stop if the six capped regional CSVs are still present in data/raw/gbif/
    found_capped = [f for f in REGIONAL_RAW_FILES if (config.GBIF_RAW_DIR / f).exists()]
    if found_capped:
        raise ValueError(
            f"GBIF cleaning stopped: Capped regional sample CSV files found in {config.GBIF_RAW_DIR}: {found_capped}. "
            "An official uncapped GBIF download for Aves (taxonKey 212) is required. "
            "Please follow MANUAL_DATA_ACQUISITION.md."
        )

    _verify_gbif_metadata_provenance()

    raw_files = list(config.GBIF_RAW_DIR.glob("*.csv")) + list(config.GBIF_RAW_DIR.glob("*.txt"))
    raw_files = [f for f in raw_files if f.name != "GBIF_METADATA.md"]
    if not raw_files:
        raise FileNotFoundError(
            f"No raw GBIF download file found in {config.GBIF_RAW_DIR}. "
            "Follow MANUAL_DATA_ACQUISITION.md to place your official GBIF Aves download."
        )

    frames = []
    for path in raw_files:
        frame = pd.read_csv(path, sep=None, engine="python", dtype={"speciesKey": "string"}, low_memory=False)
        missing = GBIF_RAW_COLUMNS.difference(frame.columns)
        if missing:
            raise ValueError(f"{path.name} is missing required columns: {sorted(missing)}")
        frames.append(frame)
    combined = pd.concat(frames, ignore_index=True)

    # ClassKey column presence check
    if "classKey" not in combined.columns:
        raise ValueError("GBIF dataset missing required 'classKey' column.")

    non_null_class = combined["classKey"].dropna().astype(str).str.strip()
    if len(non_null_class) == 0:
        raise ValueError("GBIF dataset contains no non-null 'classKey' values.")

    aves_matches = (non_null_class == str(config.TARGET_TAXON_KEY)) | (non_null_class == "212")
    share = float(aves_matches.sum()) / float(len(non_null_class))
    print(f"GBIF Dataset Validation: Found {share:.2%} of non-null classKey rows matching Aves (taxonKey 212).")

    if share < 0.99:
        raise ValueError(
            f"GBIF validation failed: Found {share:.2%} of non-null classKey rows matching Aves (taxonKey {config.TARGET_TAXON_KEY}), "
            f"which is below the required 99.00% threshold."
        )

    return combined


def _coordinates_valid(data: pd.DataFrame) -> pd.Series:
    latitude = pd.to_numeric(data["decimalLatitude"], errors="coerce")
    longitude = pd.to_numeric(data["decimalLongitude"], errors="coerce")
    return (
        np.isfinite(latitude)
        & np.isfinite(longitude)
        & latitude.between(-90, 90)
        & longitude.between(-180, 180)
    )


def _species_valid(data: pd.DataFrame) -> pd.Series:
    species_key = data["speciesKey"].astype("string").str.strip()
    return data["taxonRank"].eq("SPECIES") & species_key.notna() & species_key.ne("")


def _assign_5deg_cell(latitude: pd.Series, longitude: pd.Series) -> pd.DataFrame:
    n_lat = int(180 / GRID_DEGREES)
    n_lon = int(360 / GRID_DEGREES)
    lat_index = np.floor((latitude + 90.0) / GRID_DEGREES).astype(int)
    lon_index = np.floor((longitude + 180.0) / GRID_DEGREES).astype(int)
    lat_index = np.clip(lat_index, 0, n_lat - 1)
    lon_index = np.clip(lon_index, 0, n_lon - 1)
    lat_lower = lat_index * GRID_DEGREES - 90.0
    lon_lower = lon_index * GRID_DEGREES - 180.0
    return pd.DataFrame(
        {
            "grid_lat_index": lat_index,
            "grid_lon_index": lon_index,
            "grid_cell_id": [
                f"lat{lat_idx}_lon{lon_idx}" for lat_idx, lon_idx in zip(lat_index, lon_index)
            ],
            "cell_latitude_center": lat_lower + GRID_DEGREES / 2.0,
            "cell_longitude_center": lon_lower + GRID_DEGREES / 2.0,
        }
    )


def summarize_continent_sampling(raw: pd.DataFrame, cleaned: pd.DataFrame) -> pd.DataFrame:
    """Summarize regional sampling effort to document inequality and sampling caps."""
    cell_cols = _assign_5deg_cell(cleaned["decimalLatitude"], cleaned["decimalLongitude"])
    annotated = pd.concat([cleaned.reset_index(drop=True), cell_cols], axis=1)

    rows = []
    for filename in REGIONAL_RAW_FILES:
        continent = FILENAME_TO_CONTINENT[filename]
        raw_count = int((raw["source_file"] == filename).sum())
        clean_sub = cleaned.loc[cleaned["source_file"] == filename]
        clean_count = len(clean_sub)
        sp_count = int(clean_sub["speciesKey"].nunique())

        ann_sub = annotated.loc[annotated["source_file"] == filename]
        grouped = ann_sub.groupby("grid_cell_id").agg(
            records=("speciesKey", "size"),
            species=("speciesKey", "nunique"),
        )
        n_cells = len(grouped)
        med_rec = float(grouped["records"].median()) if n_cells > 0 else 0.0
        max_rec = int(grouped["records"].max()) if n_cells > 0 else 0
        med_sp = float(grouped["species"].median()) if n_cells > 0 else 0.0
        max_sp = int(grouped["species"].max()) if n_cells > 0 else 0

        rows.append(
            {
                "continent": continent,
                "raw_file": filename,
                "raw_records_downloaded": raw_count,
                "cleaned_records": clean_count,
                "unique_species": sp_count,
                "occupied_5deg_cells": n_cells,
                "median_records_per_cell": med_rec,
                "max_records_per_cell": max_rec,
                "median_species_per_cell": med_sp,
                "max_species_per_cell": max_sp,
                "sample_nature_note": (
                    f"Sample-based occurrence archive (capped at {raw_count} records during download); "
                    "not a complete biodiversity census."
                ),
            }
        )
    return pd.DataFrame(rows)


def export_5deg_richness_raster(grouped: pd.DataFrame, output_path: Path) -> None:
    """Export 5-degree unique species richness grid as GeoTIFF matching climate grid."""
    n_lat = int(180 / GRID_DEGREES)
    n_lon = int(360 / GRID_DEGREES)
    grid = np.full((n_lat, n_lon), np.nan, dtype="float32")
    for _, row in grouped.iterrows():
        lat_idx = int(row["grid_lat_index"])
        lon_idx = int(row["grid_lon_index"])
        raster_row = (n_lat - 1) - lat_idx
        grid[raster_row, lon_idx] = float(row["species_richness"])

    transform_5deg = Affine(GRID_DEGREES, 0, -180.0, 0, -GRID_DEGREES, 90.0)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    profile = {
        "driver": "GTiff",
        "dtype": "float32",
        "count": 1,
        "width": n_lon,
        "height": n_lat,
        "crs": "EPSG:4326",
        "transform": transform_5deg,
        "nodata": np.nan,
        "compress": "deflate",
    }
    with rasterio.open(output_path, "w", **profile) as dst:
        dst.write(grid, 1)
        dst.set_band_description(1, "GBIF unique species richness per 5-degree cell")


def clean_gbif() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Clean raw GBIF records and write processed outputs plus cleaning and sampling summaries."""
    config.GBIF_PROCESSED_DIR.mkdir(parents=True, exist_ok=True)
    config.GBIF_RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    config.TABLES_DIR.mkdir(parents=True, exist_ok=True)

    raw = _load_raw_gbif()
    raw_records = len(raw)

    coordinate_mask = _coordinates_valid(raw)
    after_coordinates = raw.loc[coordinate_mask].copy()
    records_after_coordinate_filtering = len(after_coordinates)

    species_mask = _species_valid(after_coordinates)
    after_species = after_coordinates.loc[species_mask].copy()
    records_after_species_filtering = len(after_species)

    duplicate_key_mask = after_species["key"].duplicated(keep="first")
    duplicate_row_mask = after_species.duplicated(keep="first")
    exact_duplicate_mask = duplicate_key_mask | duplicate_row_mask
    after_exact = after_species.loc[~exact_duplicate_mask].copy()
    exact_duplicates_removed = int(exact_duplicate_mask.sum())

    species_coordinate_mask = after_exact.duplicated(
        subset=["speciesKey", "decimalLatitude", "decimalLongitude"],
        keep="first",
    )
    cleaned = after_exact.loc[~species_coordinate_mask].copy()
    species_coordinate_duplicates_removed = int(species_coordinate_mask.sum())

    cleaned["decimalLatitude"] = pd.to_numeric(cleaned["decimalLatitude"], errors="coerce")
    cleaned["decimalLongitude"] = pd.to_numeric(cleaned["decimalLongitude"], errors="coerce")
    cleaned = cleaned.sort_values(["speciesKey", "decimalLatitude", "decimalLongitude", "key"]).reset_index(
        drop=True
    )

    final_records = len(cleaned)
    unique_species = int(cleaned["speciesKey"].nunique())

    output_columns = [
        "key",
        "scientificName",
        "taxonRank",
        "decimalLatitude",
        "decimalLongitude",
        "speciesKey",
        "continent",
        "source_file",
    ]
    # Save to active processed directory
    cleaned[output_columns].to_csv(config.GBIF_PROCESSED_DIR / "gbif_clean.csv", index=False)
    cleaned[["speciesKey", "scientificName", "continent", "decimalLatitude", "decimalLongitude"]].to_csv(
        config.GBIF_PROCESSED_DIR / "gbif_species_occurrences.csv",
        index=False,
    )

    # Sync to legacy data/processed/gbif directory to keep identical state
    legacy_gbif_dir = config.PROCESSED_DATA_DIR / "gbif"
    legacy_gbif_dir.mkdir(parents=True, exist_ok=True)
    cleaned[output_columns].to_csv(legacy_gbif_dir / "gbif_clean.csv", index=False)
    cleaned[["speciesKey", "scientificName", "continent", "decimalLatitude", "decimalLongitude"]].to_csv(
        legacy_gbif_dir / "gbif_species_occurrences.csv",
        index=False,
    )

    alt_legacy_dir = config.BASE_DIR.parent / "Palaeoclimate_Biodiversity_Project" / "processed_data" / "gbif"
    if alt_legacy_dir.parent.parent.exists():
        alt_legacy_dir.mkdir(parents=True, exist_ok=True)
        cleaned[output_columns].to_csv(alt_legacy_dir / "gbif_clean.csv", index=False)
        cleaned[["speciesKey", "scientificName", "continent", "decimalLatitude", "decimalLongitude"]].to_csv(
            alt_legacy_dir / "gbif_species_occurrences.csv",
            index=False,
        )

    summary = pd.DataFrame(
        [
            {
                "raw_records": raw_records,
                "records_after_coordinate_filtering": records_after_coordinate_filtering,
                "records_after_species_filtering": records_after_species_filtering,
                "exact_duplicates_removed": exact_duplicates_removed,
                "species_coordinate_duplicates_removed": species_coordinate_duplicates_removed,
                "final_records": final_records,
                "unique_species": unique_species,
            }
        ]
    )
    summary.to_csv(config.GBIF_RESULTS_DIR / "gbif_cleaning_summary.csv", index=False)

    species_summary = pd.DataFrame([{"unique_species": unique_species, "final_records": final_records}])
    species_summary.to_csv(config.GBIF_RESULTS_DIR / "gbif_species_summary.csv", index=False)

    sampling_effort = aggregate_sampling_effort_5deg(cleaned)

    # Continent sampling effort and inequality summary
    continent_summary = summarize_continent_sampling(raw, cleaned)
    continent_summary.to_csv(config.GBIF_CONTINENT_SUMMARY_CSV, index=False)
    continent_summary.to_csv(config.TABLES_DIR / "gbif_continent_sampling_summary.csv", index=False)

    # Export 5° richness raster matching climate grid
    export_5deg_richness_raster(sampling_effort, config.GBIF_RICHNESS_5DEG_RASTER)
    export_5deg_richness_raster(sampling_effort, legacy_gbif_dir / "gbif_species_richness_5deg.tif")

    return cleaned, summary, continent_summary


def _rarefy(species_counts: np.ndarray, target_k: int) -> float:
    N = int(species_counts.sum())
    if N < target_k:
        return np.nan
    probs = []
    k_arr = np.arange(target_k)
    for ns in species_counts:
        if N - ns < target_k:
            probs.append(0.0)
        else:
            p = np.exp(np.sum(np.log(N - ns - k_arr) - np.log(N - k_arr)))
            probs.append(p)
    return float(np.sum(1.0 - np.array(probs)))


def aggregate_sampling_effort_5deg(cleaned: pd.DataFrame | None = None) -> pd.DataFrame:
    """Aggregate occurrence count (sampling effort), unique species richness, and rarefied richness per 5-degree cell."""
    config.GBIF_PROCESSED_DIR.mkdir(parents=True, exist_ok=True)
    config.GBIF_RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    config.TABLES_DIR.mkdir(parents=True, exist_ok=True)

    if cleaned is None:
        cleaned = pd.read_csv(
            config.GBIF_PROCESSED_DIR / "gbif_clean.csv",
            dtype={"speciesKey": "string"},
            low_memory=False,
        )

    cell_columns = _assign_5deg_cell(cleaned["decimalLatitude"], cleaned["decimalLongitude"])
    annotated = pd.concat([cleaned.reset_index(drop=True), cell_columns], axis=1)

    # Basic counts
    grouped = (
        annotated.groupby(
            ["grid_cell_id", "grid_lat_index", "grid_lon_index", "cell_latitude_center", "cell_longitude_center"],
            sort=True,
        )
        .agg(
            occurrence_count=("speciesKey", "size"),
            species_richness=("speciesKey", "nunique"),
        )
        .reset_index()
    )

    # Rarefied species richness for k=5 and k=10
    rar_5_dict = {}
    rar_10_dict = {}
    for cell_id, group in annotated.groupby("grid_cell_id"):
        sp_counts = group["speciesKey"].value_counts().to_numpy(dtype=int)
        rar_5_dict[cell_id] = _rarefy(sp_counts, 5)
        rar_10_dict[cell_id] = _rarefy(sp_counts, 10)

    grouped["richness_rarefied_5"] = grouped["grid_cell_id"].map(rar_5_dict)
    grouped["richness_rarefied_10"] = grouped["grid_cell_id"].map(rar_10_dict)

    # Save to processed directory and results directories
    grouped.to_csv(config.GBIF_SAMPLING_EFFORT_5DEG, index=False)
    grouped.to_csv(config.GBIF_RESULTS_DIR / "gbif_sampling_effort_5deg.csv", index=False)
    grouped.to_csv(config.TABLES_DIR / "gbif_sampling_effort_5deg.csv", index=False)
    grouped.to_csv(config.RESULTS_DIR / "gbif_species_richness_5deg.csv", index=False)
    legacy_gbif_dir = config.PROCESSED_DATA_DIR / "gbif"
    legacy_gbif_dir.mkdir(parents=True, exist_ok=True)
    grouped.to_csv(legacy_gbif_dir / "gbif_sampling_effort_5deg.csv", index=False)
    return grouped


def run_gbif_pipeline() -> dict[str, float | int]:
    """Run GBIF cleaning and 5-degree sampling-effort aggregation."""
    cleaned, summary, continent_summary = clean_gbif()
    sampling_effort = aggregate_sampling_effort_5deg(cleaned)

    duplicates_removed = int(
        summary.loc[0, "exact_duplicates_removed"]
        + summary.loc[0, "species_coordinate_duplicates_removed"]
    )
    return {
        "final_record_count": int(summary.loc[0, "final_records"]),
        "final_unique_species_count": int(summary.loc[0, "unique_species"]),
        "occupied_5deg_cells": len(sampling_effort),
        "median_records_per_5deg_cell": float(sampling_effort["occurrence_count"].median()),
        "median_species_richness_per_5deg_cell": float(sampling_effort["species_richness"].median()),
        "duplicate_records_removed": duplicates_removed,
    }


def print_pipeline_report(metrics: dict[str, float | int]) -> None:
    print("GBIF cleaning pipeline complete.")
    print(f"Final record count: {metrics['final_record_count']}")
    print(f"Final unique species count: {metrics['final_unique_species_count']}")
    print(f"Occupied 5-degree cells: {metrics['occupied_5deg_cells']}")
    print(f"Median records per 5-degree cell: {metrics['median_records_per_5deg_cell']:.1f}")
    print(
        "Median species richness per 5-degree cell: "
        f"{metrics['median_species_richness_per_5deg_cell']:.1f}"
    )
    print(f"Duplicate records removed: {metrics['duplicate_records_removed']}")
