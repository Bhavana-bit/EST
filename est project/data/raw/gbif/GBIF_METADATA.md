# Local GBIF Data Notes

## Source Files

Phase 1 reads regional raw CSVs from `data/raw/gbif/`:

| File | Records | Sampling cap |
|---|---:|---:|
| `gbif_africa_5000.csv` | 5,000 | 5,000 |
| `gbif_asia_5000.csv` | 5,000 | 5,000 |
| `gbif_europe_5000.csv` | 5,000 | 5,000 |
| `gbif_north_america_1000.csv` | 1,000 | 1,000 |
| `gbif_south_america_1000.csv` | 1,000 | 1,000 |
| `gbif_oceania_1000.csv` | 1,000 | 1,000 |

These caps were applied during download and are **not equal across continents**. The combined raw table contains 18,000 occurrence records.

## Available Fields

Each raw CSV contains:

- `key` — GBIF occurrence record identifier
- `scientificName`
- `taxonRank`
- `decimalLatitude`, `decimalLongitude`
- `speciesKey`

The local files do **not** include occurrence year, basis of record, coordinate uncertainty, taxonomic class/order, or GBIF download DOI. Download query and filtering provenance are **not verified** from the local CSV alone.

## Dataset Character

The records are **mixed-taxonomic** and are **not restricted to birds or any single class**. The dataset is a **sample-based occurrence archive**, not a complete census of global biodiversity. Geographic and taxonomic reporting is uneven because of the unequal continental caps and GBIF reporting bias.

Removing duplicate occurrence keys or repeated records of the same species at identical coordinates reduces inflation in occurrence counts, but it **does not make the GBIF sample unbiased** or representative of true species richness.

## Reproducible Cleaning

Run:

```bash
python run_gbif_pipeline.py --clean-only
```

Cleaning stages (`src/gbif_clean.py`):

1. Load regional raw CSVs.
2. Keep records with valid WGS84 coordinates.
3. Keep species-level records (`taxonRank == SPECIES`) with a non-null `speciesKey`.
4. Remove exact duplicate occurrence records (duplicate `key` or fully identical rows).
5. Remove repeated records of the same `speciesKey` at identical coordinates, keeping the first record.

Outputs:

- `processed_data/gbif/gbif_clean.csv` — final cleaned occurrence records (with `continent`)
- `processed_data/gbif/gbif_species_occurrences.csv` — coordinate/species subset
- `processed_data/gbif/gbif_sampling_effort_5deg.csv` — per 5° cell `occurrence_count` (sampling effort) and `species_richness` (unique species)
- `data/processed/gbif_species_richness_5deg.tif` — 5° GeoTIFF raster of unique species richness
- `results/gbif/gbif_cleaning_summary.csv` — row counts removed at each stage
- `results/gbif/gbif_species_summary.csv` — final record and unique-species totals
- `results/gbif/gbif_continent_sampling_summary.csv` — per-continent sampling caps, records, species, and cell effort metrics

## Phase 1 Measures

Two spatial summaries are used:

- **5° sampling-effort grid** (`gbif_sampling_effort_5deg.csv`): `occurrence_count` (valid records per cell) and `species_richness` (unique `speciesKey` values per cell). This grid summarizes GBIF sampling effort at a coarse scale.
- **WorldClim 10 arc-minute analysis grid** (Phase 1 pipeline): observed species richness and occurrence count per occupied climate cell for correlation with thermal-change proxies. This keeps the biodiversity–climate analysis spatially aligned with the climate rasters.

Use the term **“GBIF-observed species richness across recorded taxa.”** This occurrence-based measure does not represent true or complete biodiversity.
