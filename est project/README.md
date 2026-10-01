# Palaeoclimate Stability and Modern Biodiversity

**Research question:** Are present biodiversity-rich regions disproportionately associated with long-term climatic stability?

## Current Work

This mid-evaluation phase uses three locally available datasets/products:

- **WorldClim:** present, LGM, and mid-Holocene BIO1. BIO1 is converted from °C × 10 to °C. Present–LGM absolute change is the primary thermal-change proxy; present–mid-Holocene change is reported separately.
- **GBIF:** species-rank occurrence records from a reproducibly cleaned local table (`run_gbif_pipeline.py`). The raw download is capped unequally by continent (5,000 or 1,000 records). Cleaning keeps valid species-level coordinates, removes exact duplicate occurrence records, and removes repeated records of the same species at identical coordinates. Richness is unique `speciesKey` count per occupied grid cell; occurrence count is retained as a sampling-effort indicator. The records are mixed-taxonomic sample data, not a complete census of global biodiversity.
- **CHELSA v1.2B:** a separate current–LGM BIO1 comparison using the locally available rasters. These files are not PaleoClim. The PaleoView late-Holocene archive is retained but not analyzed.

The WorldClim present BIO1 raster defines the target grid. Past climate rasters and CHELSA rasters are aligned by geospatial reprojection using their CRS and transforms. The WorldClim and CHELSA comparisons are calculated separately. Associations are descriptive and do not establish causation.

**Future additions:** IUCN range data and Copernicus DEM data are not implemented. A verified PaleoClim product is also not present, so the complete five-dataset project is not finished.

## Run

Regenerate cleaned GBIF tables and the 5° sampling-effort summary:

```bash
python run_gbif_pipeline.py --clean-only
```

Run the full GBIF cleaning plus Phase 1 climate-richness analysis:

```bash
python run_gbif_pipeline.py
```

Validate Phase 1 inputs without writing outputs:

```bash
python run_phase1.py --validate-only
```

Current tables and figures are saved in `results/tables/` and `results/maps/`. Cleaned GBIF tables are stored in `processed_data/gbif/`; cleaning summaries are in `results/gbif/`. Climate DeltaT rasters are stored in `data/processed/`. The pipeline uses local data only and has no synthetic-data fallback.

## Main Outputs

- `results/tables/final_analysis.csv` and `summary_statistics.csv`: WorldClim contrasts joined with GBIF cell richness and record counts.
- `results/tables/chelsa_lgm_cell_analysis.csv` and `chelsa_lgm_summary.csv`: separate CHELSA current–LGM comparison.
- `results/maps/`: current thermal-change maps, GBIF richness map, and scatter plots.

GBIF occurrences are unevenly sampled geographically and taxonomically. The 10-arc-minute geographic grid has unequal cell area by latitude. Both limit interpretation of observed cell richness.