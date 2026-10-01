# Phase 1 Data Status

## Five-Dataset Readiness

| Intended dataset | Local status | Pipeline status |
|---|---|---|
| PaleoClim | **Missing.** No verified PaleoClim product/archive is present. Files under `data/raw/paleoclim/` are CHELSA v1.2B and PaleoView products, based on their archive names and local XML lineage. | Genuine PaleoClim is not used. CHELSA current–LGM is analyzed separately and labelled CHELSA. |
| WorldClim | Present: extracted present BIO1 and LGM/mid-Holocene BIO1 rasters, with original source ZIP archives retained. | Used in Phase 1. |
| GBIF | Present: regional raw CSVs with unequal continental caps (5,000 or 1,000 records). Reproducible cleaning is implemented in `src/gbif_clean.py`. The cleaned records are mixed-taxonomic and provenance is incomplete. | Used in Phase 1 as a sample-based occurrence dataset, not a biodiversity census. |
| IUCN | **Missing.** `data/raw/iucn/` is empty; no local range-polygon data were found in either project tree. | Not used. |
| Copernicus DEM | **Missing.** `data/raw/dem/` is empty; no local DEM source raster/archive was found in either project tree. | Not used. |

**Overall five-dataset project status: NOT COMPLETE.** The current three-product phase is runnable with WorldClim, GBIF, and the locally available CHELSA current/LGM rasters. Genuine PaleoClim, IUCN, and Copernicus DEM inputs are absent. CHELSA/PaleoView are not substitutes for PaleoClim.

## Inputs Used

| Input | Local file | Use |
|---|---|---|
| WorldClim present BIO1 | `data/raw/worldclim/present/bio/bio_1/` | Authoritative analysis grid and present baseline |
| WorldClim LGM BIO1 | `data/raw/worldclim/lgm/cclgmbi1.tif` | Primary past climate comparison |
| WorldClim mid-Holocene BIO1 | `data/raw/worldclim/mid_holocene/ccmidbi1.tif` | Separate secondary comparison |
| GBIF cleaned observations | `processed_data/gbif/gbif_clean.csv` | Mixed-taxon sample-based occurrence records after reproducible cleaning |
| GBIF 5° sampling effort | `processed_data/gbif/gbif_sampling_effort_5deg.csv` | Per-cell occurrence count and species richness on a 5° grid |
| GBIF cleaning summary | `results/gbif/gbif_cleaning_summary.csv` | Row counts removed at each cleaning stage |
| CHELSA current BIO1 | `data/raw/paleoclim/paleoclim_current_BIO1.tif` | Independent current baseline for CHELSA comparison |
| CHELSA LGM BIO1 | `data/raw/paleoclim/paleoclim_LGM_BIO1.tif` | Independent CHELSA LGM comparison |

The present WorldClim BIO1 grid is EPSG:4326, 2160 columns by 900 rows, with nominal 10 arc-minute resolution and approximately −180° to 180° longitude and −60° to 90° latitude coverage. Its exact transform is read from the raster and used for every output. LGM and mid-Holocene rasters are geospatially reprojected to that grid. Source NoData is masked before BIO1 is divided by 10 to convert °C × 10 to °C.

## Definitions

- `DeltaT_LGM = abs(Present BIO1 - LGM BIO1)` in °C. Interpret only as a present–LGM thermal-change proxy for climatic stability.
- `DeltaT_MH = abs(Present BIO1 - Mid-Holocene BIO1)` in °C, reported as a separate secondary contrast.
- Per-cell observed richness is the number of unique `speciesKey` values in the cleaned GBIF records.
- Per-cell occurrence count is the number of GBIF records and is a sampling-effort indicator, not richness.
- The 5° grid reports coarse-scale sampling effort (`occurrence_count`, `species_richness`). The Phase 1 climate grid remains the WorldClim 10 arc-minute raster for spatial comparability with thermal-change proxies.

The cleaned GBIF table contains species-rank records but mixed taxa. Outputs must be called “GBIF-observed species richness across recorded taxa”; they are not bird richness, mammal richness, or a complete census of global biodiversity. Continental download caps, geographic/taxonomic reporting bias, and unequal geographic cell area limit interpretation. Removing duplicate occurrence keys or species-coordinate repeats reduces count inflation but does not make the sample unbiased. Correlations are associations, not causal evidence.

## Not Used

The CHELSA v1.2B current and LGM rasters are used only in a separate comparison; they are not combined with WorldClim inputs. The PaleoView late-Holocene raster is retained but not analyzed. No genuine PaleoClim, IUCN, or DEM analysis is implemented.

## Generated Outputs

Run `python run_phase1.py` to regenerate the WorldClim and CHELSA DeltaT rasters, GBIF per-cell tables, the WorldClim combined analysis table, separate CHELSA comparison tables, summary statistics, maps, and scatter plots. Run with `--validate-only` to check inputs without writing outputs. The `data/raw/iucn/` and `data/raw/dem/` directories are placeholders; genuine PaleoClim is also missing. Only regenerated `phase1_` outputs in `results/maps/` and the current CSVs in `results/tables/` are submission outputs.