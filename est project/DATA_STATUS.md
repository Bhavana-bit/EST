# Dataset Status & Product Integrity

> [!CAUTION]
> **PRELIMINARY: mixed-taxon, capped sample (5,000/1,000 per continent), no GBIF download key. Not valid for taxon richness.**

## Dataset Readiness Overview

| Dataset | Local Status | Pipeline Status & Product Identity |
|---|---|---|
| **Paleoclimate Rasters** | **Present** | Rasters under `data/raw/paleoclim/` (`paleoclim_current_BIO1.tif`, `paleoclim_late_holocene_BIO1.tif`, `paleoclim_LGM_BIO1_aligned.tif`). Used for 5° LGM-dominated temperature stability index (Headline Analysis) and 10′ Phase 1 $\Delta T$ comparison. |
| **WorldClim** | **Present & Supporting** | Present, LGM, and mid-Holocene BIO1 rasters at 10 arc-minute resolution under `data/raw/worldclim/`. Used for Phase 1 10′ thermal change proxies. |
| **GBIF** | **Present & Cleaned** | Regional raw occurrence downloads ($16,521$ records cleaned reproducibly via `src/gbif_clean.py`). Aggregated to 5° grid cells for sampling effort (`occurrence_count`) and species richness (`species_richness` & rarefied richness). **PRELIMINARY: mixed-taxon, capped sample (5,000/1,000 per continent), no GBIF download key. Not valid for taxon richness.** |
| **IUCN** | **Missing** | `data/raw/iucn/` is empty; no local range-polygon data are present. Marked as future extension. |
| **Copernicus DEM** | **Missing** | `data/raw/dem/` is empty; no local elevation/DEM rasters are present. Marked as future extension. |

---

## Methodological Status & Team Decisions

- **Taxon Selection & Scope:** Confirmed by team (Aves, GBIF taxonKey 212).
- **Spatial Grid Scale:** Confirmed by team (5° × 5° grid).
- **Coordinate Uncertainty Filtering Threshold:** Decision required from project team.

---

## Analysis Identification

1. **Headline Analysis (5° Climate Stability Index):**
   - Evaluates multi-period temperature variability on a 5° global grid using BIO1 across Current, Late Holocene, and LGM.
   - Index formula: `1 / (1 + SD(Current, Late Holocene, LGM))`.
   - LGM-dominated: SD index is 0.9996 correlated with Present–LGM change.
   - **GBIF Results Notice:** PRELIMINARY: mixed-taxon, capped sample (5,000/1,000 per continent), no GBIF download key. Not valid for taxon richness.

2. **Phase 1 Supporting Analysis (10′ ΔT):**
   - Evaluates absolute temperature change at 10 arc-minute resolution for WorldClim Present–LGM, WorldClim Present–mid-Holocene, and CHELSA/PaleoClim Present–LGM.
   - **GBIF Results Notice:** PRELIMINARY: mixed-taxon, capped sample (5,000/1,000 per continent), no GBIF download key. Not valid for taxon richness.

---

## File Verification & Generated Artifacts

- All cleaned GBIF occurrences are written to `processed_data/gbif/gbif_clean.csv`.
- Climate stability rasters are saved in `processed_data/climate/climate_stability_5deg.tif` and `data/processed/climate/climate_stability_5deg.tif`.
- Final statistical tables and maps are exported to `results/`, `results/tables/`, and `results/maps/`.