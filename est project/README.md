# Palaeoclimate Stability and Modern Biodiversity

**Research question:** Are present biodiversity-rich regions disproportionately associated with long-term climatic stability?

> [!CAUTION]
> **PRELIMINARY: mixed-taxon, capped sample (5,000/1,000 per continent), no GBIF download key. Not valid for taxon richness.**

---

## Executive Summary & Analysis Structure

This project investigates the relationship between multi-period paleoclimatic temperature stability and global patterns of modern biodiversity richness using reproducibly cleaned GBIF species occurrence data and high-resolution paleoclimate rasters.

The repository contains two complementary analyses:

1. **Headline Analysis (5° Climate Stability Index):**
   - **Scale:** 5° × 5° global grid (aligned with GBIF sampling effort).
   - **Climate Metric:** BIO1 (mean annual temperature) across three periods (Current, Late Holocene, and Last Glacial Maximum).
   - **Stability Index:** LGM-dominated temperature stability index defined as `1 / (1 + SD(Current, Late Holocene, LGM))`.
   - **Findings (PRELIMINARY: mixed-taxon, capped sample (5,000/1,000 per continent), no GBIF download key. Not valid for taxon richness):**
     Primary result: unadjusted bivariate correlation across all occupied 5° cells is $\rho = +0.026$ ($n = 369, p = 0.623$). Controlling for GBIF sampling effort yields a partial correlation of $r = +0.231$ ($p < 0.0001$), and fixed-sample rarefaction (10 records/cell) yields $r = +0.391$ ($p < 0.00001$).

2. **Phase 1 Supporting Analysis (10′ ΔT):**
   - **Scale:** Fine 10 arc-minute geographic grid.
   - **Climate Metric:** Absolute temperature change ($\Delta T = |\text{Present} - \text{Past}|$) for WorldClim LGM, WorldClim mid-Holocene, and CHELSA/PaleoClim LGM.
   - **Findings (PRELIMINARY: mixed-taxon, capped sample (5,000/1,000 per continent), no GBIF download key. Not valid for taxon richness):** Provides high-resolution baseline spatial proxies for local thermal displacement.

---

## Methodological Status & Team Decisions

- **Taxon Selection & Scope:** Confirmed by team (Aves, GBIF taxonKey 212).
- **Spatial Grid Scale:** Confirmed by team (5° × 5° grid).
- **Coordinate Uncertainty Filtering Threshold:** Decision required from project team.

---

## Input Data & Product Verification

- **Paleoclimatic Temperature Sources:** Rasters under `data/raw/paleoclim/` and `data/raw/worldclim/`. Stored as integer °C × 10 and converted to °C before analysis.
- **GBIF Occurrences (PRELIMINARY: mixed-taxon, capped sample (5,000/1,000 per continent), no GBIF download key. Not valid for taxon richness):** Reproducibly cleaned global occurrence database (`gbif_clean.csv`, $n = 16,521$ records across 5,895 unique species). Aggregated to 5° grid cells for sampling effort (`occurrence_count`) and species richness (`species_richness`).
- **Missing Intended Datasets:** IUCN range polygons and Copernicus DEM data are not present in this release and remain marked as future extensions.

---

## Methodological & Statistical Notes

> [!IMPORTANT]
> **1. LGM-Dominated Stability Index:**
> The multi-period standard deviation index is effectively an **LGM-dominated index**. Late Holocene temperature differs from present temperature by an average of only 0.68 °C, whereas LGM temperature differs by 11.6 °C. The SD index has a Spearman correlation of $r = 0.9996$ with absolute Present–LGM change ($|\text{LGM} - \text{Present}|$).
>
> **2. Spatial Autocorrelation & Unadjusted p-values:**
> Strong positive spatial autocorrelation is present in species richness (Moran's $I = 0.284, p \le 0.001$), climate stability ($I = 0.793, p \le 0.001$), and model regression residuals ($I = 0.233, p \le 0.001$). Standard p-values ignore spatial clustering and should be treated as unadjusted/overstated.
>
> **3. Sampling Effort Control & Rarefaction (PRELIMINARY: mixed-taxon, capped sample (5,000/1,000 per continent), no GBIF download key. Not valid for taxon richness):**
> Unadjusted richness reflects sampling density. Primary result: $\rho = +0.026$ ($p = 0.623$). Partial correlation controlling for effort yields $r = +0.231$, and fixed-sample rarefied richness (10 records/cell) yields $r = +0.391$.
>
> **4. Regional Subsets (PRELIMINARY: mixed-taxon, capped sample (5,000/1,000 per continent), no GBIF download key. Not valid for taxon richness):**
> Extratropical cells ($|\text{latitude}| > 23.5^\circ$) show $\rho = +0.218$ ($p = 0.0015$), whereas Tropical cells ($|\text{latitude}| \le 23.5^\circ$) show $\rho = -0.190$ ($p = 0.0159$).
>
> **5. Sampling Threshold Sensitivity (PRELIMINARY: mixed-taxon, capped sample (5,000/1,000 per continent), no GBIF download key. Not valid for taxon richness):**
> As minimum cell occurrence threshold increases from $\ge 1$ to $\ge 100$ records, $\rho$ increases from $+0.026$ ($n=369$) to $+0.530$ ($n=35$). At $\ge 100$ records, 30 of 35 cells are in Asia, Europe, and Africa.
>
> **6. Geographic & Taxonomic Sample Bias (PRELIMINARY: mixed-taxon, capped sample (5,000/1,000 per continent), no GBIF download key. Not valid for taxon richness):**
> Sampling is taxonomically and geographically uneven. For example, Europe has 4,795 records covering only 404 species, whereas Africa has 4,308 records covering 2,238 species.

---

## Execution Order & Reproduction

To reproduce all clean data tables, climate stability rasters, statistical models, and figures, run the runner scripts in the following order:

```bash
# Step 1: Clean raw GBIF records and aggregate 5° sampling effort
py run_gbif_pipeline.py --clean-only

# Step 2: Build 5° temperature stability raster & cell table
py run_climate_stability_pipeline.py

# Step 3: Run headline 5° statistical analysis, partial correlations, Moran's I, & threshold CIs
py run_biodiversity_climate_stats.py

# Step 4: Run Phase 1 10′ ΔT supporting analysis
py run_phase1.py
```

---

## Shipped Outputs Directory Structure

All generated outputs are saved cleanly in the following locations:

- `results/`
  - `final_biodiversity_climate_5deg.csv` — Primary merged 5° analysis dataset (richness, effort, rarefied richness, stability index).
  - `biodiversity_climate_statistics.csv` — Statistical models (Spearman, partial correlations, regional subsets, rarefaction).
  - `biodiversity_climate_sensitivity_thresholds.csv` — Record count threshold sensitivity with 95% CIs.
  - `biodiversity_climate_moran_i.csv` — Spatial autocorrelation Moran's I tests and dynamic interpretations.
  - `climate_stability_5deg.csv` & `climate_stability_summary.csv` — 5° climate stability cell data and summary.
  - `final_analysis.csv` & `summary_statistics.csv` — Phase 1 10′ cell statistics and contrasts.
- `results/tables/` — Identical tabular output CSVs formatted for documentation & publication.
- `results/maps/`
  - `biodiversity_vs_climate_stability_5deg.png` — Main headline scatter plot (log scale richness, dynamic stats text box).
  - `biodiversity_vs_stability_by_threshold.png` — Sensitivity plot across record thresholds with 95% CIs.
  - `climate_stability_5deg.png` — Global 5° climate stability index map.
  - `gbif_species_richness_5deg.png` — Global 5° GBIF species richness map (log scale).
  - `phase1_*.png` — Phase 1 10′ thermal change maps and scatter plots.
- `processed_data/`
  - `climate/climate_stability_5deg.tif` — Project 5° GeoTIFF raster for climate stability.
  - `gbif/gbif_clean.csv` — Reproducibly cleaned GBIF occurrence dataset.