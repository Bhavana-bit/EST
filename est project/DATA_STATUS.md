# DATA_STATUS.md — Phase 1 Dataset Verification Report

**Project Title**: 24. Palaeoclimate Stability and Modern Biodiversity  
**Phase**: Phase 1 Mid-Evaluation Pipeline  
**Date of Audit**: 2026-09-30  
**Target Analysis Grid**: WorldClim 10 arc-minutes (1080 rows x 2160 cols, EPSG:4326)

---

## 1. DATASET STATUS INVENTORY TABLE

| Dataset | Period / Model | Spatial Resolution | CRS | Variables Verified | Units & Scaling | File / Data Status | Verification Note |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| **WorldClim 1.4 Present** | 1960–1990 Baseline | 10 arc-minutes (~18.5 km) | EPSG:4326 | `BIO1` (Annual Mean Temp)<br>`BIO12` (Annual Precip) | `BIO1`: °C × 10 (Integer)<br>`BIO12`: mm/yr | Staged / Processing Ready | Divide `BIO1` by 10.0 to convert to °C |
| **WorldClim 1.4 LGM** | ~22,000 BP (CCSM4 GCM) | 10 arc-minutes (~18.5 km) | EPSG:4326 | `BIO1` (Annual Mean Temp)<br>`BIO12` (Annual Precip) | `BIO1`: °C × 10 (Integer)<br>`BIO12`: mm/yr | Staged / Processing Ready | Primary palaeoclimate contrast for Phase 1 |
| **WorldClim 1.4 Mid-Holocene** | ~6,000 BP (CCSM4 GCM) | 10 arc-minutes (~18.5 km) | EPSG:4326 | `BIO1` (Annual Mean Temp) | `BIO1`: °C × 10 (Integer) | Staged / Retained | Retained for Phase 2 extended analysis |
| **GBIF Occurrences** | Modern Era (1950–Present) | Georeferenced Point Records | EPSG:4326 | `speciesKey`, `species`, `decimalLatitude`, `decimalLongitude` | Lat/Lon Coordinates | Cleaned Snapshot Ready | Binned to 10' cells for UNIQUE species count ($S$) |
| **PaleoClim v1.0** | LGM (~21ka), Mid-Holocene (~6ka) | 10 arc-minutes | EPSG:4326 | `BIO1`, `BIO12` | °C, mm/yr | Phase 2 Deferred | Independent validation in Phase 2 |
| **IUCN Spatial Data** | Modern Range Boundaries | Vector Range Polygons | EPSG:4326 | Range Polygons (`presence==1`, `origin==1`) | Vector Geometry | Phase 2 Deferred | Independent range validation in Phase 2 |
| **Copernicus DEM** | Modern Topography | 90m / 10' Aggregated Grid | EPSG:4326 | Elevation ($m$), Roughness ($\sigma_z$) | Meters ($m$) | Phase 2 Deferred | Topographic heterogeneity control in Phase 2 |

---

## 2. PHASE 1 VERIFICATION & METHODOLOGICAL DECISIONS

1. **Common Spatial Grid**:
   - WorldClim 10 arc-minute global grid ($1080 \text{ rows} \times 2160 \text{ columns}$, extent $[-180, 180, -90, 90]$).
   - **Unequal-Area Limitation**: 10-arc-minute geographic grid cells decrease in physical surface area from the equator toward the poles. In Phase 1, raw cell species counts are evaluated with this limitation explicitly documented.

2. **Project-Defined Climate Stability Metric**:
   - Primary contrast: Present (1960–1990) vs. Last Glacial Maximum ($\sim 22,000\text{ BP}$ CCSM4 model).
   - Metric:
     $$\Delta T = | \text{Present BIO1 (°C)} - \text{LGM BIO1 (°C)} |$$
   - Smaller $\Delta T$ represents greater temperature stability under this project-defined metric.

3. **GBIF Species Richness Calculation**:
   - Taxon: **Mammalia** (Class Key 359).
   - Calculated as $S_i = \text{number of UNIQUE species recorded in cell } i$. Raw GBIF record count is NOT used as species richness.

4. **Phase 1 Output Deliverables**:
   - MAP 1: Present-vs-LGM BIO1 absolute temperature difference ($\Delta T$).
   - MAP 2: GBIF Mammalia unique species richness ($S$).
   - Scatter Plot: Log-transformed species richness $\log(S+1)$ vs. Temperature Difference $\Delta T$.
   - Descriptive Summary Table & Mid-Evaluation Documentation.
