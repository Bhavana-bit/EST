# Palaeoclimate Stability and Modern Biodiversity

A scientifically rigorous, reproducible research framework investigating whether modern biodiversity-rich regions are disproportionately associated with long-term climatic stability.

---

## 1. Official Dataset Inventory

| Dataset | Official Source | Version / Product Used | Spatial Resolution | CRS | Temporal Period | Key Variables | Units | Primary Research Role |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| **WorldClim** | [worldclim.org](https://www.worldclim.org/) | WorldClim v1.4 | 10 arc-minutes (~18.5 km) | EPSG:4326 | Present (1960–1990), Mid-Holocene (~6ka CCSM4), LGM (~22ka CCSM4) | `BIO1` (Annual Mean Temp), `BIO12` (Annual Precip) | `BIO1`: °C × 10<br>`BIO12`: mm/yr | Primary palaeoclimate stability source ($\Delta T = \|T_{\text{present}} - T_{\text{LGM}}\|$) |
| **PaleoClim** | [paleoclim.org](http://www.paleoclim.org/) | PaleoClim v1.0 | 10 arc-minutes (or resampled 2.5') | EPSG:4326 | Mid-Holocene (~6ka), LGM (~21ka), LIG (~130ka) | `BIO1`, `BIO12`, `BIO4`, `BIO15` | °C, mm/yr | Independent palaeoclimate model validation source |
| **GBIF** | [gbif.org](https://www.gbif.org/) | API Occurrence Snapshot (Mammalia) | Georeferenced Point Records | EPSG:4326 | Modern Era | `species`, `decimalLatitude`, `decimalLongitude`, `occurrenceStatus`, `basisOfRecord` | Latitude/Longitude coordinates | Modern empirical species richness ($S = \text{unique species / cell}$) |
| **IUCN Red List** | [iucnredlist.org](https://www.iucnredlist.org/) | IUCN Spatial Range Polygons | Vector Polygons | EPSG:4326 | Modern Range Boundaries | Range Polygons (`presence==1`, `origin==1`) | Geographic Polygons | Independent expert-validated range richness cross-validation |
| **Copernicus DEM** | [copernicus-dem](https://registry.opendata.aws/copernicus-dem/) | GLO-90 / Aggregated 10' Grid | 90m aggregated to 10 arc-min | EPSG:4326 | Modern Relief | Elevation ($m$), Topographic Heterogeneity / Roughness ($\sigma_z$) | Meters ($m$) | Topographic relief & microrefugia confounding control |

---

## 2. Directory Structure

```
est project/
├── DATA_STATUS.md
├── README.md
├── requirements.txt
├── config.py
├── run_pipeline.py
├── data/
│   ├── raw/
│   └── processed/
├── src/
│   ├── download_all_data.py
│   ├── 01_validate_data.py
│   ├── 02_prepare_worldclim.py
│   ├── 03_prepare_paleoclim.py
│   ├── 04_prepare_gbif.py
│   ├── 05_prepare_iucn.py
│   ├── 06_prepare_dem.py
│   ├── 07_calculate_stability.py
│   ├── 08_biodiversity_analysis.py
│   └── 09_visualization.py
└── results/
    ├── maps/
    ├── tables/
    └── statistics/
```

---

## 3. How to Execute

### Option A: Complete Pipeline Run
```bash
py run_pipeline.py
```

### Option B: Automated Dataset Download
```bash
py src/download_all_data.py
```

---

## 4. Key Methodological Principles

1. **GBIF Species Richness Calculation**:
   - $S_i =$ Count of **UNIQUE** species in cell $i$, NOT raw record count.
   - Sampling effort ($\log(\text{Records}_i)$) included as a covariate in regression to control spatial sampling bias.

2. **WorldClim Temperature Unit Verification**:
   - WorldClim 1.4 `BIO1` values are encoded as °C × 10. They are divided by 10 to yield °C before computing temperature differences ($\Delta T$).

3. **Climate Stability Index ($I_{\text{stability}}$)**:
   - $I_{\text{stability}} = 1 - \frac{1}{2}(\widetilde{\Delta T_{\text{LGM}}} + \widetilde{\Delta P_{\text{LGM}}})$, where $\widetilde{\Delta T}$ and $\widetilde{\Delta P}$ are 99th-percentile normalized anomalies.

4. **Observational Language**:
   - All findings report statistical **associations** and **correlations**, avoiding unfounded causal claims.
