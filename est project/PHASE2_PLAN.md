# PHASE 2 IMPLEMENTATION PLAN & ROADMAP

**Project Title**: 24. Palaeoclimate Stability and Modern Biodiversity  
**Status**: Deferred to Post-Mid-Evaluation Phase

---

## 1. OVERVIEW OF PHASE 2 TASKS

Phase 2 builds upon the working Phase 1 pipeline by adding multi-source validation, topographic controls, precipitation dynamics, spatial autocorrelation modeling, and sensitivity testing.

```
+-----------------------------------------------------------------------------+
|                               PHASE 2 ROADMAP                               |
+-------------------+-------------------+-------------------+-----------------+
|     PHASE 2A      |     PHASE 2B      |     PHASE 2C      |    PHASE 2D     |
| PaleoClim         | IUCN Range        | Copernicus DEM    | Precipitation   |
| Comparison        | Validation        | Controls          | Stability (BIO12|
+-------------------+-------------------+-------------------+-----------------+
|     PHASE 2E      |     PHASE 2F      |     PHASE 2G      |                 |
| Spatial Models    | Sensitivity &     | Multi-Taxa        |                 |
| (SAR / Moran's I) | Equal-Area Grids  | Breakdown         |                 |
+-------------------+-------------------+-------------------+-----------------+
```

---

## 2. DETAILED PHASE 2 MODULE SPECIFICATIONS

### Phase 2A: PaleoClim Independent Comparison
* **Objective**: Evaluate whether palaeoclimate stability patterns computed from WorldClim 1.4 CCSM4 are consistent with independent **PaleoClim v1.0** reconstructions.
* **Datasets**: PaleoClim LGM (~21 ka), Mid-Holocene (~6 ka), and Last Interglacial (~130 ka).
* **Protocol**: Calculate $\Delta T_{\text{PaleoClim}} = | \text{BIO1}_{\text{Present}} - \text{BIO1}_{\text{PaleoClim, LGM}} |$ separately on the 10 arc-minute grid. Cross-correlate $\Delta T_{\text{WorldClim}}$ with $\Delta T_{\text{PaleoClim}}$ without mixing raw raster values into a single metric.

### Phase 2B: IUCN Range-Based Biodiversity Validation
* **Objective**: Validate GBIF empirical occurrence richness ($S_{\text{GBIF}}$) against expert-validated **IUCN Red List** species distribution range polygons.
* **Filter Rules**: `presence == 1` (Extant), `origin == 1` (Native), `seasonal IN (1, 2)` (Resident/Breeding).
* **Protocol**: Vector-to-raster polygon overlay onto the 10 arc-minute grid to compute $S_{\text{IUCN}}$. Assess spatial agreement via grid-cell correlation $r(S_{\text{GBIF}}, S_{\text{IUCN}})$.

### Phase 2C: Copernicus DEM Topographic Heterogeneity Controls
* **Objective**: Disentangle macroclimate stability from topographic microrefugia effects.
* **Dataset**: Copernicus DEM GLO-90 / GLO-30.
* **Metrics**: For each 10' cell, compute:
  1. Mean Elevation ($\mu_z$ in meters)
  2. Elevation Standard Deviation / Terrain Roughness ($\sigma_z$ in meters)
  3. Elevation Range ($\max(z) - \min(z)$ in meters)
* **Statistical Control**: Include $\sigma_z$ as a covariate in multivariate regression models.

### Phase 2D: Annual Precipitation Stability ($\text{BIO12}$)
* **Objective**: Incorporate moisture stability ($\Delta P$) alongside thermal stability ($\Delta T$).
* **Metric**: Relative precipitation change $\Delta P_{\text{rel}} = \frac{| \text{BIO12}_{\text{Present}} - \text{BIO12}_{\text{LGM}} |}{\text{BIO12}_{\text{Present}} + 1.0}$.
* **Synthesis**: Test bivariate stability models:
  $$\log(S + 1) \sim \beta_0 + \beta_1 \Delta T + \beta_2 \Delta P_{\text{rel}} + \beta_3 \text{Roughness} + \epsilon$$

### Phase 2E: Spatial Autocorrelation Diagnostics & Spatial Regression
* **Objective**: Account for spatial autocorrelation in geographical grid residuals.
* **Diagnostics**: Compute Moran's $I$ on OLS model residuals across spatial weight matrices.
* **Models**: If spatial autocorrelation is significant ($p < 0.05$), fit **Spatial Autoregressive Models** (SAR / Spatial Error Models).

### Phase 2F: Equal-Area Equal-Earth Sensitivity Grid
* **Objective**: Address the equal-area limitation of 10-arc-minute geographic latitude/longitude cells.
* **Protocol**: Reproject data to an Equal Earth spatial projection (EPSG:8857) with constant cell area ($100 \text{ km} \times 100 \text{ km}$) and re-evaluate correlation statistics.

### Phase 2G: Multi-Taxonomic Group Comparison
* **Objective**: Expand beyond Class Mammalia to compare stability associations across Amphibia (low mobility), Aves (high mobility), and Squamata.
