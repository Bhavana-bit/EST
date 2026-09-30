# GBIF_METADATA.md — GBIF Occurrence Data Provenance & Query Record

**Project Title**: 24. Palaeoclimate Stability and Modern Biodiversity  
**Confirmed Taxon**: Class Aves (Birds)  
**GBIF Taxon Key**: `212` (Phylum Chordata, Class Aves)  
**Date of Query / Record**: 2026-09-30  
**API Endpoint**: `https://api.gbif.org/v1/occurrence/search`

---

## 1. QUERY FILTER PARAMETERS

| Parameter | Value | Rationale |
| :--- | :--- | :--- |
| `taxonKey` | `212` | Class Aves (Bird species) |
| `hasCoordinate` | `true` | Requires valid latitude and longitude coordinates |
| `hasGeospatialIssue` | `false` | Excludes records with flagged geospatial anomalies (e.g. invalid bounds, country centroids) |
| `occurrenceStatus` | `PRESENT` | Excludes species absence records |
| `basisOfRecord` | `PRESERVED_SPECIMEN`, `HUMAN_OBSERVATION`, `OBSERVATION` | Verified empirical observation and specimen records |
| `coordinateUncertaintyInMeters` | `<= 10000` (10 km) | Ensures spatial precision within the 10 arc-minute grid cell resolution (~18.5 km) |

---

## 2. BIODIVERSITY VARIABLE CALCULATION

* **Primary Metric**: **Unique Species Count ($S_i$)** per 10 arc-minute spatial grid cell.
  $$S_i = \text{Count of UNIQUE species keys / names in cell } i$$
* **Raw Record Count Exclusion**: Raw occurrence record density ($N_i$) is recorded ONLY as a diagnostic covariate for sampling effort control. It is **NEVER** used directly as species richness.

---

## 3. KNOWN SAMPLING LIMITATIONS

1. **Geographic Sampling Bias**: GBIF bird records reflect higher sampling density near roads, urban centers, and birdwatching hotspots (e.g., eBird contributions).
2. **Unequal Cell Area**: 10-arc-minute geographic grid cells shrink towards the poles.
3. **Taxonomic Coverage**: Focuses exclusively on Class Aves for Phase 1.
