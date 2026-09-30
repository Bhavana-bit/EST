# WORLDCLIM_METADATA.md — WorldClim 1.4 Provenance & Verification Record

**Project Title**: 24. Palaeoclimate Stability and Modern Biodiversity  
**Dataset Version**: WorldClim Version 1.4 (Release 1.4)  
**Spatial Resolution**: 10 arc-minutes (~18.5 km at equator)  
**CRS**: EPSG:4326 (WGS 84 Geographic Latitude/Longitude)  
**Spatial Extent**: Global [-180.0, 180.0, -90.0, 90.0]  
**Grid Shape**: 1080 rows x 2160 columns  

---

## 1. OFFICIAL SOURCE URLS & PERIOD SPECIFICATIONS

| Period / Scenario | Temporal Period | GCM Model | Official Source URL | File Format |
| :--- | :--- | :--- | :--- | :--- |
| **Present Baseline** | 1960–1990 Baseline | Observed Climatology | `https://biogeo.ucdavis.edu/data/climate/worldclim/1_4/grid/cur/bio_10m_bil.zip` | BIL / GeoTIFF archive |
| **Last Glacial Maximum (LGM)** | ~22,000 BP | CCSM4 | `https://biogeo.ucdavis.edu/data/climate/worldclim/1_4/grid/past/lgm/ccsm4_lgm_10m.zip` | BIL / GeoTIFF archive |
| **Mid-Holocene** | ~6,000 BP | CCSM4 | `https://biogeo.ucdavis.edu/data/climate/worldclim/1_4/grid/past/mid_holocene/ccsm4_mid_holocene_10m.zip` | BIL / GeoTIFF archive |

---

## 2. VARIABLE & UNIT VERIFICATION

| Variable | Official Name | Raw Data Unit / Encoding | Processed Unit | Conversion Rule |
| :--- | :--- | :--- | :--- | :--- |
| `BIO1` | Annual Mean Temperature | Integer (°C × 10) | Floating-point °C | `BIO1_degc = raw_bio1 / 10.0` |
| `BIO12` | Annual Precipitation | Integer (mm/year) | mm/year | No conversion required |

---

## 3. PROJECT-DEFINED TEMPERATURE DIFFERENCE METRIC

$$\Delta T_{\text{LGM}} = | \text{BIO1}_{\text{Present}} - \text{BIO1}_{\text{LGM}} | \quad (°\text{C})$$

* **Interpretation**: Smaller absolute temperature difference $\Delta T$ represents greater long-term thermal stability under this project-defined metric.
