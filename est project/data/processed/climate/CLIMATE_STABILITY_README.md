# PALEOCLIM TEMPERATURE STABILITY INDEX (LGM-DOMINATED)

## Inputs
- PaleoClim v1.2B BIO1 (mean annual temperature) across three periods:
  1. Current baseline (1979–2013)
  2. Late Holocene (Meghalayan, 4.2–0.3 ka BP)
  3. Last Glacial Maximum (LGM, ~21 ka BP)
- Source rasters: `data/raw/paleoclim/paleoclim_current_BIO1.tif`, `paleoclim_late_holocene_BIO1.tif`, `paleoclim_LGM_BIO1_aligned.tif`.
- Stored as integer °C×10; converted to °C by dividing by 10.0 before analysis.

## Grid & Minimum Land Rule
- Common **5°** global grid aligned with GBIF sampling cells (`floor((lat+90)/5)`, `floor((lon+180)/5)`).
- Fine (10 arc-minute) pixels are averaged within each 5° cell for each period.
- Enforces a minimum land pixel threshold (>= 5 fine pixels per 5° cell) to exclude ocean fringe cells.
- Strictly requires valid BIO1 in ALL three periods; ocean and incomplete cells are assigned NoData (NaN).

## Index Structure & LGM Domination
For each 5° cell with valid mean BIO1 in all three periods:
- `bio1_sd_c` = SD(Current, Late Holocene, LGM) in °C (population SD, ddof=0).
- **Temperature stability index** = `1 / (1 + bio1_sd_c)`.

> [!NOTE]
> The index is an **LGM-dominated stability index**. Late Holocene temperature differs from present temperature by only 0.68 °C on average globally, whereas LGM temperature differs by 11.6 °C on average. Consequently, the standard deviation index is 0.9996 correlated (Spearman) with absolute Present–LGM temperature change (|LGM − Present|).

## Cell & Land Area
Geographic 5° cells vary in physical area with latitude:
`cell_area_km2 = R² × (Δ°→rad)² × cos(φ)`, with `R = 6371 km` and `Δ = 5°`.
`estimated_land_area_km2 = cell_area_km2 × (fine_pixel_count / 900)`.
