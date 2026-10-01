"""Local paths for the real-data Phase 1 analysis."""

from pathlib import Path

BASE_DIR = Path(__file__).resolve().parent
DATA_DIR = BASE_DIR / "data"
RAW_DATA_DIR = DATA_DIR / "raw"
WORLDCLIM_DIR = RAW_DATA_DIR / "worldclim"
PALEOCLIM_DIR = RAW_DATA_DIR / "paleoclim"
IUCN_DIR = RAW_DATA_DIR / "iucn"
DEM_DIR = RAW_DATA_DIR / "dem"
PROCESSED_DATA_DIR = DATA_DIR / "processed"
GBIF_PROCESSED_DIR = BASE_DIR / "processed_data" / "gbif"
RESULTS_DIR = BASE_DIR / "results"
GBIF_RESULTS_DIR = RESULTS_DIR / "gbif"
MAPS_DIR = RESULTS_DIR / "maps"
TABLES_DIR = RESULTS_DIR / "tables"

GBIF_RAW_DIR = RAW_DATA_DIR / "gbif"
GBIF_CSV = GBIF_PROCESSED_DIR / "gbif_clean.csv"

# Target Taxon Configuration (Confirmed by team)
TARGET_TAXON_NAME = "Aves"
TARGET_TAXON_KEY = 212

WORLDCLIM_PRESENT_BIO1 = WORLDCLIM_DIR / "present" / "bio" / "bio_1"
WORLDCLIM_LGM_BIO1 = WORLDCLIM_DIR / "lgm" / "cclgmbi1.tif"
WORLDCLIM_MID_HOLOCENE_BIO1 = WORLDCLIM_DIR / "mid_holocene" / "ccmidbi1.tif"

PALEOCLIM_CURRENT_BIO1 = PALEOCLIM_DIR / "paleoclim_current_BIO1.tif"
PALEOCLIM_LATE_HOLOCENE_BIO1 = PALEOCLIM_DIR / "paleoclim_late_holocene_BIO1.tif"
PALEOCLIM_LGM_BIO1 = PALEOCLIM_DIR / "paleoclim_LGM_BIO1.tif"
PALEOCLIM_LGM_BIO1_ALIGNED = PALEOCLIM_DIR / "paleoclim_LGM_BIO1_aligned.tif"

# Backward compatibility aliases
CHELSA_CURRENT_BIO1 = PALEOCLIM_CURRENT_BIO1
CHELSA_LATE_HOLOCENE_BIO1 = PALEOCLIM_LATE_HOLOCENE_BIO1
CHELSA_LGM_BIO1 = PALEOCLIM_LGM_BIO1
CHELSA_LGM_BIO1_ALIGNED = PALEOCLIM_LGM_BIO1_ALIGNED

GRID_DEGREES = 5.0
CLIMATE_PROCESSED_DIR = PROCESSED_DATA_DIR / "climate"
CLIMATE_STABILITY_RASTER = CLIMATE_PROCESSED_DIR / "climate_stability_5deg.tif"
CLIMATE_STABILITY_CELLS_CSV = TABLES_DIR / "climate_stability_5deg.csv"
CLIMATE_STABILITY_SUMMARY_CSV = TABLES_DIR / "climate_stability_summary.csv"
CLIMATE_STABILITY_README = CLIMATE_PROCESSED_DIR / "CLIMATE_STABILITY_README.md"

GBIF_SAMPLING_EFFORT_5DEG = GBIF_PROCESSED_DIR / "gbif_sampling_effort_5deg.csv"
GBIF_RICHNESS_5DEG_RASTER = PROCESSED_DATA_DIR / "gbif_species_richness_5deg.tif"
GBIF_CONTINENT_SUMMARY_CSV = GBIF_RESULTS_DIR / "gbif_continent_sampling_summary.csv"
FINAL_ANALYSIS_CSV = TABLES_DIR / "final_biodiversity_climate_5deg.csv"
BIODIV_CLIMATE_STATS_CSV = RESULTS_DIR / "biodiversity_climate_statistics.csv"
BIODIV_CLIMATE_MORAN_CSV = RESULTS_DIR / "biodiversity_climate_moran_i.csv"
BIODIV_CLIMATE_SENSITIVITY_CSV = RESULTS_DIR / "biodiversity_climate_sensitivity_thresholds.csv"
MAIN_SCATTER_PNG = MAPS_DIR / "biodiversity_vs_climate_stability_5deg.png"


