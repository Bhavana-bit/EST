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
RESULTS_DIR = BASE_DIR / "results"
MAPS_DIR = RESULTS_DIR / "maps"
TABLES_DIR = RESULTS_DIR / "tables"

GBIF_CSV = PROCESSED_DATA_DIR / "gbif" / "gbif_clean.csv"

WORLDCLIM_PRESENT_BIO1 = WORLDCLIM_DIR / "present" / "bio" / "bio_1"
WORLDCLIM_LGM_BIO1 = WORLDCLIM_DIR / "lgm" / "cclgmbi1.tif"
WORLDCLIM_MID_HOLOCENE_BIO1 = WORLDCLIM_DIR / "mid_holocene" / "ccmidbi1.tif"

CHELSA_CURRENT_BIO1 = PALEOCLIM_DIR / "paleoclim_current_BIO1.tif"
CHELSA_LGM_BIO1 = PALEOCLIM_DIR / "paleoclim_LGM_BIO1.tif"
