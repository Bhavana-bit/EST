"""
Configuration file for Palaeoclimate Stability and Modern Biodiversity Project.
Uses OS-independent pathlib paths.
Taxon confirmed: Class Aves (Bird species, GBIF TaxonKey = 212)
"""

from pathlib import Path

# Base Directory
BASE_DIR = Path(__file__).resolve().parent

# Data Directories
DATA_DIR = BASE_DIR / "data"
RAW_DATA_DIR = DATA_DIR / "raw"
PROCESSED_DATA_DIR = DATA_DIR / "processed"

# Subdirectories for raw data
WORLDCLIM_DIR = RAW_DATA_DIR / "worldclim"
PALEOCLIM_DIR = RAW_DATA_DIR / "paleoclim"
GBIF_DIR = RAW_DATA_DIR / "gbif"
IUCN_DIR = RAW_DATA_DIR / "iucn"
DEM_DIR = RAW_DATA_DIR / "dem"

# Results Directories
RESULTS_DIR = BASE_DIR / "results"
MAPS_DIR = RESULTS_DIR / "maps"
TABLES_DIR = RESULTS_DIR / "tables"
STATS_DIR = RESULTS_DIR / "statistics"

# Grid Configuration
# WorldClim 10 arc-minute resolution
GRID_RESOLUTION_DEG = 10.0 / 60.0  # 0.16666666666666666 degrees
TARGET_CRS = "EPSG:4326"
GRID_SHAPE = (1080, 2160)  # (rows, cols) for global [-180, 180, -90, 90]

# Confirmed Taxon Specifications
TAXON_NAME = "Aves"
TAXON_KEY = 212  # GBIF TaxonKey for Class Aves

# Ensure directories exist
for directory in [
    RAW_DATA_DIR, PROCESSED_DATA_DIR,
    WORLDCLIM_DIR, PALEOCLIM_DIR, GBIF_DIR, IUCN_DIR, DEM_DIR,
    MAPS_DIR, TABLES_DIR, STATS_DIR
]:
    directory.mkdir(parents=True, exist_ok=True)
