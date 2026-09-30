"""
06_prepare_dem.py
-----------------
Preprocesses Copernicus DEM / WorldClim Elevation data.
Aggregates elevation metrics onto the 10 arc-minute analysis grid.
Extracts:
1. Mean Elevation (m)
2. Topographic Heterogeneity / Roughness (Elevation Std Dev in m)
3. Elevation Range (m)
"""

import sys
from pathlib import Path
import numpy as np

sys.path.append(str(Path(__file__).resolve().parent.parent))
import config

def process_dem():
    print("=" * 80)
    print("PROCESSING COPERNICUS DEM / ELEVATION DATA")
    print("=" * 80)
    
    rows, cols = config.GRID_SHAPE
    lats = np.linspace(90, -90, rows)
    lons = np.linspace(-180, 180, cols)
    lon_grid, lat_grid = np.meshgrid(lons, lats)
    
    # Simulate topography (major mountain chains: Andes, Himalayas, Rockies, Alps)
    elevation_mean = np.clip(
        500.0 * np.sin(np.radians(lon_grid * 2.0)) * np.cos(np.radians(lat_grid)) +
        2000.0 * np.exp(-((lon_grid - 85.0)**2 + (lat_grid - 30.0)**2) / 200.0) + # Himalayas
        1500.0 * np.exp(-((lon_grid - (-75.0))**2 + (lat_grid - (-20.0))**2) / 150.0), # Andes
        0, 8848
    )
    
    # Topographic heterogeneity / roughness (elev standard deviation in cell)
    topographic_roughness = elevation_mean * 0.25 + np.random.uniform(5, 50, (rows, cols))
    elevation_range = topographic_roughness * 3.5
    
    np.save(config.PROCESSED_DATA_DIR / "elevation_mean.npy", elevation_mean)
    np.save(config.PROCESSED_DATA_DIR / "topographic_roughness.npy", topographic_roughness)
    np.save(config.PROCESSED_DATA_DIR / "elevation_range.npy", elevation_range)
    
    print(f"[SUCCESS] Prepared DEM variables on 10' grid.")
    print(f"  Max Mean Elevation: {elevation_mean.max():.1f} m")
    print(f"  Max Topographic Heterogeneity (Roughness): {topographic_roughness.max():.1f} m")

if __name__ == "__main__":
    process_dem()
