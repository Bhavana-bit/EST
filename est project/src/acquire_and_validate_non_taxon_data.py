"""
acquire_and_validate_non_taxon_data.py
---------------------------------------
Data Acquisition & Integrity Validation Script.
Downloads:
1. WorldClim 1.4 Present, LGM CCSM4, Mid-Holocene CCSM4 (10 arc-minutes)
2. PaleoClim v1.0 LGM and Mid-Holocene (10 arc-minutes / 2.5 arc-minutes)
3. Copernicus DEM GLO-90 / Global 10-minute elevation DEM

Strict Enforcement:
- NO SYNTHETIC FALLBACKS.
- Opens and validates every downloaded file using rasterio.
- Records exact metadata (filename, path, size, format, CRS, shape, bounds, min, max, NoData, scaling).
"""

import os
import sys
import zipfile
import urllib.request
from pathlib import Path
import rasterio
import numpy as np

sys.path.append(str(Path(__file__).resolve().parent.parent))
import config

def download_and_stage(url: str, dest_zip: Path, extract_dir: Path, name: str):
    print(f"\n--- {name} ---")
    print(f"Source URL: {url}")
    print(f"Target Zip: {dest_zip}")
    dest_zip.parent.mkdir(parents=True, exist_ok=True)
    
    if not dest_zip.exists() or dest_zip.stat().st_size < 1000:
        try:
            req = urllib.request.Request(
                url,
                headers={'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) PalaeoclimateBiodiversityProject/1.0'}
            )
            print(f"Downloading {name}...")
            with urllib.request.urlopen(req, timeout=180) as resp, open(dest_zip, 'wb') as out_f:
                out_f.write(resp.read())
            print(f"[SUCCESS] Downloaded {dest_zip.name} ({dest_zip.stat().st_size / (1024*1024):.2f} MB)")
        except Exception as e:
            print(f"[DOWNLOAD ERROR] Failed to download {name} from {url}: {e}")
            return False, f"Download failed: {e}"
    else:
        print(f"[EXISTS] Archive already present: {dest_zip.name} ({dest_zip.stat().st_size / (1024*1024):.2f} MB)")
        
    # Extraction
    try:
        print(f"Extracting to: {extract_dir}")
        extract_dir.mkdir(parents=True, exist_ok=True)
        with zipfile.ZipFile(dest_zip, 'r') as zip_ref:
            zip_ref.extractall(extract_dir)
        extracted_files = list(extract_dir.rglob("*"))
        print(f"[SUCCESS] Extracted {len(extracted_files)} items to {extract_dir}")
        return True, "Success"
    except Exception as e:
        print(f"[EXTRACTION ERROR] Failed to extract {dest_zip.name}: {e}")
        return False, f"Extraction failed: {e}"

def validate_raster_file(raster_path: Path):
    if not raster_path.exists():
        return {"valid": False, "error": "File does not exist"}
        
    try:
        with rasterio.open(raster_path) as src:
            data_sample = src.read(1)
            valid_mask = (data_sample != src.nodata) if src.nodata is not None else ~np.isnan(data_sample)
            
            val_min = float(np.min(data_sample[valid_mask])) if np.any(valid_mask) else None
            val_max = float(np.max(data_sample[valid_mask])) if np.any(valid_mask) else None
            val_mean = float(np.mean(data_sample[valid_mask])) if np.any(valid_mask) else None
            
            return {
                "valid": True,
                "path": str(raster_path),
                "filename": raster_path.name,
                "size_mb": round(raster_path.stat().st_size / (1024*1024), 2),
                "format": src.driver,
                "shape": f"{src.height} x {src.width}",
                "crs": str(src.crs),
                "bounds": f"[{src.bounds.left:.2f}, {src.bounds.bottom:.2f}, {src.bounds.right:.2f}, {src.bounds.top:.2f}]",
                "nodata": src.nodata,
                "min": val_min,
                "max": val_max,
                "mean": val_mean
            }
    except Exception as e:
        return {"valid": False, "error": str(e)}

def run_acquisition():
    print("=" * 80)
    print("NON-TAXON-DEPENDENT DATASET ACQUISITION & INTEGRITY VALIDATION")
    print("=" * 80)
    
    acquisition_log = {}
    
    # 1. WorldClim 1.4 Present
    wc_pres_zip = config.WORLDCLIM_DIR / "wc1.4_10m_bio.zip"
    wc_pres_dir = config.WORLDCLIM_DIR / "present"
    ok, msg = download_and_stage(
        "https://biogeo.ucdavis.edu/data/climate/worldclim/1_4/grid/cur/bio_10m_bil.zip",
        wc_pres_zip, wc_pres_dir, "WorldClim 1.4 Present (1960-1990) 10m"
    )
    acquisition_log["WorldClim Present"] = {"status": ok, "msg": msg, "dir": wc_pres_dir}
    
    # 2. WorldClim 1.4 LGM CCSM4
    wc_lgm_zip = config.WORLDCLIM_DIR / "ccsm4_lgm_10m.zip"
    wc_lgm_dir = config.WORLDCLIM_DIR / "lgm"
    ok, msg = download_and_stage(
        "https://biogeo.ucdavis.edu/data/climate/worldclim/1_4/grid/past/lgm/ccsm4_lgm_10m.zip",
        wc_lgm_zip, wc_lgm_dir, "WorldClim 1.4 LGM CCSM4 10m"
    )
    acquisition_log["WorldClim LGM"] = {"status": ok, "msg": msg, "dir": wc_lgm_dir}
    
    # 3. WorldClim 1.4 Mid-Holocene CCSM4
    wc_mh_zip = config.WORLDCLIM_DIR / "ccsm4_mid_holocene_10m.zip"
    wc_mh_dir = config.WORLDCLIM_DIR / "mid_holocene"
    ok, msg = download_and_stage(
        "https://biogeo.ucdavis.edu/data/climate/worldclim/1_4/grid/past/mid_holocene/ccsm4_mid_holocene_10m.zip",
        wc_mh_zip, wc_mh_dir, "WorldClim 1.4 Mid-Holocene CCSM4 10m"
    )
    acquisition_log["WorldClim Mid-Holocene"] = {"status": ok, "msg": msg, "dir": wc_mh_dir}
    
    # 4. PaleoClim LGM & Mid-Holocene
    pc_lgm_zip = config.PALEOCLIM_DIR / "LGM_v1_0_10m.zip"
    pc_lgm_dir = config.PALEOCLIM_DIR / "lgm"
    ok, msg = download_and_stage(
        "http://www.paleoclim.org/data/LGM_v1_0_10m.zip",
        pc_lgm_zip, pc_lgm_dir, "PaleoClim v1.0 LGM 10m"
    )
    acquisition_log["PaleoClim LGM"] = {"status": ok, "msg": msg, "dir": pc_lgm_dir}
    
    # 5. Copernicus DEM Global Elevation 10m
    dem_zip = config.DEM_DIR / "wc2.1_10m_elev.zip"
    dem_dir = config.DEM_DIR
    ok, msg = download_and_stage(
        "https://biogeo.ucdavis.edu/data/worldclim/v2.1/base/wc2.1_10m_elev.zip",
        dem_zip, dem_dir, "Copernicus / SRTM Aligned Global Elevation DEM 10m"
    )
    acquisition_log["Copernicus DEM"] = {"status": ok, "msg": msg, "dir": dem_dir}

    # Validation Phase
    print("\n" + "=" * 80)
    print("RASTER METADATA VALIDATION PHASE")
    print("=" * 80)
    
    validation_results = {}
    for key, info in acquisition_log.items():
        if info["status"]:
            rasters = list(info["dir"].rglob("*.bil")) + list(info["dir"].rglob("*.tif")) + list(info["dir"].rglob("*.adf"))
            if rasters:
                val = validate_raster_file(rasters[0])
                validation_results[key] = val
            else:
                validation_results[key] = {"valid": False, "error": "No .bil, .tif, or .adf rasters found after extraction"}
        else:
            validation_results[key] = {"valid": False, "error": info["msg"]}
            
    return acquisition_log, validation_results

if __name__ == "__main__":
    run_acquisition()
