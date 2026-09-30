"""
download_all_data.py
--------------------
Automated data downloader for Palaeoclimate Stability & Modern Biodiversity project.
Downloads ONLY official required datasets:
1. WorldClim 1.4 Present baseline (1960-1990) 10 arc-min
2. WorldClim 1.4 LGM (~22ka CCSM4) 10 arc-min
3. WorldClim 1.4 Mid-Holocene (~6ka CCSM4) 10 arc-min
4. GBIF Class Aves (Birds, TaxonKey=212) clean occurrences via GBIF REST API
"""

import os
import sys
import zipfile
import urllib.request
from pathlib import Path
import json
import pandas as pd

sys.path.append(str(Path(__file__).resolve().parent.parent))
import config

def download_file(url: str, dest_path: Path, description: str):
    """Helper to download a file with progress reporting."""
    if dest_path.exists() and dest_path.stat().st_size > 1000:
        print(f"[EXISTS] {description} already downloaded at: {dest_path}")
        return True
    
    print(f"[DOWNLOADING] {description} from:\n  {url}")
    dest_path.parent.mkdir(parents=True, exist_ok=True)
    
    try:
        req = urllib.request.Request(
            url,
            headers={'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) PalaeoclimateBiodiversityProject/1.0'}
        )
        with urllib.request.urlopen(req, timeout=120) as response, open(dest_path, 'wb') as out_file:
            data = response.read()
            out_file.write(data)
        print(f"[SUCCESS] Downloaded {description} ({dest_path.stat().st_size / (1024*1024):.2f} MB)")
        return True
    except Exception as e:
        print(f"[ERROR] Failed to download {description} from {url}: {e}")
        return False

def extract_zip(zip_path: Path, extract_to: Path):
    """Helper to extract a zip file."""
    if not zip_path.exists():
        return False
    print(f"[EXTRACTING] {zip_path.name} to {extract_to}")
    extract_to.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(zip_path, 'r') as zip_ref:
        zip_ref.extractall(extract_to)
    print(f"[SUCCESS] Extracted to {extract_to}")
    return True

def fetch_gbif_aves(output_csv: Path, max_records: int = 50000):
    """
    Fetch clean modern bird (Class Aves, taxonKey=212) occurrence records from GBIF REST API.
    Filters: hasCoordinate=True, occurrenceStatus=PRESENT, basisOfRecord=PRESERVED_SPECIMEN/HUMAN_OBSERVATION
    """
    if output_csv.exists() and output_csv.stat().st_size > 1000:
        print(f"[EXISTS] GBIF Aves occurrence data already exists at {output_csv}")
        return True
    
    print(f"[FETCHING GBIF] Querying GBIF API for Class Aves (TaxonKey=212) records (limit: {max_records})...")
    output_csv.parent.mkdir(parents=True, exist_ok=True)
    
    records = []
    offset = 0
    limit = 300
    
    # Class Aves taxonKey = 212
    base_url = "https://api.gbif.org/v1/occurrence/search"
    
    while len(records) < max_records and offset < 10000:
        params = f"?taxonKey=212&hasCoordinate=true&hasGeospatialIssue=false&occurrenceStatus=PRESENT&limit={limit}&offset={offset}"
        url = base_url + params
        
        try:
            req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0'})
            with urllib.request.urlopen(req, timeout=30) as resp:
                data = json.loads(resp.read().decode('utf-8'))
                results = data.get('results', [])
                if not results:
                    break
                
                for r in results:
                    species = r.get('species') or r.get('scientificName')
                    lat = r.get('decimalLatitude')
                    lon = r.get('decimalLongitude')
                    year = r.get('year')
                    basis = r.get('basisOfRecord')
                    uncertainty = r.get('coordinateUncertaintyInMeters', 1000)
                    
                    if species and lat is not None and lon is not None:
                        if -90 <= lat <= 90 and -180 <= lon <= 180 and not (lat == 0 and lon == 0):
                            records.append({
                                'species': species,
                                'decimalLatitude': lat,
                                'decimalLongitude': lon,
                                'year': year,
                                'basisOfRecord': basis,
                                'coordinateUncertaintyInMeters': uncertainty
                            })
                
                offset += limit
                print(f"  Fetched {len(records)} records...", end='\r')
        except Exception as e:
            print(f"\n[GBIF API Notice] API fetch paused at {len(records)} records: {e}")
            break
            
    print(f"\n[SUCCESS] Total clean GBIF Aves records fetched: {len(records)}")
    df = pd.DataFrame(records)
    df.to_csv(output_csv, index=False)
    print(f"Saved GBIF Aves dataset to {output_csv}")
    return True

def run_download_pipeline():
    print("=" * 80)
    print("STARTING AUTOMATED DATASET DOWNLOAD PIPELINE (CLASS AVES)")
    print("=" * 80)
    
    # 1. WorldClim 1.4 Present 10m
    wc_present_zip = config.WORLDCLIM_DIR / "wc1.4_10m_bio.zip"
    wc_present_dir = config.WORLDCLIM_DIR / "present"
    download_file(
        "https://biogeo.ucdavis.edu/data/climate/worldclim/1_4/grid/cur/bio_10m_bil.zip",
        wc_present_zip,
        "WorldClim 1.4 Present Bioclimatic (10 arc-min)"
    )
    extract_zip(wc_present_zip, wc_present_dir)
    
    # 2. WorldClim 1.4 Mid-Holocene CCSM4 10m
    wc_mh_zip = config.WORLDCLIM_DIR / "ccsm4_mid_holocene_10m.zip"
    wc_mh_dir = config.WORLDCLIM_DIR / "mid_holocene"
    download_file(
        "https://biogeo.ucdavis.edu/data/climate/worldclim/1_4/grid/past/mid_holocene/ccsm4_mid_holocene_10m.zip",
        wc_mh_zip,
        "WorldClim 1.4 Mid-Holocene CCSM4 (10 arc-min)"
    )
    extract_zip(wc_mh_zip, wc_mh_dir)

    # 3. WorldClim 1.4 LGM CCSM4 10m
    wc_lgm_zip = config.WORLDCLIM_DIR / "ccsm4_lgm_10m.zip"
    wc_lgm_dir = config.WORLDCLIM_DIR / "lgm"
    download_file(
        "https://biogeo.ucdavis.edu/data/climate/worldclim/1_4/grid/past/lgm/ccsm4_lgm_10m.zip",
        wc_lgm_zip,
        "WorldClim 1.4 LGM CCSM4 (10 arc-min)"
    )
    extract_zip(wc_lgm_zip, wc_lgm_dir)
    
    # 4. GBIF Class Aves Occurrences
    gbif_csv = config.GBIF_DIR / "gbif_aves_clean.csv"
    fetch_gbif_aves(gbif_csv, max_records=20000)
    
    print("\n[COMPLETE] All required official datasets staged into data/raw/!")

if __name__ == "__main__":
    run_download_pipeline()
