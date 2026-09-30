"""
inspect_est_repo.py
-------------------
Inspects the actual datasets and processed files in EST/Palaeoclimate_Biodiversity_Project/
"""

from pathlib import Path
import zipfile
import pandas as pd

base_path = Path("EST/Palaeoclimate_Biodiversity_Project")

print("=" * 80)
print("INSPECTING EST PROJECT REPOSITORY DATASETS")
print("=" * 80)

# 1. GBIF files
gbif_dir = base_path / "data" / "gbif"
print("\n--- 1. GBIF Occurrence Files ---")
for f in gbif_dir.glob("*.csv"):
    df = pd.read_csv(f, nrows=5)
    size_mb = f.stat().st_size / (1024*1024)
    total_lines = sum(1 for _ in open(f, encoding='utf-8', errors='ignore')) - 1
    print(f"\nFile: {f.name} ({size_mb:.2f} MB)")
    print(f"  Total Rows: {total_lines}")
    print(f"  Columns: {list(df.columns)}")

# 2. WorldClim files
wc_dir = base_path / "data" / "worldclim"
print("\n\n--- 2. WorldClim 1.4 Zip Archives ---")
for f in wc_dir.glob("*.zip"):
    size_mb = f.stat().st_size / (1024*1024)
    print(f"\nFile: {f.name} ({size_mb:.2f} MB)")
    with zipfile.ZipFile(f, 'r') as z:
        members = z.namelist()
        print(f"  Contains {len(members)} entries (Sample: {members[:5]})")

# 3. PaleoClim files
pc_dir = base_path / "data" / "paleoclim"
print("\n\n--- 3. PaleoClim Files ---")
for f in pc_dir.glob("*"):
    size_mb = f.stat().st_size / (1024*1024)
    print(f"  - {f.name} ({size_mb:.2f} MB)")

# 4. Processed Data
proc_dir = base_path / "processed_data"
print("\n\n--- 4. Processed Data Directory ---")
for root, dirs, files in os.walk(proc_dir):
    for f in files:
        fp = Path(root) / f
        print(f"  - {fp.relative_to(proc_dir)} ({fp.stat().st_size / 1024:.1f} KB)")
