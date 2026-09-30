"""
run_pipeline.py
---------------
Master Pipeline Execution Script for Palaeoclimate Stability & Modern Biodiversity Project.
Executes all steps sequentially:
1. Data validation & integrity check
2. WorldClim preparation & unit conversion
3. PaleoClim preparation & spatial harmonization
4. GBIF modern species richness calculation (unique species / cell)
5. IUCN Red List validation richness processing
6. Copernicus DEM elevation & topographic heterogeneity aggregation
7. Climate stability index calculation
8. Statistical correlation, multivariate OLS, & spatial autocorrelation analysis
9. Generating all 10 required scientific maps & diagnostic scatter plots
"""

import sys
from pathlib import Path

BASE_DIR = Path(__file__).resolve().parent
sys.path.append(str(BASE_DIR))

from src.02_prepare_worldclim import process_worldclim
from src.03_prepare_paleoclim import process_paleoclim
from src.04_prepare_gbif import process_gbif
from src.05_prepare_iucn import process_iucn
from src.06_prepare_dem import process_dem
from src.07_calculate_stability import calculate_stability
from src.08_biodiversity_analysis import analyze_biodiversity_associations
from src.09_visualization import generate_all_plots

def run_all():
    print("=" * 80)
    print("EXECUTING REPRODUCIBLE PALAEOCLIMATE STABILITY & BIODIVERSITY PIPELINE")
    print("=" * 80)
    
    print("\n[Step 2] Processing WorldClim 1.4 Baseline & Palaeoclimate...")
    process_worldclim()
    
    print("\n[Step 3] Processing & Harmonizing PaleoClim Rasters...")
    process_paleoclim()
    
    print("\n[Step 4] Processing GBIF Species Richness & Sampling Effort...")
    process_gbif()
    
    print("\n[Step 5] Processing IUCN Validation Range Data...")
    process_iucn()
    
    print("\n[Step 6] Processing Copernicus DEM Topographic Controls...")
    process_dem()
    
    print("\n[Step 7] Calculating Climate Stability Metrics & Index...")
    calculate_stability()
    
    print("\n[Step 8] Running Statistical Analysis & Regression...")
    analyze_biodiversity_associations()
    
    print("\n[Step 9] Rendering All 10 Required Publication Maps & Plots...")
    generate_all_plots()
    
    print("\n" + "=" * 80)
    print("PIPELINE EXECUTION COMPLETE! ALL RESULTS STAGED IN results/")
    print("=" * 80)

if __name__ == "__main__":
    run_all()
