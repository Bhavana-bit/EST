"""
run_phase1.py
-------------
Master Execution Script for PHASE 1 (Mid-Evaluation Deliverables).
Executes:
1. phase1_01_grid_and_climate.py  -> Common 10' grid & Present vs LGM BIO1 ΔT
2. phase1_02_gbif_richness.py      -> GBIF Mammalia unique species richness (S)
3. phase1_03_analysis_and_plots.py -> Statistics, MAP 1, MAP 2, & Scatter Plot
"""

import sys
from pathlib import Path

BASE_DIR = Path(__file__).resolve().parent
sys.path.append(str(BASE_DIR))

from src.phase1_01_grid_and_climate import process_phase1_climate
from src.phase1_02_gbif_richness import process_phase1_gbif
from src.phase1_03_analysis_and_plots import run_phase1_analysis

def execute_phase1():
    print("=" * 80)
    print("EXECUTING PHASE 1 MID-EVALUATION PIPELINE")
    print("=" * 80)
    
    print("\n[Step 1/3] Processing WorldClim 10' Grid & Temperature Difference (|ΔT|)...")
    process_phase1_climate()
    
    print("\n[Step 2/3] Processing GBIF Mammalia Unique Species Richness (S)...")
    process_phase1_gbif()
    
    print("\n[Step 3/3] Running Statistical Analysis & Rendering Phase 1 Maps/Plots...")
    run_phase1_analysis()
    
    print("\n" + "=" * 80)
    print("PHASE 1 EXECUTION COMPLETE!")
    print("All deliverables ready for mid-evaluation in results/maps/ and results/tables/")
    print("=" * 80)

if __name__ == "__main__":
    execute_phase1()
