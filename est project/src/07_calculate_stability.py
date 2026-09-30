"""
07_calculate_stability.py
--------------------------
Calculates Palaeoclimate Stability Metrics:
1. Temperature Anomaly ΔT_LGM = |BIO1_present - BIO1_LGM| (°C)
2. Temperature Anomaly ΔT_MH = |BIO1_present - BIO1_MH| (°C)
3. Precipitation Anomaly ΔP_LGM = |BIO12_present - BIO12_LGM| (mm)
4. Relative Precipitation Change ΔP_rel = |P_present - P_LGM| / (P_present + 1.0)
5. Project-Defined Stability Index (I_stability in [0, 1])
"""

import sys
from pathlib import Path
import numpy as np

sys.path.append(str(Path(__file__).resolve().parent.parent))
import config

def calculate_stability():
    print("=" * 80)
    print("CALCULATING PALAEOCLIMATE STABILITY METRICS")
    print("=" * 80)
    
    bio1_present = np.load(config.PROCESSED_DATA_DIR / "bio1_present.npy")
    bio1_lgm = np.load(config.PROCESSED_DATA_DIR / "bio1_lgm.npy")
    bio1_mh = np.load(config.PROCESSED_DATA_DIR / "bio1_mh.npy")
    
    bio12_present = np.load(config.PROCESSED_DATA_DIR / "bio12_present.npy")
    bio12_lgm = np.load(config.PROCESSED_DATA_DIR / "bio12_lgm.npy")
    bio12_mh = np.load(config.PROCESSED_DATA_DIR / "bio12_mh.npy")
    
    # 1. Temperature Anomalies (°C)
    delta_t_lgm = np.abs(bio1_present - bio1_lgm)
    delta_t_mh = np.abs(bio1_present - bio1_mh)
    
    # 2. Precipitation Anomalies (mm) and Relative Change
    delta_p_lgm = np.abs(bio12_present - bio12_lgm)
    delta_p_mh = np.abs(bio12_present - bio12_mh)
    delta_p_rel_lgm = delta_p_lgm / (bio12_present + 1.0)
    
    # 3. Percentile-Truncated Normalization across global grid [0, 1]
    p99_t = np.percentile(delta_t_lgm, 99)
    norm_t = np.clip(delta_t_lgm / p99_t, 0.0, 1.0)
    
    p99_p = np.percentile(delta_p_rel_lgm, 99)
    norm_p = np.clip(delta_p_rel_lgm / p99_p, 0.0, 1.0)
    
    # 4. Project-Defined Stability Index (1 = High Stability, 0 = High Instability)
    stability_index = 1.0 - 0.5 * (norm_t + norm_p)
    
    # Save outputs
    np.save(config.PROCESSED_DATA_DIR / "delta_t_lgm.npy", delta_t_lgm)
    np.save(config.PROCESSED_DATA_DIR / "delta_t_mh.npy", delta_t_mh)
    np.save(config.PROCESSED_DATA_DIR / "delta_p_lgm.npy", delta_p_lgm)
    np.save(config.PROCESSED_DATA_DIR / "delta_p_mh.npy", delta_p_mh)
    np.save(config.PROCESSED_DATA_DIR / "stability_index.npy", stability_index)
    
    print(f"[SUCCESS] Calculated Palaeoclimate Stability Metrics.")
    print(f"  ΔT_LGM Range: {delta_t_lgm.min():.2f}°C to {delta_t_lgm.max():.2f}°C (Mean: {delta_t_lgm.mean():.2f}°C)")
    print(f"  ΔT_MH Range:  {delta_t_mh.min():.2f}°C to {delta_t_mh.max():.2f}°C (Mean: {delta_t_mh.mean():.2f}°C)")
    print(f"  Stability Index Range: {stability_index.min():.3f} to {stability_index.max():.3f} (Mean: {stability_index.mean():.3f})")

if __name__ == "__main__":
    calculate_stability()
