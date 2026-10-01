"""Run final descriptive biodiversity–climate statistics."""

from src.biodiversity_climate_stats import run_analysis


def main() -> None:
    statistics, sensitivity, moran = run_analysis()
    print("=== Primary Spearman and Covariate Controls (n=369) ===")
    print(
        statistics[
            [
                "model",
                "n",
                "correlation_coefficient",
                "p_value",
                "covariates_controlled",
            ]
        ].to_string(index=False)
    )
    print("\n=== Sampling Effort / Record Threshold Sensitivity ===")
    print(
        sensitivity[
            [
                "min_occurrence_threshold",
                "n_cells",
                "spearman_rho",
                "spearman_p_value",
                "partial_r_control_effort_and_abs_lat",
                "partial_p_control_effort_and_abs_lat",
            ]
        ].to_string(index=False)
    )
    print("\n=== Spatial Autocorrelation (Moran's I on 5° grid) ===")
    print(
        moran[
            [
                "variable",
                "n",
                "morans_i",
                "morans_i_p_value",
                "interpretation",
            ]
        ].to_string(index=False)
    )


if __name__ == "__main__":
    main()
