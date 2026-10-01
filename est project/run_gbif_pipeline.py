"""Run the reproducible GBIF cleaning pipeline."""

import argparse

from src.gbif_clean import print_pipeline_report, run_gbif_pipeline


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--clean-only",
        action="store_true",
        help="Run GBIF cleaning and 5-degree sampling-effort aggregation only.",
    )
    args = parser.parse_args()

    metrics = run_gbif_pipeline()
    print_pipeline_report(metrics)

    if not args.clean_only:
        from src.phase1_pipeline import run_pipeline

        print("\nRunning Phase 1 climate-richness analysis on cleaned GBIF records...")
        run_pipeline()


if __name__ == "__main__":
    main()
