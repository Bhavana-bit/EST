"""Run the real-data Phase 1 analysis."""

import argparse

from src.phase1_pipeline import run_pipeline, validate_inputs


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--validate-only",
        action="store_true",
        help="Validate local inputs without creating analysis outputs.",
    )
    args = parser.parse_args()

    if args.validate_only:
        validate_inputs()
    else:
        run_pipeline()


if __name__ == "__main__":
    main()
