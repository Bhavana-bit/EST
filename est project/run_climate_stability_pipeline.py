"""Run the CHELSA 5-degree temperature stability pipeline."""

from src.climate_stability_pipeline import run_pipeline


def main() -> None:
    summary = run_pipeline()
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
