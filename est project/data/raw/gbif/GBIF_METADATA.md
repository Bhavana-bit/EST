# Local GBIF Data Notes

## Files Used

The Phase 1 pipeline reads `data/processed/gbif/gbif_clean.csv`. The copied file contains 16,764 rows, all marked `SPECIES`, with `speciesKey`, scientific name, latitude, and longitude fields. Its scientific names show mixed taxa, including plants and multiple animal groups.

The file does not contain class/order hierarchy, occurrence year, basis of record, coordinate uncertainty, or a GBIF download DOI. The local claim that this dataset is Aves with taxon key 212 is not supported by the records and must not be used to describe the analyzed data. Download query and filtering provenance are NOT VERIFIED from the local CSV.

## Phase 1 Measures

- Observed species richness per grid cell: number of unique non-null `speciesKey` values.
- Occurrence count per grid cell: number of records assigned to the cell; used as a sampling-effort indicator.

Use the term “GBIF-observed species richness across recorded taxa.” This occurrence-based measure does not represent true or complete biodiversity. Records are unevenly distributed geographically and taxonomically, and occurrence counts reflect observation and reporting effort as well as species presence.