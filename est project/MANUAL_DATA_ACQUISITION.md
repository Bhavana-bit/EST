# Manual GBIF Data Acquisition Procedure for Aves (taxonKey 212)

> **Important Notice:** The automated cleaning pipeline strictly enforces that raw GBIF occurrences must belong to a confirmed download for class **Aves (taxonKey 212)**. Direct API searching or unconfirmed sample fallback data generation is prohibited.

---

## Step-by-Step GBIF Download Instructions

1. **Log in to GBIF:**
   - Go to [GBIF.org](https://www.gbif.org) and log into your account (create an account if necessary to generate an official download DOI).

2. **Navigate to Occurrence Search:**
   - Go to **Data** $\rightarrow$ **Occurrences** $\rightarrow$ **Search** (`https://www.gbif.org/occurrence/search`).

3. **Apply Required Download Query Filters:**
   - **Scientific Name / Taxon:** Select **Aves** (GBIF taxonKey: `212`).
   - **Location:** Check `hasCoordinate = true` (only records with valid decimal coordinates).
   - **Occurrence Status:** Select `occurrenceStatus = PRESENT`.
   - **Geospatial Quality:** Set `hasGeospatialIssue = false` (exclude flagged coordinate errors).

4. **Request Official Download:**
   - Click **Download** at the top right of the search interface.
   - Choose either **Darwin Core Archive (DwC-A)** or **Simple CSV** format.
   - Submit the download request.

5. **Place Downloaded Files in Project Repository:**
   - Once GBIF processes your request, download the zip archive or extracted CSV.
   - Move/copy the downloaded raw file(s) into:
     ```
     data/raw/gbif/
     ```

6. **Record Metadata in `GBIF_METADATA.md`:**
   - Open `data/raw/gbif/GBIF_METADATA.md` and record the following official provenance details:
     - **GBIF Download Key** (e.g., `0012345-240101120000000`)
     - **GBIF DOI** (e.g., `https://doi.org/10.15468/dl.xxxxxx`)
     - **Download Date** (YYYY-MM-DD)
     - **Query Filters** (e.g., `taxonKey=212`, `hasCoordinate=true`, `occurrenceStatus=PRESENT`, `hasGeospatialIssue=false`)
     - **Total Record Count** (number of raw occurrence records downloaded)
