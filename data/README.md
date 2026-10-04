# Data

## Files read by the scripts

| File | Used by | Contents |
|---|---|---|
| `Day_cleaned.csv` | 01, 03, 04 | Day tomograms: one block of rows per tomogram (8 tomograms, 102 rows). Columns are per-tomogram features; measurement columns (globule/vesicle diameters) continue on following rows, so most cells in a tomogram's later rows are blank. |
| `Night_cleaned.csv` | 01, 03, 04 | Night tomograms (64 rows). Same feature columns as Day, plus `Nickname`, `Grid`, `Tomogram`, and `NumHat`. **`NumHat` is incomplete** (23 tomograms / 38 hats vs 27 / 49 in `hat_counts_by_tomogram.csv`) and is no longer used by any script. |
| `lamellae_thickthincombined.csv` | 02 | 322 lamellae width measurements (nm): `Thick` (n = 213), `Thin` (n = 109), `Combined` (all 322). |
| `hat_counts_by_tomogram.csv` | 01, 03, 04 | One row per tomogram (72 total): `state` (Day / Evening / Night / Morning), `tomogram_position` (order within that period), `hats` (hat structures counted), `vol` (matrix volume, µm³). Totals: Day 1 / 1.601, Evening 24 / 2.957, Night 49 / 4.979, Morning 31 / 5.859. |

### Columns used (names as R's `read.csv` sees them)

| Column | Meaning | Values |
|---|---|---|
| `Internal.Lamellae.Type` | Lamellae class | `Thin`, `Thick`, `Thin and Thick`, `Thick and Thin`, `None`, blank |
| `Globules.Observed` / `Globule.Diameters..Widest.` | Globule count per tomogram / widest diameter (nm), one per row | numeric |
| `Spiky.Vesicles.Observed` / `Vesicle.Diameter..widest.` | Spiky-vesicle count per tomogram / widest diameter (nm) | numeric |
| `Hat.Structure.Observed` (Day) / `Hat.Structures.Observed.` (Night) | Whether a hat structure was seen | `Yes`, `No`, `Maybe`, blank |
| `NumHat` (Night only) | Number of hat structures (superseded by `hat_counts_by_tomogram.csv`) | numeric |

**Nitroplasts vs tomograms:** several tomograms can cover the same nitroplast. The per-nitroplast normalizations in
script 03 divide by distinct *nitroplast* counts (Day 8, Night 33, Morning 28; set in `scripts/00_config.R`), which are
not the same as the tomogram counts in `hat_counts_by_tomogram.csv` (Day 8, Evening 18, Night 27, Morning 19).

**Edit note:** in `Day_cleaned.csv`, Day tomogram 1's `Hat Structure Observed` was changed from `No` to `Yes` (the single Day hat structure).

**Encoding:** `lamellae_thickthincombined.csv` and `Night_cleaned.csv` begin with a UTF-8 byte-order mark
(typical of Excel exports). The scripts read all data through `read_data()` in `scripts/00_config.R`, which
handles this. If you read these files directly with `read.csv()`, the first column is misnamed (e.g. `X...Thick`).

## Data hardcoded in scripts (not in `data/`)

| Data | Script | Note |
|---|---|---|
| Evening vesicle diameters (39 values), Morning vesicle diameters (3 values: 54, 50, 50) | 04 (03 for Morning) | |
| Vesicle counts and matrix volumes per period (Day 88 / 1.601 µm³, Evening 39 / 2.957, Night 8 / 4.979, Morning 3 / 5.859) | 04 | Day and Night counts are checked against the CSVs in the script |
| Thylakoid widths: Day (81), Night (113), Morning (125), Evening (68) | 05 | |
| Membrane contact site percentages | 06 | |

Moving these into CSVs under `data/` would make the repo easier to audit; ask if you'd like that done.

## Hosting

TODO: decide where data live (this repo / Zenodo / institutional storage / not shareable) and update this file and the main README.
