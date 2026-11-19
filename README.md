# Analyzing Microplastic Pollution in the Northwestern Mediterranean

CSE 163 Part 2 EDA. Author: Celine Wu.

Public repository:
[ywu222-ctrl/cse163-microplastics](https://github.com/ywu222-ctrl/cse163-microplastics).

## Data

The original data is NOAA NCEI accession **0253923**, collected by the
University of Barcelona team and citizen scientists. The Excel file in this
project is an unchanged copy of that archived dataset. All reported
measurements and figures use these records.

- [NOAA metadata and sampling information](https://www.ncei.noaa.gov/access/metadata/landing-page/bin/iso?id=gov.noaa.nodc%3A0253923)
- [Download the Excel file from the GitHub mirror](https://github.com/ywu222-ctrl/cse163-microplastics/raw/refs/heads/main/NOAA_0253923_original.xlsx)
- [Original NOAA Excel file (backup download)](https://www.ncei.noaa.gov/data/oceans/archive/arc0203/0253923/1.1/data/0-data/MP_data_Sanchez_Vidal_S2N.xlsx)

The GitHub mirror contains the same original Excel file; it is not another
dataset. Use the NOAA metadata page for the study and sampling details.

1. A copy is included at `data/NOAA_0253923_original.xlsx`. If that file is
   missing, download it using the GitHub link above and save it at that path.
2. Read the workbook's `Data_Dictionary` sheet for variable definitions.
3. `eda.py` reads the `MP_data_S2N` worksheet. It never modifies the Excel file.

There are 124 observations, 10 source columns, and 12 coordinate locations.
The main repeated-site subset has 65 observations at 8 sites observed in
every month from January through April 2021. Small artificial data tables in
`testing.py` are for tests only, never for research results.

## Run

Use Python 3.10 or later. In a terminal in this folder:

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
python eda.py
python testing.py
python -m flake8 eda.py testing.py
```

On Windows, activate with `.venv\Scripts\activate` instead.

The analysis saves CSV summaries in `tables/` and PNG plots in `figures/`.
It runs a 9,999-resample SciPy test, so it can take a short time. The seed is
163. The reported conditional result was produced with SciPy 1.18.1.

## Files

- `eda.py`: data checks, summaries, plots, and exploratory SciPy analysis.
- `testing.py`: 16 plain assertion tests with small, hand-checkable examples.
- `requirements.txt`: Python packages needed to run the project.
- `data/NOAA_0253923_original.xlsx`: the original data copy read by the code.
- `eda.pdf`: numbered PDF report, named as required by Gradescope.
- `tables/`: supporting statistics, all category counts, and sensitivity checks.
- `figures/`: site plot, concentration boxplots, monthly lines, and rank heatmap.

## Methods and interpretation

Coordinates sorted by latitude and longitude define sites S01-S12.
The common-month sites are S01, S02, S03, S04, S05, S07, S10, and S12.
Site-month medians give each site one summary per month. They do not make
uncertainty equal: some cells contain only one observation.

Rank 1 means the highest observed concentration within a month. All ties at
the second rank are included in the top-two count, so a tied month could have
more than two selected sites. No ties occur between sites in the main table.

The Friedman statistic is calibrated with SciPy's paired-sample permutation
function. Month values stay grouped by site. This avoids relying on the
small-sample chi-square approximation but still assumes month-label
exchangeability and independent site blocks. Nearby sites and consecutive
months make these assumptions uncertain. Its p-value is conditional and
exploratory, not a proof of no month effect. Spearman coefficients describe
ranking agreement; their default p-values are not used.

Mean-versus-median and leave-one-month-out checks are saved. Removing the
largest observation in the common subset leaves S04 without a February
value. That missing cell is not filled with zero, and the original 8-by-4
test is not repeated on that incomplete table.

No claim is made about year-round seasonality, long-term trends, pollution
sources, ecological safety thresholds, or harm to animals or people.

## Resources and assistance

- Celine Wu's submitted project proposal.
- NOAA metadata, data dictionary, and original workbook.
- CSE 163 course materials and SciPy documentation.
- OpenAI Codex assisted with code, tests, interpretation checks, and drafting.
