"""Explore NOAA 0253923 coastal microplastic measurements for CSE 163.

Run this file to save summary tables and figures beside the script.
The original Excel workbook is never changed.
"""

from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from scipy import stats


plt.switch_backend("Agg")
MONTHS = ["2021-01", "2021-02", "2021-03", "2021-04"]
BASE = Path(__file__).resolve().parent
DATA_FILE = BASE / "data" / "NOAA_0253923_original.xlsx"
COLUMN_NAMES = {
    "Date": "date",
    "Latitude": "latitude",
    "Longitude": "longitude",
    "Microplastics measurement": "concentration",
    "Unit": "unit",
    "Sampling Method": "sampling_method",
    "Short Reference": "short_reference",
    "Long Reference": "long_reference",
    "DOI": "doi",
    "ORGANIZATION": "organization"
}


def load_data(file_path):
    """Read the source worksheet and give its columns shorter names."""
    data = pd.read_excel(file_path, sheet_name="MP_data_S2N")
    data.columns = data.columns.str.strip()
    return data.rename(columns=COLUMN_NAMES)


def clean_data(data):
    """Check the records and add a site label and month to each row.

    Remove exact duplicates and rows with invalid required values.
    Keep real zero concentrations. The input table is not changed.
    """
    data = data.copy()
    data["date"] = pd.to_datetime(data["date"], errors="coerce")
    numeric = ["latitude", "longitude", "concentration"]
    for column in numeric:
        data[column] = pd.to_numeric(data[column], errors="coerce")
    data[numeric] = data[numeric].replace(
        [float("inf"), -float("inf")], float("nan")
    )
    data = data.drop_duplicates()
    data = data.dropna(subset=["date"] + numeric)
    valid_latitude = data["latitude"].between(-90, 90)
    valid_longitude = data["longitude"].between(-180, 180)
    valid_concentration = data["concentration"] >= 0
    data = data[valid_latitude & valid_longitude & valid_concentration].copy()
    sites = data[["latitude", "longitude"]].drop_duplicates()
    sites = sites.sort_values(["latitude", "longitude"])
    site_labels = []
    for number in range(1, len(sites) + 1):
        site_labels.append("S" + str(number).zfill(2))
    sites["site"] = site_labels
    data = data.merge(sites, on=["latitude", "longitude"], how="left")
    data["month"] = data["date"].dt.strftime("%Y-%m")
    return data


def select_common_data(data, months):
    """Keep requested months and sites with a record in every such month."""
    selected = data[data["month"].isin(months)]
    counts = selected.groupby("site")["month"].nunique()
    sites = counts[counts == len(months)].index
    return selected[selected["site"].isin(sites)].copy()


def monthly_concentrations(data, months, method="median"):
    """Calculate one median (or mean) for each site in each month."""
    grouped = data.groupby(["site", "month"])["concentration"].agg(method)
    # Put the months in columns so the same sites line up for comparisons.
    table = grouped.unstack("month")
    table = table.reindex(columns=months)
    return table.sort_index()


def rank_sites(table):
    """Rank each month's sites high to low; use the minimum rank for ties."""
    return table.rank(axis=0, ascending=False, method="min")


def hotspot_counts(table):
    """Count months in the top two ranks, including all ties at the cutoff."""
    ranks = rank_sites(table)
    in_top_two = ranks <= 2
    counts = in_top_two.sum(axis=1)
    return counts.rename("top_two_months")


def ranking_correlations(table):
    """Return descriptive Spearman coefficients without default p-values."""
    result = pd.DataFrame(index=table.columns, columns=table.columns,
                          dtype=float)
    for first in table.columns:
        for second in table.columns:
            coefficient = stats.spearmanr(table[first], table[second])
            result.loc[first, second] = coefficient.statistic
    return result


def friedman_statistic(january, february, march, april):
    """Calculate only the repeated-measures Friedman statistic."""
    return stats.friedmanchisquare(january, february, march, april).statistic


def month_test(table, resamples=9999):
    """Use SciPy to compare the four months at the same sites.

    Each array has one monthly median per site in identical site order.
    The test swaps month labels within sites. Its assumptions and limitations
    are discussed in the report, so the p-value is not a firm conclusion.
    """
    samples = []
    for month in MONTHS:
        samples.append(table[month].to_numpy())
    return stats.permutation_test(
        samples, friedman_statistic, permutation_type="samples",
        vectorized=False, n_resamples=resamples, alternative="greater",
        rng=163
    )


def save_data_checks(raw, data, common, table, folder):
    """Save missing-value counts and the main data checks."""
    raw.isna().sum().rename("missing").to_csv(folder / "missing_values.csv")
    checks = pd.Series({
        "source_rows": len(raw), "source_columns": len(raw.columns),
        "missing_cells": int(raw.isna().sum().sum()),
        "duplicate_rows": int(raw.duplicated().sum()),
        "duplicate_location_dates": int(raw.duplicated(
            ["date", "latitude", "longitude"]).sum()),
        "removed_rows": len(raw) - len(data),
        "sites": data["site"].nunique(), "common_rows": len(common),
        "common_sites": len(table), "common_months": len(table.columns),
        "missing_monthly_medians": int(table.isna().sum().sum())
    }, name="value")
    checks.to_csv(folder / "data_checks.csv")


def save_summaries(data, folder):
    """Save numerical summaries, category counts, and site information."""
    numeric = data[["concentration", "latitude", "longitude"]]
    numeric.describe().T.to_csv(folder / "numeric_summary.csv")
    categories = []
    category_columns = ["site", "month", "unit", "sampling_method",
                        "short_reference", "long_reference", "doi",
                        "organization"]
    for column in category_columns:
        for value, count in data[column].value_counts().sort_index().items():
            categories.append([column, value, count])
    category_counts = pd.DataFrame(
        categories, columns=["variable", "value", "count"]
    )
    category_counts.to_csv(folder / "category_counts.csv", index=False)
    grouped = data.groupby("site")
    sites = grouped[["latitude", "longitude"]].first()
    sites["samples"] = grouped["concentration"].count()
    sites["median"] = grouped["concentration"].median()
    sites.to_csv(folder / "site_summary.csv")
    coverage = pd.crosstab(data["site"], data["month"])
    coverage.to_csv(folder / "sampling_coverage.csv")
    grouped["concentration"].describe().to_csv(
        folder / "site_concentrations.csv"
    )


def save_monthly_results(common, table, folder):
    """Save the repeated-site records, monthly values, and rankings."""
    common.to_csv(folder / "common_month_records.csv", index=False)
    pd.crosstab(common["site"], common["month"]).to_csv(
        folder / "common_sampling_counts.csv"
    )
    table.to_csv(folder / "monthly_medians.csv")
    rank_sites(table).to_csv(folder / "monthly_ranks.csv")
    ranking_correlations(table).to_csv(folder / "rank_correlations.csv")


def save_sensitivity_results(common, table, folder):
    """Check results using means and by leaving out one month or high value."""
    sensitivity = hotspot_counts(table).to_frame("median_top_two_months")
    means = monthly_concentrations(common, MONTHS, "mean")
    means.to_csv(folder / "monthly_means.csv")
    sensitivity["mean_top_two_months"] = hotspot_counts(means)
    for month in MONTHS:
        sensitivity["omit_" + month] = hotspot_counts(
            table.drop(columns=month)
        )
    sensitivity.to_csv(folder / "hotspot_sensitivity.csv")
    largest = common["concentration"].idxmax()
    reduced = monthly_concentrations(common.drop(index=largest), MONTHS)
    reduced.to_csv(folder / "largest_sample_removed.csv")
    common.loc[[largest], ["site", "date", "concentration"]].to_csv(
        folder / "largest_common_sample.csv", index=False
    )


def save_tables(raw, data, common, table, folder):
    """Save all the tables used in the report."""
    save_data_checks(raw, data, common, table, folder)
    save_summaries(data, folder)
    save_monthly_results(common, table, folder)
    save_sensitivity_results(common, table, folder)


def plot_locations(data, folder):
    """Make a coordinate plot; points identify sites, not a coastline map."""
    sites = data.groupby("site")[["latitude", "longitude"]].first()
    fig, axes = plt.subplots(1, 2, figsize=(9, 4))
    axes[0].scatter(sites["longitude"], sites["latitude"],
                    color="#23778a", s=40)
    cluster = sites[sites["longitude"].between(2.18, 2.23)]
    axes[1].scatter(cluster["longitude"], cluster["latitude"],
                    color="#23778a", s=40)
    for site, row in sites.iterrows():
        if site in cluster.index:
            ax = axes[1]
        else:
            ax = axes[0]
        ax.annotate(site, (row["longitude"], row["latitude"]),
                    xytext=(5, 4), textcoords="offset points", fontsize=11)
    axes[0].annotate("S04-S08", (2.21, 41.40), xytext=(-50, 15),
                     textcoords="offset points", fontsize=11)
    for ax in axes:
        ax.set(xlabel="Longitude (degrees east)",
               ylabel="Latitude (degrees north)")
        ax.margins(0.18)
    axes[0].set_title("All 12 sampling locations")
    axes[1].set_title("Detail of the central cluster")
    fig.tight_layout()
    fig.savefig(folder / "sampling_locations.png", dpi=180)
    plt.close(fig)


def plot_concentrations(data, common, folder):
    """Compare raw concentrations and the equal-month median summaries."""
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.6), sharey=True)
    sns.boxplot(data=data, x="site", y="concentration", ax=axes[0],
                order=sorted(data["site"].unique()), color="#88b7c0")
    table = monthly_concentrations(common, MONTHS)
    # Make a long table so seaborn can draw a separate box for each site.
    points = table.reset_index().melt(
        id_vars="site", value_name="concentration"
    )
    sns.boxplot(data=points, x="site", y="concentration", ax=axes[1],
                color="#b6cda1")
    sns.stripplot(data=points, x="site", y="concentration", ax=axes[1],
                  color="#25343b", size=4, jitter=False)
    axes[0].set_title("All samples (12 sites)")
    axes[1].set_title("Monthly medians (8 sites)")
    for ax in axes:
        ax.set(xlabel="Site", ylabel="Microplastics (pieces/m³)")
        ax.set_yscale("log")
        ax.tick_params(axis="x", rotation=45)
    fig.suptitle("Concentration differs within and between sites")
    fig.tight_layout()
    fig.savefig(folder / "site_concentrations.png", dpi=180)
    plt.close(fig)


def plot_months(table, folder):
    """Plot each site's four monthly medians without pooling sites."""
    fig, ax = plt.subplots(figsize=(8, 4.5))
    for site in table.index:
        ax.plot(table.columns, table.loc[site], marker="o", label=site)
    ax.set(xlabel="Month", ylabel="Microplastics (pieces/m³)",
           title="Monthly median concentration at eight repeated sites")
    ax.legend(title="Site", bbox_to_anchor=(1.02, 1), loc="upper left")
    fig.tight_layout()
    fig.savefig(folder / "monthly_concentrations.png", dpi=180)
    plt.close(fig)


def plot_ranks(table, folder):
    """Show relative ranks; rank 1 is highest and is not a safety threshold."""
    fig, ax = plt.subplots(figsize=(7, 4.8))
    sns.heatmap(rank_sites(table), annot=True, fmt=".0f", cmap="YlGnBu_r",
                vmin=1, vmax=len(table), linewidths=0.5, ax=ax,
                cbar_kws={"label": "Rank (1 = highest)"})
    ax.set(xlabel="Month", ylabel="Site",
           title="High-concentration locations change across months")
    fig.tight_layout()
    fig.savefig(folder / "monthly_ranks.png", dpi=180)
    plt.close(fig)


def main():
    """Run the EDA on the supplied NOAA workbook and save all outputs."""
    tables = BASE / "tables"
    figures = BASE / "figures"
    tables.mkdir(exist_ok=True)
    figures.mkdir(exist_ok=True)
    raw = load_data(DATA_FILE)
    data = clean_data(raw)
    common = select_common_data(data, MONTHS)
    table = monthly_concentrations(common, MONTHS)
    if table.isna().any().any() or len(table) < 2:
        raise ValueError("Month comparisons need complete repeated sites.")
    save_tables(raw, data, common, table, tables)
    result = month_test(table)
    tests = pd.DataFrame({
        "measure": ["friedman_statistic", "conditional_permutation_p",
                    "permutations", "seed"],
        "value": [result.statistic, result.pvalue, 9999, 163]
    })
    tests.to_csv(tables / "exploratory_month_test.csv", index=False)
    sns.set_theme(style="whitegrid", context="notebook")
    plot_locations(data, figures)
    plot_concentrations(data, common, figures)
    plot_months(table, figures)
    plot_ranks(table, figures)
    print("Raw rows and columns:", raw.shape)
    print("Clean rows:", len(data), "Sites:", data["site"].nunique())
    print("Missing values:", int(raw.isna().sum().sum()))
    print("Common-month records:", len(common), "Table:", table.shape)
    print("Friedman statistic:", result.statistic)
    print("Conditional permutation p-value:", result.pvalue)
    print("Interpretation depends on assumptions; see the EDA report.")
    print("Tables and figures saved beside eda.py.")


if __name__ == "__main__":
    main()
