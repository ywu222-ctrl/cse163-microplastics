"""Test data-processing functions with small, hand-checkable examples."""

import pandas as pd

from eda import clean_data
from eda import friedman_statistic
from eda import hotspot_counts
from eda import monthly_concentrations
from eda import rank_sites
from eda import ranking_correlations
from eda import select_common_data


def make_data(rows):
    """Build a small raw table from date, latitude, longitude, and value."""
    data = pd.DataFrame(
        rows, columns=["date", "latitude", "longitude", "concentration"]
    )
    data["unit"] = "pieces/m3"
    data["sampling_method"] = "Manta net"
    return data


def monthly_example():
    """Return observations with known monthly medians and means."""
    return pd.DataFrame(
        [
            ["S01", "2021-01", 1],
            ["S01", "2021-01", 3],
            ["S01", "2021-01", 9],
            ["S01", "2021-02", 2],
            ["S01", "2021-02", 8],
            ["S02", "2021-01", 4],
            ["S02", "2021-01", 6],
            ["S02", "2021-02", 10],
        ],
        columns=["site", "month", "concentration"],
    )


def test_clean_missing_values_and_zero():
    """Keep real zeros, parse numbers, and remove missing required data."""
    data = make_data(
        [
            ["2021-01-01", 41, 2, 0],
            ["2021-01-02", 41, 2, "2.5"],
            ["2021-01-03", 41, 2, None],
            ["not-a-date", 41, 2, 1],
            ["2021-01-05", None, 2, 1],
            ["2021-01-06", 41, None, 1],
            [None, 41, 2, 1],
        ]
    )
    result = clean_data(data)
    assert len(result) == 2
    assert sorted(result["concentration"].tolist()) == [0, 2.5]
    assert pd.api.types.is_datetime64_any_dtype(result["date"])
    assert pd.api.types.is_numeric_dtype(result["concentration"])


def test_clean_invalid_numbers_and_coordinates():
    """Remove invalid text, infinities, negatives, and invalid coordinates."""
    data = make_data(
        [
            ["2021-01-01", 41, 2, 1],
            ["2021-01-02", "bad", 2, 1],
            ["2021-01-03", 41, "bad", 1],
            ["2021-01-04", 41, 2, "bad"],
            ["2021-01-05", 41, 2, float("inf")],
            ["2021-01-06", 41, 2, float("-inf")],
            ["2021-01-07", float("inf"), 2, 1],
            ["2021-01-08", 41, float("inf"), 1],
            ["2021-01-09", 41, 2, -1],
            ["2021-01-10", 91, 2, 1],
            ["2021-01-11", -91, 2, 1],
            ["2021-01-12", 41, 181, 1],
            ["2021-01-13", 41, -181, 1],
        ]
    )
    result = clean_data(data)
    assert len(result) == 1
    assert result.iloc[0]["concentration"] == 1
    assert result.iloc[0]["date"] == pd.Timestamp("2021-01-01")


def test_clean_coordinate_boundaries():
    """Latitude and longitude boundary values are valid."""
    data = make_data(
        [
            ["2021-01-01", -90, -180, 1],
            ["2021-01-02", 90, 180, 1],
        ]
    )
    assert len(clean_data(data)) == 2


def test_clean_duplicate_rows():
    """Remove an exact duplicate without losing a different observation."""
    data = make_data(
        [
            ["2021-01-01", 41, 2, 3],
            ["2021-01-01", 41, 2, 3],
            ["2021-01-02", 41, 2, 4],
        ]
    )
    result = clean_data(data)
    assert len(result) == 2
    assert sorted(result["concentration"].tolist()) == [3, 4]


def test_clean_does_not_change_input():
    """Cleaning returns a new table and leaves the original untouched."""
    data = make_data(
        [
            ["2021-01-01", "41", "2", "3"],
            ["bad date", 41, 2, None],
        ]
    )
    original = data.copy(deep=True)
    clean_data(data)
    pd.testing.assert_frame_equal(data, original)


def test_site_numbers_are_deterministic():
    """Number sites by coordinates, not the order of input rows."""
    data = make_data(
        [
            ["2021-01-01", 42, 1, 1],
            ["2021-01-02", 41, 3, 1],
            ["2021-01-03", 41, 2, 1],
        ]
    )
    expected = {(41, 2): "S01", (41, 3): "S02", (42, 1): "S03"}
    result = clean_data(data)
    site_names = result.groupby(["latitude", "longitude"])["site"].first()
    assert site_names.to_dict() == expected
    reordered = clean_data(data.sort_values("date", ascending=False))
    other_names = reordered.groupby(["latitude", "longitude"])["site"].first()
    assert other_names.to_dict() == expected


def test_month_date_boundaries():
    """Assign dates at year and month boundaries to the correct month."""
    dates = [
        "2020-12-31", "2021-01-01", "2021-01-31",
        "2021-02-01", "2021-04-30", "2021-05-01",
    ]
    data = make_data([[date, 41, 2, 1] for date in dates])
    result = clean_data(data)
    expected = ["2020-12", "2021-01", "2021-01",
                "2021-02", "2021-04", "2021-05"]
    assert result["month"].tolist() == expected


def test_select_common_sites_and_months():
    """Keep all selected-month rows only for sites with complete coverage."""
    data = pd.DataFrame(
        [
            ["S01", "2020-12", 9],
            ["S01", "2021-01", 1],
            ["S01", "2021-01", 3],
            ["S01", "2021-02", 2],
            ["S02", "2021-01", 4],
            ["S03", "2021-01", 5],
            ["S03", "2021-02", 6],
            ["S03", "2021-03", 8],
        ],
        columns=["site", "month", "concentration"],
    )
    result = select_common_data(data, ["2021-01", "2021-02"])
    expected = data.iloc[[1, 2, 3, 5, 6]]
    pd.testing.assert_frame_equal(
        result.reset_index(drop=True), expected.reset_index(drop=True)
    )


def test_select_common_no_complete_site():
    """Return no observations when a requested month is absent everywhere."""
    data = monthly_example()
    result = select_common_data(data, ["2021-01", "2021-03"])
    assert result.empty


def test_monthly_medians():
    """Check medians for odd and even sample sizes."""
    result = monthly_concentrations(
        monthly_example(), ["2021-01", "2021-02"]
    )
    expected = pd.DataFrame(
        {"2021-01": [3, 5], "2021-02": [5, 10]},
        index=["S01", "S02"],
    )
    pd.testing.assert_frame_equal(
        result, expected, check_dtype=False, check_names=False
    )


def test_monthly_means():
    """Check the mean sensitivity calculation against hand calculations."""
    result = monthly_concentrations(
        monthly_example(), ["2021-01", "2021-02"], method="mean"
    )
    expected = pd.DataFrame(
        {"2021-01": [13 / 3, 5], "2021-02": [5, 10]},
        index=["S01", "S02"],
    )
    pd.testing.assert_frame_equal(
        result, expected, check_dtype=False, check_names=False
    )


def test_monthly_incomplete_sites_and_column_order():
    """Keep incomplete sites and follow the requested month order."""
    data = pd.DataFrame(
        [
            ["S01", "2021-01", 2],
            ["S01", "2021-02", 4],
            ["S02", "2021-01", 8],
        ],
        columns=["site", "month", "concentration"],
    )
    result = monthly_concentrations(data, ["2021-02", "2021-01"])
    expected = pd.DataFrame(
        {"2021-02": [4, float("nan")], "2021-01": [2, 8]},
        index=["S01", "S02"],
    )
    pd.testing.assert_frame_equal(
        result, expected, check_dtype=False, check_names=False
    )


def test_rank_sites_with_ties():
    """Give higher concentrations better ranks and use minimum tie ranks."""
    table = pd.DataFrame(
        {"2021-01": [10, 8, 8, 2], "2021-02": [4, 4, 2, 1]},
        index=["S01", "S02", "S03", "S04"],
    )
    expected = pd.DataFrame(
        {"2021-01": [1, 2, 2, 4], "2021-02": [1, 1, 3, 4]},
        index=table.index,
    )
    pd.testing.assert_frame_equal(
        rank_sites(table), expected, check_dtype=False
    )


def test_hotspot_counts_include_tied_second_place():
    """Include every site tied at rank two instead of breaking ties."""
    table = pd.DataFrame(
        {"2021-01": [10, 8, 8, 2], "2021-02": [4, 4, 2, 1]},
        index=["S01", "S02", "S03", "S04"],
    )
    expected = pd.Series(
        [2, 2, 1, 0], index=table.index, name="top_two_months"
    )
    pd.testing.assert_series_equal(hotspot_counts(table), expected)


def test_friedman_statistic_known_order():
    """Check the statistic for three sites with identical month ordering."""
    result = friedman_statistic(
        [1, 1, 1], [2, 2, 2], [3, 3, 3], [4, 4, 4]
    )
    assert result == 9


def test_ranking_correlations_known_orders():
    """Check perfectly matching and reversed rankings without p-values."""
    table = pd.DataFrame(
        {
            "increasing": [1, 2, 3, 4],
            "also_increasing": [10, 20, 30, 40],
            "decreasing": [40, 30, 20, 10],
        },
        index=["S01", "S02", "S03", "S04"],
    )
    result = ranking_correlations(table)
    assert abs(result.loc["increasing", "also_increasing"] - 1) < 1e-12
    assert abs(result.loc["increasing", "decreasing"] + 1) < 1e-12


def main():
    """Run every test and print a message if all assertions pass."""
    test_clean_missing_values_and_zero()
    test_clean_invalid_numbers_and_coordinates()
    test_clean_coordinate_boundaries()
    test_clean_duplicate_rows()
    test_clean_does_not_change_input()
    test_site_numbers_are_deterministic()
    test_month_date_boundaries()
    test_select_common_sites_and_months()
    test_select_common_no_complete_site()
    test_monthly_medians()
    test_monthly_means()
    test_monthly_incomplete_sites_and_column_order()
    test_rank_sites_with_ties()
    test_hotspot_counts_include_tied_second_place()
    test_friedman_statistic_known_order()
    test_ranking_correlations_known_orders()
    print("All 16 tests passed.")


if __name__ == "__main__":
    main()
