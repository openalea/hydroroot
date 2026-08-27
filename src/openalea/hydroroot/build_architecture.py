"""
Build a root architecture in AQUA format from SmartRoot base/cut CSV exports.

2026-08-24: Refactored by codex from an original script from FB.
F. Bauget 2024-04-08: added possibility to add segments on LRII.
F. Bauget 2025-03-25: added comments.

SmartRoot exports are expected to contain one base file and one file per cut,
named with the same plant prefix: ``<plant>-base.csv``, ``<plant>-cut1.csv``,
...

AQUA team format is a tab-separated file with columns:
    distance_from_base_(mm)
    lateral_root_length_(mm)
    averaged_diameter_(mm)

Console diagnostics, per plant:
    cut 1 xxx
        Difference between the expected total length from CSV files and the
        built total length after cut 1. Ideally zero; a large value indicates
        unused or wrongly pasted roots.
    plant primary length: zzz total length: www lengths <>: vvv
        zzz is the primary root length, www the total CSV length, vvv the
        difference between CSV and built architecture lengths.
    plant cut lengths: [lPR, lcut1, lcut2, ...]
        Cut lengths in m, sorted from longest to shortest.
"""

from __future__ import annotations

import glob
import os
from dataclasses import dataclass, field
from typing import Dict, Iterable, List, Optional, Tuple

import numpy as np
import pandas as pd

CM_TO_MM = 10.0
MM_TO_M = 1.0e-3
CSV_SEPARATORS = (";", ",", "\t")


@dataclass
class Config:
    """User-editable settings for one reconstruction run."""

    path: str = "directory"
    list_fn: List[str] = field(default_factory=lambda: ["p15r1"])
    extension: str = "csv"
    primary_prefix: str = "PR"

    write_architecture: bool = True
    random_segment_order: bool = False
    sort_by_diameter: bool = False
    reverse_cut_filename_order: bool = True

    flag_verbose: bool = False # some prints

    # Tolerances around the cut position. tol1 selects laterals close enough to
    # a cut, tol2 rejects laterals that extend too far past the cut.
    tol1: float = 0.1
    tol2: float = 0.1

@dataclass
class Results:
    """output from reconstruction run."""
    d_archi: dict[str, pd.DataFrame] = field(default_factory=dict) # dict with DataFrame of the archi in Aqua format
    l_cut: dict[str, list[float]] = field(default_factory=dict) # list of the distance from the tip to the base at each cut
    plant_length: dict[str, float] = field(default_factory=dict) # total length of the plant root
    output_name: dict[str, str] = field(default_factory=dict) # name of the output archi

def drop_mislabeled_primary_records(df: pd.DataFrame, cut_index: int, flag_verbose: bool = False) -> pd.DataFrame:
    """
    Drop rows exported as parentless primary roots whose name identifies a lateral.

    SmartRoot files sometimes contain rows with parent_name == '-1' although the
    root name contains 'lat' or 'LRII'. Those rows are treated as export mistakes.
    """
    rows_to_drop = []
    for index, row in df.iterrows():
        root_name = str(row[" root_name"])
        if str(row[" parent_name"]) == "-1" and ("lat" in root_name or "lR" in root_name or "LRII" in root_name):
            rows_to_drop.append(index)

    if rows_to_drop:
        if flag_verbose: print(
            "cut nb ",
            cut_index,
            " ",
            "the following rows will be dropped from base because it is a primary "
            "root with 'lat' in its name: ",
        )
        for row_index in rows_to_drop:
            if flag_verbose: print(df[" root_name"][row_index])
        df = df.drop(rows_to_drop)

    return df.reset_index(drop=True)

def move_primary_row_to_end(df: pd.DataFrame, primary_prefix: str) -> pd.DataFrame:
    """Ensure the primary-root tip marker is the last row."""
    if df.empty or primary_prefix in str(df.iloc[-1][" root_name"]):
        return df

    primary_mask = df[" root_name"].str.contains(primary_prefix, na=False)
    return pd.concat([df[~primary_mask], df[primary_mask]], ignore_index=True)


def rename_segment_laterals(df: pd.DataFrame) -> pd.DataFrame:
    """
    Rename children of SmartRoot 'seg' objects as LRII-N.

    This keeps detached laterals consistent with the order naming expected by
    the reconstruction logic, and updates child parent names accordingly.
    """
    segment_child_mask = df[" parent_name"].str.contains("seg", case=False, na=False)
    for i, row_index in enumerate(df.index[segment_child_mask], start=1):
        old_root_id = df.at[row_index, " root"]
        new_root_name = f"LRII-{i}"

        df.at[row_index, " root_name"] = new_root_name
        df.loc[df[" parent"] == old_root_id, " parent_name"] = new_root_name

    return df


def read_smartroot_csv(filename: str) -> pd.DataFrame:
    """Read a SmartRoot CSV, normalize column names, and convert cm to mm."""
    df = None
    for separator in CSV_SEPARATORS:
        candidate = pd.read_csv(filename, sep=separator, keep_default_na=True)
        if candidate.shape[1] > 1:
            df = candidate
            break

    if df is None:
        raise ValueError(f"Unable to read SmartRoot CSV with known separators: {filename}")

    df = df.rename(
        columns={
            " insertion_position": "x",
            " root_order": "order",
            " length": "l",
            " diameter": "d",
        }
    )
    df = df.sort_values(by=["x"], ascending=True)
    df[" root_name"] = df[" root_name"].astype("string")
    df[" parent_name"] = df[" parent_name"].astype("string")
    df["order"] = df["order"].astype("string")
    df = rename_segment_laterals(df)
    df[["x", "l", "d"]] = df[["x", "l", "d"]] * CM_TO_MM
    return df


def find_plant_files(directory: str, prefix_fn: str, extension: str, reverse_order: bool) -> List[str]:
    """Return base and cut files for one plant, ordered in reconstruction order."""
    filenames = glob.glob(os.path.join(directory, f"{prefix_fn}-*.{extension}"))
    if not filenames:
        return []

    filenames.sort(reverse=reverse_order)
    if "base" in os.path.basename(filenames[-1]).lower():
        filenames = filenames[-1:] + filenames[:-1]
    return filenames


def load_plant_files(filenames: Iterable[str]) -> Dict[int, pd.DataFrame]:
    """Load one plant's ordered CSV files into {cut_index: dataframe}."""
    return {index: read_smartroot_csv(filename) for index, filename in enumerate(filenames)}


def build_base_dataframe(df: pd.DataFrame, primary_prefix: str) -> pd.DataFrame:
    """Create the primary-root dataframe from the base SmartRoot export."""
    primary_mask = df[" root_name"].str.contains(primary_prefix, case=False, na=False)
    primary_child_mask = df[" parent_name"].str.contains(primary_prefix, case=False, na=False)
    base = df[primary_mask | primary_child_mask].reset_index(drop=True)
    base = drop_mislabeled_primary_records(base, 0)

    base.loc[0, "x"] = base["l"].loc[0]
    base.loc[0, "l"] = 0.0
    base.loc[0, "order"] = "1"
    base = base.sort_values(by=["x"], ascending=True).reset_index(drop=True)
    base["Cut_nb"] = 0
    return base


def lrii_for_parent_positions(
    cut_lrii: pd.DataFrame,
    parent_positions: pd.DataFrame,
    first_order_index: int,
) -> pd.DataFrame:
    """
    Select LRII rows belonging to parent_positions and assign AQUA order labels.

    first_order_index is the 1-based order number assigned to the first row in
    parent_positions.
    """
    if cut_lrii.empty or " root" not in parent_positions.columns:
        return pd.DataFrame(columns=cut_lrii.columns)

    chunks = []
    for offset, (_, parent_row) in enumerate(parent_positions.iterrows()):
        root_id = parent_row[" root"]
        lrii = cut_lrii[cut_lrii[" parent"] == root_id].copy()
        if lrii.empty:
            continue
        lrii["order"] = f"1-{first_order_index + offset}"
        chunks.append(lrii)

    if not chunks:
        return pd.DataFrame(columns=cut_lrii.columns)
    return pd.concat(chunks, ignore_index=True)


def extract_cut_dataframes(
    df: pd.DataFrame,
    primary_prefix: str,
    cut_index: int,
    random_segment_order: bool,
) -> Optional[Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]]:
    """
    Extract primary-root data, detached segments, and LRII rows for one cut.

    Returns None when the cut contains no primary-root row for the selected plant.
    """
    if df[df[" root_name"].str.contains(primary_prefix, case=False, na=False)].empty:
        return None

    attached = df[~df[" root_name"].str.contains("seg", case=False, na=False)]
    primary_mask = attached[" root_name"].str.contains(primary_prefix, case=False, na=False)
    primary_child_mask = attached[" parent_name"].str.contains(primary_prefix, case=False, na=False)
    cut_pr = attached[primary_mask | primary_child_mask].reset_index(drop=True)
    cut_pr = drop_mislabeled_primary_records(cut_pr, cut_index)

    if cut_pr.empty:
        cut_pr = pd.DataFrame({"x": 0.0, "l": 0.0, "order": "1", "d": 0.0}, index=[0])
    else:
        for row_index in range(1, len(cut_pr)):
            if cut_pr.loc[row_index, "x"] > cut_pr["l"].loc[0]:
                cut_pr.loc[row_index, "x"] = cut_pr["l"].loc[0] - 1.0
        cut_pr.loc[0, "x"] = cut_pr["l"].loc[0]
        cut_pr.loc[0, "l"] = 0.0
        cut_pr.loc[0, "order"] = "1"

    cut_pr = cut_pr.sort_values(by=["x"], ascending=True)
    cut_pr = move_primary_row_to_end(cut_pr, primary_prefix)
    cut_pr["Cut_nb"] = cut_index

    parentless = df[df[" parent_name"].str.contains("-1", case=False, na=False)]
    cut_seg = parentless[parentless[" root_name"].str.contains("seg", case=False, na=False)].copy()
    cut_seg = cut_seg.sort_values("l", ascending=False)
    if not cut_seg.empty and random_segment_order:
        cut_seg = cut_seg.sample(len(cut_seg))
    cut_seg["x"] = cut_seg["l"]
    cut_seg = cut_seg.reset_index(drop=True)

    cut_lrii = df[df[" root_name"].str.contains("II", case=False, na=False)].reset_index(drop=True)
    cut_lrii = drop_mislabeled_primary_records(cut_lrii, cut_index)

    return cut_pr, cut_seg, cut_lrii


def total_length(df: pd.DataFrame) -> float:
    """Length contribution of a root dataframe in mm."""
    if df.empty:
        return 0.0
    return float(df["l"].sum())


def built_length(d_pr: pd.DataFrame, d_lrii: pd.DataFrame, d_lriii: pd.DataFrame) -> float:
    """Current total length of the reconstructed architecture in mm."""
    length = float(d_pr["x"].max() + d_pr["l"].sum())
    length += total_length(d_lrii)
    length += total_length(d_lriii)
    return length


def sort_laterals_and_segments(
    d_pr: pd.DataFrame,
    cut_seg: pd.DataFrame,
    cut_lrii: pd.DataFrame,
    config: Config,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Sort candidate laterals and detached segments before segment pasting."""
    if config.random_segment_order:
        return d_pr, cut_seg

    if config.sort_by_diameter:
        return (
            d_pr.sort_values(by=["d"], ascending=False, ignore_index=True),
            cut_seg.sort_values("d", ascending=False, ignore_index=True),
        )

    lrii_length_by_parent = cut_lrii.groupby(" parent")["l"].sum()

    d_pr = d_pr.copy()
    d_pr["total_lat_length"] = d_pr["l"] + d_pr[" root"].map(lrii_length_by_parent).fillna(0.0)
    d_pr = d_pr.sort_values("total_lat_length", ascending=False, ignore_index=True)

    cut_seg = cut_seg.copy()
    cut_seg["total_length"] = cut_seg["l"] + cut_seg[" root"].map(lrii_length_by_parent).fillna(0.0)
    cut_seg = cut_seg.sort_values("total_length", ascending=False, ignore_index=True)

    return d_pr, cut_seg


def weighted_average_diameter(first: pd.Series, second: pd.Series) -> float:
    """Average diameter weighted by segment lengths."""
    return float((first["l"] * first["d"] + second["l"] * second["d"]) / (first["l"] + second["l"]))


def attach_segments_to_primary_laterals(
    d_pr: pd.DataFrame,
    cut_seg: pd.DataFrame,
    cut_lrii: pd.DataFrame,
    lcut_previous: float,
    config: Config,
) -> Tuple[pd.DataFrame, pd.DataFrame, int]:
    """
    Paste detached segments onto first-order laterals that reached the previous cut.

    Returns the updated primary dataframe, LRII rows carried by pasted segments,
    and the number of consumed detached segments.
    """
    lrii_chunks = []
    segment_index = 0
    segment_count = len(cut_seg)

    for row_index in range(len(d_pr)):
        lateral_tip = d_pr.x[row_index] + d_pr.l[row_index]
        if lateral_tip > lcut_previous * (1.0 + config.tol2):
            continue
        if lateral_tip < lcut_previous * (1.0 - config.tol1) or d_pr.l[row_index] <= 0.0:
            continue
        if segment_index >= segment_count:
            continue

        segment = cut_seg.loc[segment_index]
        d_pr.loc[row_index, "d"] = weighted_average_diameter(d_pr.loc[row_index], segment)

        carried_lrii = cut_lrii[cut_lrii[" parent"] == segment[" root"]].copy()
        if not carried_lrii.empty:
            carried_lrii["order"] = f"1-{d_pr['order'][row_index]}"
            carried_lrii["x"] += d_pr.l[row_index]
            lrii_chunks.append(carried_lrii)

        d_pr.loc[row_index, "l"] = d_pr.l[row_index] + segment["l"]
        segment_index += 1

    d_lrii_current = (
        pd.concat(lrii_chunks, ignore_index=True)
        if lrii_chunks
        else pd.DataFrame(columns=cut_lrii.columns)
    )
    return d_pr, d_lrii_current, segment_index


def attach_segments_to_lrii(
    d_lrii: pd.DataFrame,
    cut_seg: pd.DataFrame,
    cut_lrii: pd.DataFrame,
    d_pr: pd.DataFrame,
    first_unused_segment: int,
    lcut_previous: float,
    config: Config,
) -> Tuple[pd.DataFrame, pd.DataFrame, int]:
    """Paste remaining detached segments onto LRII, producing LRIII rows."""
    lriii_chunks = []
    segment_index = first_unused_segment
    segment_count = len(cut_seg)

    for row_index in range(len(d_lrii)):
        if segment_index >= segment_count:
            break

        primary_lateral_index = int(str(d_lrii.order[row_index]).replace("1-", "")) - 1
        lateral_tip = d_pr.x[primary_lateral_index] + d_lrii.x[row_index] + d_lrii.l[row_index]
        if lateral_tip < lcut_previous * (1.0 - config.tol1) or d_lrii.l[row_index] <= 0.0:
            continue

        segment = cut_seg.loc[segment_index]
        d_lrii.loc[row_index, "d"] = weighted_average_diameter(d_lrii.loc[row_index], segment)

        carried_lriii = cut_lrii[cut_lrii[" parent"] == segment[" root"]].copy()
        if not carried_lriii.empty:
            carried_lriii["order"] = f"{d_lrii['order'][row_index]}-{row_index + 1}"
            carried_lriii["x"] += d_lrii.l[row_index]
            lriii_chunks.append(carried_lriii)

        d_lrii.loc[row_index, "l"] = d_lrii.l[row_index] + segment["l"]
        segment_index += 1

    d_lriii_current = (
        pd.concat(lriii_chunks, ignore_index=True)
        if lriii_chunks
        else pd.DataFrame(columns=cut_lrii.columns)
    )
    return d_lrii, d_lriii_current, segment_index


def merge_cut_into_primary(d_pr: pd.DataFrame, cut_pr: pd.DataFrame) -> pd.DataFrame:
    """
    Append a cut's attached primary-root part to the growing primary dataframe.

    If the current architecture ends with a zero-length primary tip marker, the
    marker is replaced by the cut data and its diameter is averaged across the
    junction.
    """
    cut_pr = cut_pr.copy()
    junction_diameter = None

    if d_pr.l.iloc[-1] <= 0.0:
        junction_diameter = (
            d_pr.x.iloc[-1] * d_pr.d.iloc[-1] + cut_pr.x.iloc[-1] * cut_pr.d.iloc[-1]
        ) / (d_pr.x.iloc[-1] + cut_pr.x.iloc[-1])
        cut_pr.x += d_pr.x.iloc[-1]
        d_pr = d_pr[:-1]

    d_pr = pd.concat([d_pr, cut_pr], ignore_index=True)
    if junction_diameter is not None:
        d_pr.iloc[-1, d_pr.columns.get_loc("d")] = junction_diameter
    return d_pr


def append_laterals(base: pd.DataFrame, extra: pd.DataFrame) -> pd.DataFrame:
    """Append laterals while preserving an empty base dataframe."""
    if extra.empty:
        return base
    if base.empty:
        return extra
    return pd.concat([base, extra], ignore_index=True)


def add_cut_lrii_to_architecture(
    d_lrii: pd.DataFrame,
    cut_lrii: pd.DataFrame,
    cut_pr: pd.DataFrame,
    first_order_index: int,
) -> pd.DataFrame:
    """Add LRII rows attached to the still-connected primary part of a cut."""
    return append_laterals(
        d_lrii,
        lrii_for_parent_positions(cut_lrii, cut_pr, first_order_index=first_order_index),
    )


def drop_redundant_zero_length_rows(d_pr: pd.DataFrame) -> pd.DataFrame:
    """Drop zero-length rows immediately followed by another row of the same order."""
    next_order = d_pr["order"].shift(-1)
    redundant = d_pr["l"].eq(0.0) & d_pr["order"].eq(next_order)
    return d_pr.loc[~redundant].copy()


def append_lrii_tail_markers(d_r: pd.DataFrame, d_lrii: pd.DataFrame) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Append terminal l=0 rows to each LRII branch."""
    if d_lrii.empty:
        return d_r, d_lrii

    d_lrii = d_lrii.sort_values(by=["order", "x"], ascending=[True, True], ignore_index=True)
    chunks = []
    for row_index in range(len(d_r)):
        lrii = d_lrii[d_lrii["order"].eq(f"1-{row_index + 1}")]
        if lrii.empty:
            continue

        tail = lrii.tail(1).copy()
        tail.loc[:, "x"] = d_r.loc[row_index, "l"]
        tail.loc[:, "l"] = 0.0
        chunks.extend([lrii, tail])

    if chunks:
        d_r = pd.concat([d_r, *chunks], ignore_index=True)
    return d_r, d_lrii


def append_lriii_tail_markers(
    d_r: pd.DataFrame,
    d_lrii: pd.DataFrame,
    d_lriii: pd.DataFrame,
) -> pd.DataFrame:
    """Append terminal l=0 rows to each LRIII branch."""
    if d_lriii.empty:
        return d_r

    d_lriii = d_lriii.sort_values(by=["order", "x"], ascending=[True, True], ignore_index=True)
    chunks = []
    previous_order = ""
    same_order_rank = 1

    for row_index in range(len(d_lrii)):
        if d_lrii["order"].loc[row_index] != previous_order:
            previous_order = d_lrii["order"].loc[row_index]
            same_order_rank = 1
        else:
            same_order_rank += 1

        lriii_order = f"{d_lrii['order'].loc[row_index]}-{same_order_rank}"
        lriii = d_lriii[d_lriii["order"] == lriii_order]
        if lriii.empty:
            continue

        tail = lriii.tail(1).copy()
        tail.loc[:, "x"] = d_r.loc[row_index, "l"]
        tail.loc[:, "l"] = 0.0
        chunks.extend([lriii, tail])

    if chunks:
        d_r = pd.concat([d_r, *chunks], ignore_index=True)
    return d_r


def to_aqua_format(d_pr: pd.DataFrame, d_lrii: pd.DataFrame, d_lriii: pd.DataFrame) -> pd.DataFrame:
    """Return the final AQUA-format architecture dataframe."""
    d_r = drop_redundant_zero_length_rows(d_pr)
    d_r, d_lrii = append_lrii_tail_markers(d_r, d_lrii)
    d_r = append_lriii_tail_markers(d_r, d_lrii, d_lriii)
    d_r = d_r.reset_index(drop=True)
    d_r = d_r[["x", "l", "order", "d"]]
    return d_r.rename(
        columns={
            "x": "distance_from_base_(mm)",
            "l": "lateral_root_length_(mm)",
            "d": "averaged_diameter_(mm)",
        }
    )


def extract_all_cuts(
    dfs: Dict[int, pd.DataFrame],
    config: Config,
    d_pr: pd.DataFrame,
    totals: Dict[str, float],
) -> Tuple[Dict[int, pd.DataFrame], Dict[int, pd.DataFrame], Dict[int, pd.DataFrame], List[int], List[float]]:
    """Extract all non-base cuts and update expected lengths."""
    cut_pr: Dict[int, pd.DataFrame] = {}
    cut_seg: Dict[int, pd.DataFrame] = {}
    cut_lrii: Dict[int, pd.DataFrame] = {}
    cuts_without_data = [0]
    lcut = [float(d_pr["x"].iloc[-1])]

    for cut_index in range(1, len(dfs)):
        extracted = extract_cut_dataframes(
            dfs[cut_index],
            config.primary_prefix,
            cut_index,
            config.random_segment_order,
        )
        if extracted is None:
            cuts_without_data.append(cut_index)
            continue

        cut_pr[cut_index], cut_seg[cut_index], cut_lrii[cut_index] = extracted
        lcut.append(cut_pr[cut_index]["x"].iloc[-1] + lcut[-1])

        totals["exp"] += cut_pr[cut_index]["x"].iloc[-1] + total_length(cut_pr[cut_index])
        totals["exp"] += total_length(cut_seg[cut_index])
        totals["exp"] += total_length(cut_lrii[cut_index])
        totals[f"cut-{cut_index}"] = totals["exp"]

    return cut_pr, cut_seg, cut_lrii, cuts_without_data, lcut


def build_plant_architecture(
    dfs: Dict[int, pd.DataFrame],
    config: Config,
    plant_name: str,
) -> Tuple[pd.DataFrame, Dict[str, float], List[float]]:
    """Reconstruct one plant architecture from its base and cut dataframes."""
    d_base = build_base_dataframe(dfs[0], config.primary_prefix)
    d_pr = d_base.copy(deep=True)

    cut_lrii: Dict[int, pd.DataFrame] = {
        0: dfs[0][dfs[0][" root_name"].str.contains("II", case=False, na=False)].reset_index(drop=True)
    }
    d_lrii = lrii_for_parent_positions(cut_lrii[0], d_base, first_order_index=1)
    d_lriii = pd.DataFrame(columns=[])

    totals = {
        "exp": float(d_pr["x"].iloc[-1] + d_pr["l"].sum() + cut_lrii[0]["l"].sum()),
    }
    totals["base"] = totals["exp"]

    cut_pr, cut_seg, later_cut_lrii, cuts_without_data, lcut = extract_all_cuts(
        dfs,
        config,
        d_pr,
        totals,
    )
    cut_lrii.update(later_cut_lrii)

    for cut_index in range(1, len(dfs)):
        if cut_index in cuts_without_data:
            continue

        free_segment_count = len(cut_seg[cut_index])
        d_pr["order"] = np.arange(1, len(d_pr) + 1).astype(str)

        d_pr, cut_seg[cut_index] = sort_laterals_and_segments(
            d_pr,
            cut_seg[cut_index],
            cut_lrii[cut_index],
            config,
        )

        d_pr, d_lrii_current, used_segment_count = attach_segments_to_primary_laterals(
            d_pr,
            cut_seg[cut_index],
            cut_lrii[cut_index],
            lcut[cut_index - 1],
            config,
        )

        d_lrii = d_lrii.reset_index(drop=True)
        d_lrii, d_lriii_current, used_segment_count = attach_segments_to_lrii(
            d_lrii,
            cut_seg[cut_index],
            cut_lrii[cut_index],
            d_pr,
            used_segment_count,
            lcut[cut_index - 1],
            config,
        )
        d_lriii = append_laterals(d_lriii, d_lriii_current)
        d_lrii = append_laterals(d_lrii, d_lrii_current)

        if used_segment_count < free_segment_count and config.flag_verbose:
            print("somme segments not used", "cut", cut_index, free_segment_count - used_segment_count)

        cut_seg[cut_index] = cut_seg[cut_index].sort_values("l", ascending=False, ignore_index=True)
        d_pr = d_pr.sort_values(by=["x"], ascending=True, ignore_index=True)
        d_pr = move_primary_row_to_end(d_pr, config.primary_prefix)

        d_lrii = add_cut_lrii_to_architecture(
            d_lrii,
            cut_lrii[cut_index],
            cut_pr[cut_index],
            first_order_index=len(d_pr),
        )
        d_pr = merge_cut_into_primary(d_pr, cut_pr[cut_index])

        totals[f"build-cut-{cut_index}"] = built_length(d_pr, d_lrii, d_lriii)
        if config.flag_verbose: print("cut ", cut_index, totals[f"cut-{cut_index}"] - totals[f"build-cut-{cut_index}"])

    d_pr["order"] = "1"
    totals["build"] = built_length(d_pr, d_lrii, d_lriii)

    lcut_m = sorted((cut_length * MM_TO_M for cut_length in lcut), reverse=True)
    d_r = to_aqua_format(d_pr.copy(), d_lrii, d_lriii)

    if config.flag_verbose:
        print(
        "plant: ", plant_name, "primary length: ", d_pr.x.max(), "total length:", totals["exp"], "lengths <> :", totals["exp"] - totals["build"],
        )
        print("cut lengths:", [float(cut_length) for cut_length in lcut_m])

    return d_r, totals, lcut_m


def write_architecture_file(
    d_r: pd.DataFrame,
    directory: str,
    prefix_fn: str,
    should_write: bool,
    flag_verbose: bool = False,
) -> Optional[str]:
    """Write, plot, and re-read the AQUA architecture file when requested."""
    if not should_write:
        return None

    output_name = os.path.join(directory, f"{prefix_fn}.txt")
    d_r.to_csv(output_name, sep="\t", index=False)

    return output_name


def process_plant(config: Optional[Config] = None, results: Optional[Results] = None):
    """Load, reconstruct, and optionally write one plant architecture."""
    config = config or Config()
    results = results or Results()
    directory = config.path

    for prefix_fn in config.list_fn:
        filenames = find_plant_files(
            directory,
            prefix_fn,
            config.extension,
            config.reverse_cut_filename_order,
        )
        if not filenames:
            print('no files to process for: ', prefix_fn)
            return None

        dfs = load_plant_files(filenames)
        results.d_archi[prefix_fn], results.plant_length[prefix_fn], results.l_cut[prefix_fn] = (
            build_plant_architecture(dfs, config, plant_name=prefix_fn))
        results.output_name[prefix_fn] = write_architecture_file(results.d_archi[prefix_fn], directory, prefix_fn,
                                                      config.write_architecture, flag_verbose = config.flag_verbose)

if __name__ == "__main__":
    process_plant()
