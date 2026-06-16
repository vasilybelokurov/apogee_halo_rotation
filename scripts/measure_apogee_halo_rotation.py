#!/usr/bin/env python3
"""Measure APOGEE DR17 in-situ/accreted halo rotation versus Galactocentric radius."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import zipfile

import numpy as np
import pandas as pd
from astropy.table import Table

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D


ROOT = Path(__file__).resolve().parents[1]
DATA_DIR = ROOT / "data" / "external"
ALLSTAR_URL = "https://data.sdss.org/sas/dr17/apogee/spectro/aspcap/dr17/synspec_rev1/allStarLite-dr17-synspec_rev1.fits"
ASTRONN_URL = "https://data.sdss.org/sas/dr17/apogee/vac/apogee-astronn/apogee_astroNN-DR17.fits"
VASILIEV_GC_URL = "https://zenodo.org/records/4891252/files/clusters.zip?download=1"
DEFAULT_ALLSTAR = DATA_DIR / "allStarLite-dr17-synspec_rev1.fits"
DEFAULT_ASTRONN = DATA_DIR / "apogee_astroNN-DR17.fits"
DEFAULT_VASILIEV_GC_MEMBERS = DATA_DIR / "vasiliev_baumgardt_2021_clusters.zip"
DEFAULT_OUT_CSV = ROOT / "products" / "apogee_halo_rotation_vs_radius.csv"
DEFAULT_OUT_JSON = ROOT / "products" / "apogee_halo_rotation_summary.json"
DEFAULT_OUT_PLOT = ROOT / "plots" / "apogee_halo_rotation_vs_radius.png"
DEFAULT_OUT_CHEM_PLOT = ROOT / "plots" / "apogee_halo_chemical_selection.png"
DEFAULT_OUT_RZ_PLOT = ROOT / "plots" / "apogee_halo_rz_z_distribution.png"
DEFAULT_OUT_ZSPLIT_PLOT = ROOT / "plots" / "apogee_halo_rotation_z_slices.png"
DEFAULT_OUT_RHIST_PLOT = ROOT / "plots" / "apogee_halo_radius_histograms_z_slices.png"
VELOCITY_YLIM = (-50.0, 125.0)

BAD_ASPCAPFLAGS = (
    "STAR_BAD",
    "TEFF_BAD",
    "LOGG_BAD",
    "VERY_BRIGHT_NEIGHBOR",
    "LOW_SNR",
    "PERSIST_HIGH",
    "PERSIST_JUMP_POS",
    "PERSIST_JUMP_NEG",
    "SUSPECT_RV_COMBINATION",
)
ABUNDANCE_ERR_COLUMNS = (
    "C_FE_ERR",
    "N_FE_ERR",
    "O_FE_ERR",
    "MG_FE_ERR",
    "AL_FE_ERR",
    "SI_FE_ERR",
    "MN_FE_ERR",
    "FE_H_ERR",
    "NI_FE_ERR",
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Build an APOGEE DR17 in-situ/accreted halo sample and measure median "
            "azimuthal velocity versus Galactocentric radius."
        )
    )
    parser.add_argument("--allstar", type=Path, default=DEFAULT_ALLSTAR)
    parser.add_argument("--astronn", type=Path, default=DEFAULT_ASTRONN)
    parser.add_argument("--vasiliev-gc-members", type=Path, default=DEFAULT_VASILIEV_GC_MEMBERS)
    parser.add_argument(
        "--vasiliev-prob-min",
        type=float,
        default=0.8,
        help="Remove Vasiliev GC candidates with membership probability greater than this value.",
    )
    parser.add_argument("--out-csv", type=Path, default=DEFAULT_OUT_CSV)
    parser.add_argument("--out-json", type=Path, default=DEFAULT_OUT_JSON)
    parser.add_argument("--out-plot", type=Path, default=DEFAULT_OUT_PLOT)
    parser.add_argument("--out-pdf", type=Path, default=None)
    parser.add_argument("--out-chem-plot", type=Path, default=DEFAULT_OUT_CHEM_PLOT)
    parser.add_argument("--out-chem-pdf", type=Path, default=None)
    parser.add_argument("--out-rz-plot", type=Path, default=DEFAULT_OUT_RZ_PLOT)
    parser.add_argument("--out-rz-pdf", type=Path, default=None)
    parser.add_argument("--out-zsplit-plot", type=Path, default=DEFAULT_OUT_ZSPLIT_PLOT)
    parser.add_argument("--out-zsplit-pdf", type=Path, default=None)
    parser.add_argument("--out-rhist-plot", type=Path, default=DEFAULT_OUT_RHIST_PLOT)
    parser.add_argument("--out-rhist-pdf", type=Path, default=None)
    parser.add_argument(
        "--zsplit-comparison-feh-max",
        type=float,
        nargs="+",
        default=[-1.3, -1.5],
        help=(
            "Overlay comparison curves in the z-slice plots for stars below these [Fe/H] values. "
            "Defaults: -1.3 dashed and -1.5 dotted. Use nan to disable."
        ),
    )
    parser.add_argument("--abundance-err-max", type=float, default=0.2)
    parser.add_argument("--velocity-err-max", type=float, default=50.0)
    parser.add_argument("--distance-max-kpc", type=float, default=15.0)
    parser.add_argument(
        "--measurement-feh-max",
        type=float,
        default=-1.0,
        help="Metallicity ceiling for the halo rotation measurement. Use nan to keep all BK-classified stars.",
    )
    parser.add_argument(
        "--radius-bins",
        type=float,
        nargs="+",
        default=None,
        help=(
            "Spherical Galactocentric radius bin edges in kpc. "
            "Default: automatic bins extending to the largest measured radius."
        ),
    )
    parser.add_argument("--min-count", type=int, default=3)
    parser.add_argument("--bootstrap", type=int, default=500)
    parser.add_argument("--random-seed", type=int, default=42)
    return parser.parse_args()


def _decode_string_array(values: np.ndarray) -> np.ndarray:
    out = values.astype(str)
    return np.char.strip(out)


def resolve_data_path(path: Path, url: str) -> Path:
    expanded = path.expanduser()
    if expanded.exists():
        return expanded
    raise FileNotFoundError(
        f"Missing data file: {expanded}. Download from {url} or run scripts/download_data.py."
    )


def read_allstar(path: Path) -> pd.DataFrame:
    columns = [
        "APOGEE_ID",
        "GAIAEDR3_SOURCE_ID",
        "PROGRAMNAME",
        "ASPCAPFLAG",
        "ASPCAPFLAGS",
        "EXTRATARG",
        "TEFF",
        "LOGG",
        "FE_H",
        "MG_FE",
        "AL_FE",
        *ABUNDANCE_ERR_COLUMNS,
    ]
    table = Table.read(path, format="fits", hdu=1, memmap=True)[columns]
    gaia_source_id = np.asarray(table["GAIAEDR3_SOURCE_ID"]).astype("int64")
    df = table.to_pandas()
    df["APOGEE_ID"] = _decode_string_array(df["APOGEE_ID"].to_numpy())
    df["GAIAEDR3_SOURCE_ID"] = gaia_source_id
    df["PROGRAMNAME"] = _decode_string_array(df["PROGRAMNAME"].to_numpy())
    df["ASPCAPFLAGS"] = _decode_string_array(df["ASPCAPFLAGS"].to_numpy())
    return df.drop_duplicates(subset="APOGEE_ID", keep="first")


def read_astronn(path: Path) -> pd.DataFrame:
    columns = [
        "APOGEE_ID",
        "dist",
        "galr",
        "galz",
        "galvr",
        "galvt",
        "galvz",
        "galvr_err",
        "galvt_err",
        "galvz_err",
        "Lz",
        "e",
        "rap",
        "zmax",
    ]
    table = Table.read(path, format="fits", hdu=1, memmap=True)[columns]
    df = table.to_pandas()
    df["APOGEE_ID"] = _decode_string_array(df["APOGEE_ID"].to_numpy())
    return df.drop_duplicates(subset="APOGEE_ID", keep="first")


def contains_any(series: pd.Series, substrings: tuple[str, ...]) -> pd.Series:
    mask = pd.Series(False, index=series.index)
    text = series.fillna("").astype(str)
    for substring in substrings:
        mask |= text.str.contains(substring, case=False, regex=False)
    return mask


def finite_all(df: pd.DataFrame, columns: list[str] | tuple[str, ...]) -> pd.Series:
    mask = pd.Series(True, index=df.index)
    for column in columns:
        mask &= np.isfinite(pd.to_numeric(df[column], errors="coerce"))
    return mask


def read_vasiliev_gc_source_ids(path: Path, prob_min: float) -> np.ndarray:
    if path.suffix.lower() == ".zip":
        return read_vasiliev_gc_source_ids_from_zip(path, prob_min=prob_min)

    table = Table.read(path, format="fits", hdu=1, memmap=True)
    name_by_lower = {name.lower(): name for name in table.colnames}
    source_col = name_by_lower.get("source_id")
    prob_col = name_by_lower.get("prob")
    if source_col is None or prob_col is None:
        raise KeyError(f"{path} must contain SOURCE_ID and prob columns")

    source_id = np.asarray(table[source_col]).astype("int64")
    prob = np.asarray(table[prob_col]).astype(float)
    keep = np.isfinite(prob) & (prob > prob_min) & (source_id > 0)
    return np.unique(source_id[keep])


def read_vasiliev_gc_source_ids_from_zip(path: Path, prob_min: float) -> np.ndarray:
    source_ids: list[np.ndarray] = []
    with zipfile.ZipFile(path) as archive:
        names = sorted(
            name
            for name in archive.namelist()
            if name.startswith("catalogues/") and name.endswith(".txt")
        )
        if not names:
            raise ValueError(f"No catalogues/*.txt files found in {path}")
        for name in names:
            with archive.open(name) as handle:
                table = pd.read_csv(
                    handle,
                    sep=r"\s+",
                    comment="#",
                    header=None,
                    usecols=[0, 16],
                    names=["source_id", "prob"],
                    dtype={"source_id": "int64", "prob": "float64"},
                )
            source_id = table["source_id"].to_numpy(dtype="int64")
            prob = table["prob"].to_numpy(dtype="float64")
            keep = np.isfinite(prob) & (source_id > 0) & (prob > prob_min)
            source_ids.append(source_id[keep])

    if not source_ids:
        return np.array([], dtype="int64")
    return np.unique(np.concatenate(source_ids))


def apply_bk22_clean_selection(
    df: pd.DataFrame,
    vasiliev_gc_source_ids: np.ndarray,
    abundance_err_max: float,
    velocity_err_max: float,
    distance_max_kpc: float,
) -> tuple[pd.DataFrame, dict[str, int]]:
    work = df.copy()
    # AstroNN DR17 stores heliocentric distance in parsec.
    work["dist_kpc"] = pd.to_numeric(work["dist"], errors="coerce") / 1000.0
    work["r_gc_kpc"] = np.sqrt(
        np.square(pd.to_numeric(work["galr"], errors="coerce"))
        + np.square(pd.to_numeric(work["galz"], errors="coerce"))
    )
    gaia_source_id = pd.to_numeric(work["GAIAEDR3_SOURCE_ID"], errors="coerce").fillna(0).astype("int64")

    required = [
        "GAIAEDR3_SOURCE_ID",
        "TEFF",
        "LOGG",
        "FE_H",
        "MG_FE",
        "AL_FE",
        "dist_kpc",
        "galr",
        "galz",
        "galvt",
        "r_gc_kpc",
        "Lz",
        *ABUNDANCE_ERR_COLUMNS,
        "galvr_err",
        "galvt_err",
        "galvz_err",
    ]
    masks: dict[str, pd.Series] = {
        "finite_required": finite_all(work, required),
        "giants_logg_lt_3": pd.to_numeric(work["LOGG"], errors="coerce") < 3.0,
        "velocity_errors_lt_limit": (
            (pd.to_numeric(work["galvr_err"], errors="coerce") < velocity_err_max)
            & (pd.to_numeric(work["galvt_err"], errors="coerce") < velocity_err_max)
            & (pd.to_numeric(work["galvz_err"], errors="coerce") < velocity_err_max)
        ),
        "abundance_errors_lt_limit": pd.Series(True, index=work.index),
        "bad_aspcapflags_removed": ~contains_any(work["ASPCAPFLAGS"], BAD_ASPCAPFLAGS),
        "aspcapflag_star_bad_bit_removed": (
            work["ASPCAPFLAG"].fillna(0).astype("int64") & (2**23)
        )
        == 0,
        "duplicate_extratarg_bit_removed": (
            work["EXTRATARG"].fillna(0).astype("int64") & (2**4)
        )
        == 0,
        "telluric_extratarg_bit_removed": (
            work["EXTRATARG"].fillna(0).astype("int64") & (2**2)
        )
        == 0,
        "magcloud_program_removed": ~work["PROGRAMNAME"]
        .fillna("")
        .astype(str)
        .str.contains("magcloud", case=False, regex=False),
        "gaia_source_id_present": gaia_source_id > 0,
        "vasiliev_gc_members_removed": ~gaia_source_id.isin(vasiliev_gc_source_ids),
        "distance_lt_limit": work["dist_kpc"] < distance_max_kpc,
    }
    for column in ABUNDANCE_ERR_COLUMNS:
        masks["abundance_errors_lt_limit"] &= (
            pd.to_numeric(work[column], errors="coerce") < abundance_err_max
        )

    combined = pd.Series(True, index=work.index)
    counts: dict[str, int] = {"input_merged": int(len(work))}
    for name, mask in masks.items():
        combined &= mask
        counts[name] = int(mask.sum())
        counts[f"cumulative_after_{name}"] = int(combined.sum())

    selected = work.loc[combined].copy()
    selected["population"] = np.where(
        (selected["FE_H"] > -0.4) | (selected["AL_FE"] >= -0.075),
        "in_situ",
        "accreted",
    )
    selected.loc[
        (selected["FE_H"] > -0.4) & (selected["AL_FE"] < -0.075),
        "population",
    ] = "in_situ"
    return selected, counts


def bootstrap_interval(values: np.ndarray, rng: np.random.Generator, n_boot: int) -> tuple[float, float]:
    if len(values) < 2 or n_boot <= 0:
        return np.nan, np.nan
    draws = rng.choice(values, size=(n_boot, len(values)), replace=True)
    med = np.nanmedian(draws, axis=1)
    return tuple(np.nanpercentile(med, [16.0, 84.0]))


def automatic_radius_bins(df: pd.DataFrame) -> list[float]:
    finite_radius = pd.to_numeric(df["r_gc_kpc"], errors="coerce").to_numpy(dtype=float)
    finite_radius = finite_radius[np.isfinite(finite_radius) & (finite_radius >= 0.0)]
    if len(finite_radius) == 0:
        return [0.0, 3.0, 5.0, 7.0, 9.0, 11.0, 13.0, 15.0]

    max_radius = float(np.nanmax(finite_radius))
    inner_edges = [0.0, 3.0, 5.0, 7.0, 9.0, 11.0, 13.0, 15.0]
    if max_radius <= inner_edges[-1]:
        edges = [edge for edge in inner_edges if edge < max_radius]
        edges.append(np.ceil(max_radius * 10.0) / 10.0 + 1.0e-6)
        return sorted(set(edges))

    outer_stop = np.ceil(max_radius / 5.0) * 5.0
    outer_edges = list(np.arange(20.0, outer_stop + 5.0, 5.0))
    edges = inner_edges + [edge for edge in outer_edges if edge > inner_edges[-1]]
    if edges[-1] <= max_radius:
        edges.append(max_radius + 1.0e-6)
    return edges


def summarize_rotation(
    df: pd.DataFrame,
    radius_bins: list[float],
    min_count: int,
    n_boot: int,
    random_seed: int,
) -> pd.DataFrame:
    rng = np.random.default_rng(random_seed)
    rows: list[dict[str, float | int | str]] = []
    edges = np.asarray(radius_bins, dtype=float)
    if np.any(np.diff(edges) <= 0):
        raise ValueError("--radius-bins must be strictly increasing")

    for population in ("in_situ", "accreted"):
        pop = df.loc[df["population"] == population]
        for lo, hi in zip(edges[:-1], edges[1:]):
            sub = pop.loc[(pop["r_gc_kpc"] >= lo) & (pop["r_gc_kpc"] < hi)]
            n = len(sub)
            values = pd.to_numeric(sub["galvt"], errors="coerce").to_numpy(dtype=float)
            values = values[np.isfinite(values)]
            if len(values) >= min_count:
                boot_lo, boot_hi = bootstrap_interval(values, rng, n_boot)
                rows.append(
                    {
                        "population": population,
                        "r_min_kpc": lo,
                        "r_max_kpc": hi,
                        "r_mid_kpc": 0.5 * (lo + hi),
                        "n": int(len(values)),
                        "median_vphi_kms": float(np.nanmedian(values)),
                        "mean_vphi_kms": float(np.nanmean(values)),
                        "std_vphi_kms": float(np.nanstd(values, ddof=1)) if len(values) > 1 else np.nan,
                        "median_vphi_boot16_kms": float(boot_lo),
                        "median_vphi_boot84_kms": float(boot_hi),
                        "median_feh": float(np.nanmedian(sub["FE_H"])),
                        "median_alfe": float(np.nanmedian(sub["AL_FE"])),
                    }
                )
            else:
                rows.append(
                    {
                        "population": population,
                        "r_min_kpc": lo,
                        "r_max_kpc": hi,
                        "r_mid_kpc": 0.5 * (lo + hi),
                        "n": int(len(values)),
                        "median_vphi_kms": np.nan,
                        "mean_vphi_kms": np.nan,
                        "std_vphi_kms": np.nan,
                        "median_vphi_boot16_kms": np.nan,
                        "median_vphi_boot84_kms": np.nan,
                        "median_feh": np.nan,
                        "median_alfe": np.nan,
                    }
                )
    return pd.DataFrame(rows)


def plot_rotation(binned: pd.DataFrame, out_plot: Path, out_pdf: Path | None) -> None:
    colors = {"in_situ": "#b33a3a", "accreted": "#2f6fb0"}
    labels = {"in_situ": "In-situ halo", "accreted": "Accreted halo"}

    fig, ax = plt.subplots(figsize=(7.0, 4.8))
    for population in ("in_situ", "accreted"):
        sub = binned.loc[binned["population"] == population].copy()
        ok = np.isfinite(sub["median_vphi_kms"])
        sub = sub.loc[ok]
        if sub.empty:
            continue
        low = sub["median_vphi_kms"] - sub["median_vphi_boot16_kms"]
        high = sub["median_vphi_boot84_kms"] - sub["median_vphi_kms"]
        low = low.fillna(0.0)
        high = high.fillna(0.0)
        yerr = np.vstack(
            [
                low,
                high,
            ]
        )
        ax.errorbar(
            sub["r_mid_kpc"],
            sub["median_vphi_kms"],
            yerr=yerr,
            fmt="o-",
            lw=1.8,
            ms=5.0,
            capsize=3.0,
            color=colors[population],
            label=labels[population],
        )

    ax.axhline(0.0, color="0.25", lw=1.0, ls="--", zorder=0)
    ax.set_xlabel(r"Galactocentric spherical radius $r_{\rm GC}$ [kpc]")
    ax.set_ylabel(r"Median azimuthal velocity $V_\phi$ [km s$^{-1}$]")
    ax.legend(frameon=False)
    ax.grid(True, color="0.9", lw=0.8)
    if not binned.empty:
        ax.set_xlim(float(binned["r_min_kpc"].min()), float(binned["r_max_kpc"].max()))
    else:
        ax.set_xlim(left=0)
    ax.set_ylim(*VELOCITY_YLIM)
    fig.tight_layout()

    out_plot.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_plot, dpi=220)
    if out_pdf is not None:
        out_pdf.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(out_pdf)
    plt.close(fig)


def plot_rotation_z_slices(
    measurement: pd.DataFrame,
    radius_bins: list[float],
    min_count: int,
    n_boot: int,
    random_seed: int,
    comparison_feh_max: list[float],
    out_plot: Path,
    out_pdf: Path | None,
) -> dict[str, dict[str, int]]:
    colors = {"in_situ": "#b33a3a", "accreted": "#2f6fb0"}
    labels = {"in_situ": "In-situ halo", "accreted": "Accreted halo"}
    comparison_styles = ("--", ":")
    comparison_thresholds = [value for value in comparison_feh_max if np.isfinite(value)]
    panels = [
        ("all", "All", measurement),
        ("abs_z_lt_3", r"$|z| < 3$ kpc", measurement.loc[np.abs(measurement["galz"]) < 3.0]),
        ("abs_z_gt_3", r"$|z| > 3$ kpc", measurement.loc[np.abs(measurement["galz"]) > 3.0]),
    ]

    binned_by_panel: dict[str, pd.DataFrame] = {}
    comparison_binned_by_panel: dict[tuple[str, float], pd.DataFrame] = {}
    counts: dict[str, dict[str, int]] = {}
    for offset, (key, _title, subset) in enumerate(panels):
        binned = summarize_rotation(
            subset,
            radius_bins=radius_bins,
            min_count=min_count,
            n_boot=n_boot,
            random_seed=random_seed + 1000 + offset,
        )
        binned_by_panel[key] = binned
        for threshold in comparison_thresholds:
            comparison_subset = subset.loc[pd.to_numeric(subset["FE_H"], errors="coerce") < threshold]
            comparison_binned = summarize_rotation(
                comparison_subset,
                radius_bins=radius_bins,
                min_count=min_count,
                n_boot=0,
                random_seed=random_seed + 2000 + offset,
            )
            comparison_binned_by_panel[(key, threshold)] = comparison_binned
        counts[key] = {
            "total": int(len(subset)),
            "in_situ": int((subset["population"] == "in_situ").sum()),
            "accreted": int((subset["population"] == "accreted").sum()),
        }

    fig, axes = plt.subplots(1, 3, figsize=(14.0, 4.6), sharey=True, constrained_layout=True)
    xmin = float(np.min(radius_bins))
    xmax = float(np.max(radius_bins))

    for ax, (key, title, _subset) in zip(axes, panels):
        binned = binned_by_panel[key]
        for population in ("in_situ", "accreted"):
            sub = binned.loc[binned["population"] == population].copy()
            sub = sub.loc[np.isfinite(sub["median_vphi_kms"])]
            if sub.empty:
                continue
            low = (sub["median_vphi_kms"] - sub["median_vphi_boot16_kms"]).fillna(0.0)
            high = (sub["median_vphi_boot84_kms"] - sub["median_vphi_kms"]).fillna(0.0)
            yerr = np.vstack([low, high])
            ax.errorbar(
                sub["r_mid_kpc"],
                sub["median_vphi_kms"],
                yerr=yerr,
                fmt="o-",
                lw=1.6,
                ms=4.5,
                capsize=3.0,
                color=colors[population],
                label=labels[population],
            )
            for threshold, linestyle in zip(comparison_thresholds, comparison_styles):
                comparison_binned = comparison_binned_by_panel.get((key, threshold))
                if comparison_binned is None:
                    continue
                comparison_sub = comparison_binned.loc[comparison_binned["population"] == population].copy()
                comparison_sub = comparison_sub.loc[np.isfinite(comparison_sub["median_vphi_kms"])]
                if comparison_sub.empty:
                    continue
                ax.plot(
                    comparison_sub["r_mid_kpc"],
                    comparison_sub["median_vphi_kms"],
                    ls=linestyle,
                    lw=1.7,
                    color=colors[population],
                    alpha=0.95,
                )

        panel_counts = counts[key]
        ax.set_title(
            f"{title}\n"
            f"N={panel_counts['total']:,} "
            f"({panel_counts['in_situ']:,} in-situ, {panel_counts['accreted']:,} accreted)"
        )
        ax.axhline(0.0, color="0.25", lw=1.0, ls="--", zorder=0)
        ax.set_xlim(xmin, xmax)
        ax.set_ylim(*VELOCITY_YLIM)
        ax.set_xlabel(r"$r_{\rm GC}$ [kpc]")
        ax.grid(True, color="0.9", lw=0.8)

    axes[0].set_ylabel(r"Median azimuthal velocity $V_\phi$ [km s$^{-1}$]")
    legend_handles = [
        Line2D([0], [0], color=colors["in_situ"], marker="o", lw=1.6, label="In-situ halo"),
        Line2D([0], [0], color=colors["accreted"], marker="o", lw=1.6, label="Accreted halo"),
    ]
    for threshold, linestyle in zip(comparison_thresholds, comparison_styles):
        legend_handles.append(
            Line2D([0], [0], color="0.2", lw=1.7, ls=linestyle, label=fr"$[Fe/H] < {threshold:g}$")
        )
    axes[0].legend(handles=legend_handles, frameon=False, loc="upper right")

    out_plot.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_plot, dpi=220)
    if out_pdf is not None:
        out_pdf.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(out_pdf)
    plt.close(fig)
    return counts


def plot_radius_histograms_z_slices(
    measurement: pd.DataFrame,
    radius_bins: list[float],
    comparison_feh_max: list[float],
    out_plot: Path,
    out_pdf: Path | None,
) -> None:
    colors = {"in_situ": "#b33a3a", "accreted": "#2f6fb0"}
    labels = {"in_situ": "In-situ halo", "accreted": "Accreted halo"}
    comparison_styles = ("--", ":")
    comparison_thresholds = [value for value in comparison_feh_max if np.isfinite(value)]
    panels = [
        ("All", measurement),
        (r"$|z| < 3$ kpc", measurement.loc[np.abs(measurement["galz"]) < 3.0]),
        (r"$|z| > 3$ kpc", measurement.loc[np.abs(measurement["galz"]) > 3.0]),
    ]
    bins = np.asarray(radius_bins, dtype=float)
    if len(bins) < 2:
        raise ValueError("radius_bins must contain at least two edges")

    fig, axes = plt.subplots(1, 3, figsize=(14.0, 4.3), sharey=True, constrained_layout=True)
    for ax, (title, subset) in zip(axes, panels):
        for population in ("in_situ", "accreted"):
            pop = subset.loc[subset["population"] == population]
            radius = pd.to_numeric(pop["r_gc_kpc"], errors="coerce").to_numpy(dtype=float)
            radius = radius[np.isfinite(radius)]
            if len(radius) == 0:
                continue

            counts, edges = np.histogram(radius, bins=bins)
            ax.stairs(counts, edges, color=colors[population], lw=2.0, label=labels[population])

            for threshold, linestyle in zip(comparison_thresholds, comparison_styles):
                comparison = pop.loc[pd.to_numeric(pop["FE_H"], errors="coerce") < threshold]
                comparison_radius = pd.to_numeric(comparison["r_gc_kpc"], errors="coerce").to_numpy(dtype=float)
                comparison_radius = comparison_radius[np.isfinite(comparison_radius)]
                if len(comparison_radius):
                    comparison_counts, _ = np.histogram(comparison_radius, bins=bins)
                    ax.stairs(comparison_counts, edges, color=colors[population], lw=1.7, ls=linestyle)

        ax.set_title(f"{title}\nN={len(subset):,}")
        ax.set_xlim(float(bins[0]), float(bins[-1]))
        ax.set_yscale("log")
        ax.set_ylim(bottom=0.8)
        ax.set_xlabel(r"$r_{\rm GC}$ [kpc]")
        ax.grid(True, color="0.9", lw=0.8)

    axes[0].set_ylabel("Star count")
    legend_handles = [
        Line2D([0], [0], color=colors["in_situ"], lw=2.0, label="In-situ halo"),
        Line2D([0], [0], color=colors["accreted"], lw=2.0, label="Accreted halo"),
    ]
    for threshold, linestyle in zip(comparison_thresholds, comparison_styles):
        legend_handles.append(
            Line2D([0], [0], color="0.2", lw=1.7, ls=linestyle, label=fr"$[Fe/H] < {threshold:g}$")
        )
    axes[0].legend(handles=legend_handles, frameon=False, loc="upper right")

    out_plot.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_plot, dpi=220)
    if out_pdf is not None:
        out_pdf.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(out_pdf)
    plt.close(fig)


def plot_chemical_selection(
    selected: pd.DataFrame,
    measurement_feh_max: float,
    out_plot: Path,
    out_pdf: Path | None,
) -> None:
    plot_df = selected.loc[
        np.isfinite(selected["FE_H"]) & np.isfinite(selected["AL_FE"]),
        ["FE_H", "AL_FE", "population"],
    ].copy()

    fig, ax = plt.subplots(figsize=(7.0, 5.2))
    hb = ax.hexbin(
        plot_df["FE_H"],
        plot_df["AL_FE"],
        gridsize=90,
        extent=(-2.2, 0.6, -0.65, 0.65),
        bins="log",
        mincnt=1,
        cmap="Greys",
    )
    cbar = fig.colorbar(hb, ax=ax, pad=0.015)
    cbar.set_label("log star count")

    accreted = plot_df.loc[plot_df["population"] == "accreted"]
    in_situ = plot_df.loc[plot_df["population"] == "in_situ"]
    if not accreted.empty:
        ax.scatter(
            accreted["FE_H"],
            accreted["AL_FE"],
            s=5,
            color="#2f6fb0",
            alpha=0.22,
            linewidths=0,
            label="Accreted",
        )
    if not in_situ.empty:
        in_situ_sample = in_situ.sample(n=min(len(in_situ), 12000), random_state=1)
        ax.scatter(
            in_situ_sample["FE_H"],
            in_situ_sample["AL_FE"],
            s=3,
            color="#b33a3a",
            alpha=0.08,
            linewidths=0,
            label="In-situ",
        )

    xmin, xmax = -2.2, 0.6
    ymin, ymax = -0.65, 0.65
    ax.hlines(-0.075, xmin, -0.4, color="black", lw=2.0)
    ax.vlines(-0.4, -0.075, ymax, color="black", lw=2.0)
    if np.isfinite(measurement_feh_max):
        ax.axvline(measurement_feh_max, color="#6b2e8a", lw=1.7, ls="--")

    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    ax.set_xlabel("[Fe/H]")
    ax.set_ylabel("[Al/Fe]")
    legend_handles = [
        Line2D([0], [0], marker="o", color="none", markerfacecolor="#2f6fb0", markersize=6, label="Accreted"),
        Line2D([0], [0], marker="o", color="none", markerfacecolor="#b33a3a", markersize=6, label="In-situ"),
        Line2D([0], [0], color="black", lw=2.0, label="BK22 selection boundary"),
    ]
    if np.isfinite(measurement_feh_max):
        legend_handles.append(
            Line2D([0], [0], color="#6b2e8a", lw=1.7, ls="--", label="Measurement [Fe/H] boundary")
        )
    ax.legend(handles=legend_handles, frameon=False, loc="upper left")
    ax.grid(True, color="0.9", lw=0.8)
    fig.tight_layout()

    out_plot.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_plot, dpi=220)
    if out_pdf is not None:
        out_pdf.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(out_pdf)
    plt.close(fig)


def plot_rz_and_z_distribution(
    measurement: pd.DataFrame,
    out_plot: Path,
    out_pdf: Path | None,
) -> None:
    plot_df = measurement.loc[
        np.isfinite(measurement["galr"]) & np.isfinite(measurement["galz"]),
        ["galr", "galz", "population"],
    ].copy()
    colors = {"in_situ": "#b33a3a", "accreted": "#2f6fb0"}
    labels = {"in_situ": "In-situ", "accreted": "Accreted"}

    fig, (ax_rz, ax_z) = plt.subplots(
        1,
        2,
        figsize=(9.5, 5.0),
        sharey=True,
        constrained_layout=True,
        gridspec_kw={"width_ratios": [2.2, 1.0], "wspace": 0.05},
    )

    finite_r = pd.to_numeric(plot_df["galr"], errors="coerce").to_numpy(dtype=float)
    finite_z = pd.to_numeric(plot_df["galz"], errors="coerce").to_numpy(dtype=float)
    finite_r = finite_r[np.isfinite(finite_r)]
    finite_z = finite_z[np.isfinite(finite_z)]
    r_max = float(np.nanmax(finite_r)) if len(finite_r) else 15.0
    z_abs_max = float(np.nanmax(np.abs(finite_z))) if len(finite_z) else 8.0
    z_lim = min(max(np.ceil(z_abs_max), 4.0), 16.0)
    r_lim = max(np.ceil(r_max), 8.0)

    for population in ("accreted", "in_situ"):
        sub = plot_df.loc[plot_df["population"] == population]
        ax_rz.scatter(
            sub["galr"],
            sub["galz"],
            s=8,
            alpha=0.35,
            linewidths=0,
            color=colors[population],
            label=f"{labels[population]} (N={len(sub):,})",
        )

    bins = np.linspace(-z_lim, z_lim, 45)
    for population in ("accreted", "in_situ"):
        sub = plot_df.loc[plot_df["population"] == population]
        z = pd.to_numeric(sub["galz"], errors="coerce").to_numpy(dtype=float)
        z = z[np.isfinite(z)]
        if len(z) == 0:
            continue
        density, edges = np.histogram(z, bins=bins, density=True)
        centers = 0.5 * (edges[:-1] + edges[1:])
        ax_z.step(density, centers, where="mid", color=colors[population], lw=2.0, label=labels[population])

    ax_rz.axhline(0.0, color="0.25", lw=1.0, ls="--", zorder=0)
    ax_rz.set_xlim(0.0, r_lim)
    ax_rz.set_ylim(-z_lim, z_lim)
    ax_rz.set_xlabel(r"Galactocentric cylindrical radius $R$ [kpc]")
    ax_rz.set_ylabel(r"Galactocentric height $z$ [kpc]")
    ax_rz.legend(frameon=False, loc="upper right")
    ax_rz.grid(True, color="0.9", lw=0.8)

    ax_z.axhline(0.0, color="0.25", lw=1.0, ls="--", zorder=0)
    ax_z.set_xlabel("Normalized density")
    ax_z.grid(True, color="0.9", lw=0.8)
    ax_z.tick_params(labelleft=False)

    out_plot.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_plot, dpi=220)
    if out_pdf is not None:
        out_pdf.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(out_pdf)
    plt.close(fig)


def main() -> None:
    args = parse_args()
    allstar_path = resolve_data_path(args.allstar, ALLSTAR_URL)
    astronn_path = resolve_data_path(args.astronn, ASTRONN_URL)
    vasiliev_path = resolve_data_path(args.vasiliev_gc_members, VASILIEV_GC_URL)

    allstar = read_allstar(allstar_path)
    astronn = read_astronn(astronn_path)
    vasiliev_gc_source_ids = read_vasiliev_gc_source_ids(
        vasiliev_path,
        prob_min=args.vasiliev_prob_min,
    )
    merged = allstar.merge(astronn, on="APOGEE_ID", how="inner", validate="1:1")

    selected, cut_counts = apply_bk22_clean_selection(
        merged,
        vasiliev_gc_source_ids=vasiliev_gc_source_ids,
        abundance_err_max=args.abundance_err_max,
        velocity_err_max=args.velocity_err_max,
        distance_max_kpc=args.distance_max_kpc,
    )

    measurement = selected.copy()
    measurement_feh_max = args.measurement_feh_max
    if np.isfinite(measurement_feh_max):
        measurement = measurement.loc[measurement["FE_H"] < measurement_feh_max].copy()

    radius_bins = args.radius_bins if args.radius_bins is not None else automatic_radius_bins(measurement)

    binned = summarize_rotation(
        measurement,
        radius_bins=radius_bins,
        min_count=args.min_count,
        n_boot=args.bootstrap,
        random_seed=args.random_seed,
    )

    args.out_csv.parent.mkdir(parents=True, exist_ok=True)
    binned.to_csv(args.out_csv, index=False)
    plot_rotation(binned, args.out_plot, args.out_pdf)
    plot_chemical_selection(selected, measurement_feh_max, args.out_chem_plot, args.out_chem_pdf)
    plot_rz_and_z_distribution(measurement, args.out_rz_plot, args.out_rz_pdf)
    zsplit_counts = plot_rotation_z_slices(
        measurement,
        radius_bins=radius_bins,
        min_count=args.min_count,
        n_boot=args.bootstrap,
        random_seed=args.random_seed,
        comparison_feh_max=args.zsplit_comparison_feh_max,
        out_plot=args.out_zsplit_plot,
        out_pdf=args.out_zsplit_pdf,
    )
    plot_radius_histograms_z_slices(
        measurement,
        radius_bins=radius_bins,
        comparison_feh_max=args.zsplit_comparison_feh_max,
        out_plot=args.out_rhist_plot,
        out_pdf=args.out_rhist_pdf,
    )

    summary = {
        "allstar_path": str(allstar_path),
        "astronn_path": str(astronn_path),
        "vasiliev_gc_members_path": str(vasiliev_path),
        "source_urls": {
            "allstar": ALLSTAR_URL,
            "astronn": ASTRONN_URL,
            "vasiliev_gc_members": VASILIEV_GC_URL,
        },
        "vasiliev_prob_min": args.vasiliev_prob_min,
        "vasiliev_gc_source_ids": int(len(vasiliev_gc_source_ids)),
        "merged_rows": int(len(merged)),
        "clean_rows": int(len(selected)),
        "measurement_rows": int(len(measurement)),
        "measurement_feh_max": None if not np.isfinite(measurement_feh_max) else measurement_feh_max,
        "population_counts_clean": selected["population"].value_counts().to_dict(),
        "population_counts_measurement": measurement["population"].value_counts().to_dict(),
        "population_counts_z_slices": zsplit_counts,
        "zsplit_comparison_feh_max": [
            float(value) for value in args.zsplit_comparison_feh_max if np.isfinite(value)
        ],
        "cut_counts": cut_counts,
        "radius_bins_kpc": radius_bins,
        "min_count_per_bin": args.min_count,
        "bootstrap_resamples": args.bootstrap,
        "output_csv": str(args.out_csv),
        "output_plot": str(args.out_plot),
        "output_pdf": None if args.out_pdf is None else str(args.out_pdf),
        "output_chemical_selection_plot": str(args.out_chem_plot),
        "output_chemical_selection_pdf": None if args.out_chem_pdf is None else str(args.out_chem_pdf),
        "output_rz_z_distribution_plot": str(args.out_rz_plot),
        "output_rz_z_distribution_pdf": None if args.out_rz_pdf is None else str(args.out_rz_pdf),
        "output_zsplit_rotation_plot": str(args.out_zsplit_plot),
        "output_zsplit_rotation_pdf": None if args.out_zsplit_pdf is None else str(args.out_zsplit_pdf),
        "output_radius_histograms_plot": str(args.out_rhist_plot),
        "output_radius_histograms_pdf": None if args.out_rhist_pdf is None else str(args.out_rhist_pdf),
        "notes": [
            "APOGEE chemistry and quality flags are from DR17 allStarLite.",
            "Coordinates, velocities, and velocity errors are from AstroNN DR17.",
            "Known GC candidates are removed by Gaia EDR3 source-id matches to the Zenodo Vasiliev & Baumgardt 2021 cluster-membership catalogue.",
            "No orbital-energy window is applied.",
            "Population split follows BK22 text: [Fe/H] > -0.4 is in-situ; at lower [Fe/H], [Al/Fe] >= -0.075 is in-situ and [Al/Fe] < -0.075 is accreted.",
            "The plotted radius is spherical r_gc = sqrt(galr^2 + galz^2); galr itself is the AstroNN cylindrical radius.",
            "AstroNN dist is converted from pc to kpc before applying the heliocentric distance cut.",
        ],
    }
    args.out_json.parent.mkdir(parents=True, exist_ok=True)
    args.out_json.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")

    print(f"Merged rows: {len(merged):,}")
    print(f"allStarLite path: {allstar_path}")
    print(f"AstroNN path: {astronn_path}")
    print(f"Vasiliev GC path: {vasiliev_path}")
    print(f"Vasiliev GC source ids with prob > {args.vasiliev_prob_min}: {len(vasiliev_gc_source_ids):,}")
    print(f"Clean rows: {len(selected):,}")
    print(f"Measurement rows: {len(measurement):,}")
    print("Measurement population counts:")
    print(measurement["population"].value_counts().to_string())
    print(f"Saved {args.out_csv}")
    print(f"Saved {args.out_json}")
    print(f"Saved {args.out_plot}")
    print(f"Saved {args.out_chem_plot}")
    print(f"Saved {args.out_rz_plot}")
    print(f"Saved {args.out_zsplit_plot}")
    print(f"Saved {args.out_rhist_plot}")


if __name__ == "__main__":
    main()
