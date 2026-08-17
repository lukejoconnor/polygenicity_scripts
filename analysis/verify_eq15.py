#!/usr/bin/env python3
"""Reproduce Figure 3/Table S2 estimates and compare them with literal Eq. 15.

This is an independent numerical check of MATLAB/estimate_polygenicity.m. It
does not replace the MATLAB analysis. It is useful on systems without MATLAB.
"""

import argparse
import csv
from pathlib import Path

import numpy as np
from scipy.io import loadmat
from scipy.special import logsumexp


LIPID_ROWS = np.array([0, 1, 2, 4])
LIPID_NEFF = {
    "Triglycerides": 431952,
    "LDL": 434347,
    "HbA1c": 411840,
    "Hair color": 452720,
}


def load_inputs(matfiles):
    opts = dict(squeeze_me=True, struct_as_record=False)
    fmr1 = loadmat(matfiles / "FMRestimates_32traits.mat", **opts)
    fmr2 = loadmat(matfiles / "FMR_lipids.mat", **opts)
    names = loadmat(matfiles / "name_abbreviations.mat", **opts)
    sigma2_z = np.vstack([fmr1["ss_est"], fmr2["ss_est"][LIPID_ROWS]])
    omega_jk = np.concatenate(
        [fmr1["ww_jk"], fmr2["ww_jk"][LIPID_ROWS]], axis=0
    )
    ld4m = list(fmr1["LD4Mout"]) + [fmr2["LD4Mout"][i] for i in LIPID_ROWS]
    lipid_traits = [str(fmr2["traits"][i]) for i in LIPID_ROWS]
    neff = np.r_[fmr1["Neff"], [LIPID_NEFF[name] for name in lipid_traits]]
    traits = [str(x) for x in names["traits"]]
    return sigma2_z, omega_jk, ld4m, neff, traits, float(fmr1["mm"])


def calculate(matfiles):
    sigma2_z, omega_jk, ld4m, neff, traits, marker_count = load_inputs(matfiles)
    n_traits, _, n_jk = omega_jk.shape
    legacy = np.empty((n_traits, n_jk, 3))
    eq15 = np.empty_like(legacy)
    normalized_softmax_stable = np.empty((n_traits, n_jk))

    for trait in range(n_traits):
        omega = omega_jk[trait].T
        legacy_weight_sum = omega.sum(axis=1)
        omega = omega / legacy_weight_sum[:, None]
        nh2 = marker_count * np.mean(np.asarray(ld4m[trait].cov, dtype=float))
        h2 = nh2 / neff[trait]

        # Published implementation: normalized component variances and a
        # transformed generator. Preserve its 1e-256 softmax floor exactly.
        normalized_sigma2 = sigma2_z[trait] / nh2
        legacy[trait, :, 0] = legacy_weight_sum / np.exp(
            np.sum(omega * np.log(normalized_sigma2), axis=1)
        )
        legacy[trait, :, 1] = legacy_weight_sum / np.sum(
            omega * normalized_sigma2, axis=1
        )
        legacy[trait, :, 2] = -legacy_weight_sum * np.log(
            np.sum(
                omega
                * np.maximum(1e-256, np.exp(-1 / normalized_sigma2)),
                axis=1,
            )
        )
        normalized_softmax_stable[trait] = -legacy_weight_sum * logsumexp(
            -1 / normalized_sigma2, axis=1, b=omega
        )

        # Literal Eq. 15. sigma2 is in phenotypic-variance units and the
        # absolute component weights are h2*omega.
        sigma2 = sigma2_z[trait] / neff[trait]
        eq15[trait, :, 0] = h2 * np.exp(
            np.sum(omega * np.log(1 / sigma2), axis=1)
        )
        eq15[trait, :, 1] = h2 / np.sum(omega * sigma2, axis=1)
        eq15[trait, :, 2] = -h2 * logsumexp(
            -1 / sigma2, axis=1, b=omega
        )

    return traits, legacy, normalized_softmax_stable, eq15


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--matfiles", type=Path, required=True)
    parser.add_argument("--output", type=Path, default=Path("eq15_comparison.csv"))
    parser.add_argument(
        "--published-table",
        type=Path,
        help="Optional Table S2 xlsx used to verify the legacy reproduction",
    )
    args = parser.parse_args()

    traits, legacy, normalized_softmax_stable, eq15 = calculate(args.matfiles)
    old_log = np.mean(np.log10(legacy), axis=1)
    new_log = np.mean(np.log10(eq15), axis=1)
    old_raw = 10**old_log
    new_raw = 10**new_log
    old_se = np.std(np.log10(legacy), axis=1, ddof=1) * np.sqrt(102)
    new_se = np.std(np.log10(eq15), axis=1, ddof=1) * np.sqrt(102)
    stable_old_log = np.mean(np.log10(normalized_softmax_stable), axis=1)
    stable_old_raw = 10**stable_old_log

    if args.published_table:
        import openpyxl

        sheet = openpyxl.load_workbook(
            args.published_table, data_only=True, read_only=True
        )["S2"]
        published = np.array(
            [
                [row[1], row[2], row[3], row[4], row[5], row[6]]
                for row in sheet.iter_rows(min_row=4, max_row=39, values_only=True)
            ],
            dtype=float,
        )
        print(
            "published reproduction: max |point change| = "
            f"{np.max(np.abs(old_log - published[:, :3])):.6g}; "
            "max |SE change| = "
            f"{np.max(np.abs(old_se - published[:, 3:])):.6g}"
        )
        np.testing.assert_allclose(old_log, published[:, :3], atol=1e-12, rtol=0)
        np.testing.assert_allclose(old_se, published[:, 3:], atol=1e-12, rtol=0)

    np.testing.assert_allclose(
        eq15[:, :, :2], legacy[:, :, :2], atol=1e-9, rtol=1e-12
    )
    if not np.all(eq15[:, :, 0] >= eq15[:, :, 1]):
        raise AssertionError("Entropy polygenicity must be >= effective polygenicity")
    if not np.all(eq15[:, :, 1] >= eq15[:, :, 2]):
        raise AssertionError("Effective polygenicity must be >= softmax polygenicity")

    with args.output.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            [
                "trait",
                "legacy_softmax",
                "legacy_softmax_without_floor",
                "eq15_softmax",
                "raw_change",
                "relative_change",
                "legacy_log10_softmax",
                "eq15_log10_softmax",
                "log10_change",
                "legacy_log10_se",
                "eq15_log10_se",
                "log10_se_change",
            ]
        )
        for i, trait in enumerate(traits):
            writer.writerow(
                [
                    trait,
                    old_raw[i, 2],
                    stable_old_raw[i],
                    new_raw[i, 2],
                    new_raw[i, 2] - old_raw[i, 2],
                    new_raw[i, 2] / old_raw[i, 2] - 1,
                    old_log[i, 2],
                    new_log[i, 2],
                    new_log[i, 2] - old_log[i, 2],
                    old_se[i, 2],
                    new_se[i, 2],
                    new_se[i, 2] - old_se[i, 2],
                ]
            )

    for column, name in [(0, "entropy"), (1, "effective")]:
        print(
            f"{name}: max |jackknife raw change| = "
            f"{np.max(np.abs(eq15[:, :, column] - legacy[:, :, column])):.6g}; "
            f"max |log10 point change| = "
            f"{np.max(np.abs(new_log[:, column] - old_log[:, column])):.6g}"
        )
    dlog = new_log[:, 2] - old_log[:, 2]
    relative = new_raw[:, 2] / old_raw[:, 2] - 1
    print(
        "softmax log10 change (min/median/max): "
        f"{dlog.min():.6g} / {np.median(dlog):.6g} / {dlog.max():.6g}"
    )
    print(
        "softmax relative change (min/median/max): "
        f"{relative.min():.6g} / {np.median(relative):.6g} / {relative.max():.6g}"
    )
    scaling_only = new_log[:, 2] - stable_old_log
    print(
        "softmax scaling-only log10 change after removing legacy floor "
        "(min/median/max): "
        f"{scaling_only.min():.6g} / {np.median(scaling_only):.6g} / "
        f"{scaling_only.max():.6g}"
    )
    print(f"wrote {args.output}")


if __name__ == "__main__":
    main()
