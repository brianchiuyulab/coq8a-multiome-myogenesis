"""Define accessibility timing from GSE109828 only, before COQ8A testing.

Rates are standardized across five pooled accessible-site-depth strata within
each experiment. A three-degree-of-freedom Wald test evaluates the 24/48/72 h
profile against 0 h, with BH correction across eligible overlapping external
regions. Opening regions require a positive 72 h change in both experiments.
Half-rise time is linearly interpolated from actual sampling times, not called
pseudotime. Endpoint replication does not establish replicate onset timing.
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse, stats
from statsmodels.stats.multitest import multipletests


def profile(root, exp, promoter_width):
    cells = pd.read_csv(root / f"HSMM{exp}_cells.tsv.gz", sep="\t")
    regions = pd.read_csv(root / f"HSMM{exp}_selected_regions.tsv", sep="\t")
    matrix = sparse.load_npz(root / f"HSMM{exp}_selected_binary.npz")
    mask = cells[f"pass_{promoter_width}kb"].to_numpy()
    cells = cells[mask].reset_index(drop=True)
    matrix = matrix[mask].astype(float)
    cells["depth_stratum"] = pd.qcut(
        cells.open_sites, 5, labels=False, duplicates="drop"
    )
    hours = sorted(cells.hour.unique())
    table = pd.crosstab(cells.hour, cells.depth_stratum)
    strata = table.columns[(table > 0).all(axis=0)].tolist()
    weights = cells.depth_stratum.value_counts().reindex(strata).to_numpy(dtype=float)
    weights /= weights.sum()
    rates = []
    variances = []
    audit = []
    for hour in hours:
        rate = np.zeros(matrix.shape[1])
        variance = rate.copy()
        for stratum, weight in zip(strata, weights):
            chosen = (
                (cells.hour == hour) & (cells.depth_stratum == stratum)
            ).to_numpy()
            n = int(chosen.sum())
            hits = np.asarray(matrix[chosen].sum(axis=0)).ravel()
            raw = hits / n
            # Jeffreys smoothing stabilizes variance at all-zero/all-one strata;
            # plotted/compared accessibility proportions remain unsmoothed.
            smooth = (hits + 0.5) / (n + 1)
            rate += weight * raw
            variance += weight**2 * smooth * (1 - smooth) / (n + 1)
            audit.append(
                dict(
                    experiment=exp,
                    promoter_kb=promoter_width,
                    hour=hour,
                    stratum=stratum,
                    n_cells=n,
                    standardization_weight=weight,
                    min_open_sites=int(cells.loc[chosen, "open_sites"].min()),
                    max_open_sites=int(cells.loc[chosen, "open_sites"].max()),
                )
            )
        rates.append(rate)
        variances.append(variance)
    result = regions.copy()
    for hour, rate, var in zip(hours, rates, variances):
        result[f"p{hour}"] = rate
        result[f"v{hour}"] = var
    result["detected_cells"] = np.asarray(matrix.sum(axis=0)).ravel().astype(int)
    result["eligible"] = result.detected_cells >= max(10, 0.01 * len(cells))
    return result, pd.DataFrame(audit)


def half_rise(values):
    delta = np.asarray(values) - values[0]
    height = delta[1:].max()
    if height <= 0:
        return np.nan
    threshold = 0.5 * height
    for i in range(1, 4):
        if delta[i] >= threshold:
            return 24 * (i - 1) + 24 * (threshold - delta[i - 1]) / (
                delta[i] - delta[i - 1]
            )
    return np.nan


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for arg in ["prepared", "tables", "out"]:
        p.add_argument("--" + arg, type=Path, required=True)
    a = p.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)
    candidates = pd.read_csv(a.tables / "candidate_peak_gene.tsv.gz", sep="\t")
    overlaps = pd.read_csv(a.prepared / "HSMM1_overlaps.tsv.gz", sep="\t")
    # A target region is assigned the best geometrical external overlap, with
    # no use of its expression, association p value or accessibility effect.
    best = overlaps.sort_values(
        ["reciprocal_overlap", "overlap_bp", "external_peak"],
        ascending=[False, False, True],
    ).drop_duplicates("peak")
    all_sets = []
    definitions = []
    audits = []
    mapped_rows = []
    for promoter in [2, 1]:
        first, a1 = profile(a.prepared, 1, promoter)
        second, a2 = profile(a.prepared, 2, promoter)
        audits.extend([a1, a2])
        ext = first.merge(
            second[["external_peak", "p0", "p72", "eligible"]].rename(
                columns={"p0": "rep_p0", "p72": "rep_p72", "eligible": "rep_eligible"}
            ),
            on="external_peak",
            how="left",
        )
        d = ext[["p24", "p48", "p72"]].to_numpy() - ext.p0.to_numpy()[:, None]
        v = ext[["v24", "v48", "v72"]].to_numpy()
        v0 = ext.v0.to_numpy()
        wald = (d * d / v).sum(1) - (d / v).sum(1) ** 2 / (1 / v0 + (1 / v).sum(1))
        ext["p_dynamic"] = stats.chi2.sf(np.maximum(wald, 0), 3)
        ext["q_dynamic"] = np.nan
        eligible = ext.eligible.to_numpy()
        ext.loc[eligible, "q_dynamic"] = multipletests(
            ext.loc[eligible, "p_dynamic"], method="fdr_bh"
        )[1]
        ext["delta72"] = ext.p72 - ext.p0
        ext["rep_delta72"] = ext.rep_p72 - ext.rep_p0
        ext["opening"] = (
            (ext.q_dynamic < 0.05)
            & (ext.delta72 > 0)
            & (ext.rep_delta72 > 0)
            & ext.rep_eligible.fillna(False)
        )
        ext["closing"] = (
            (ext.q_dynamic < 0.05)
            & (ext.delta72 < 0)
            & (ext.rep_delta72 < 0)
            & ext.rep_eligible.fillna(False)
        )
        ext["half_rise_hour"] = [
            half_rise(x) for x in ext[["p0", "p24", "p48", "p72"]].to_numpy()
        ]
        # Low-change is an observed-profile comparator, not a claim of equivalence.
        ratio = ext[["p0", "p24", "p48", "p72"]].max(axis=1) / ext[
            ["p0", "p24", "p48", "p72"]
        ].min(axis=1).replace(0, np.nan)
        ext["low_change_profile"] = ext.eligible & (ratio <= 1.2)
        ext["promoter_kb"] = promoter
        ext.to_csv(
            a.out / f"external_profiles_promoter{promoter}kb.tsv.gz",
            sep="\t",
            index=False,
        )
        mapped = best.merge(ext, on="external_peak", how="left", validate="many_to_one")
        mapped_rows.append(mapped)
        for overlap in [0.5, 0.25]:
            eligiblemap = mapped[
                (mapped.reciprocal_overlap >= overlap) & mapped.eligible
            ].copy()
            opening = eligiblemap[eligiblemap.opening].copy()
            conditions = {
                "All_external_mapped": eligiblemap,
                "All_opening": opening,
                "Early_by24h": opening[opening.half_rise_hour <= 24],
                "Middle_24to48h": opening[
                    (opening.half_rise_hour > 24) & (opening.half_rise_hour <= 48)
                ],
                "Late_after48h": opening[opening.half_rise_hour > 48],
                "Closing": eligiblemap[eligiblemap.closing],
                "Low_change_profile": eligiblemap[eligiblemap.low_change_profile],
            }
            # Quantiles use unique external regions, not repeated target mappings.
            unique_open = opening.drop_duplicates("external_peak")
            for fraction in [0.2, 0.25, 1 / 3]:
                threshold = unique_open.half_rise_hour.quantile(fraction)
                conditions[f"Earliest_{round(100*fraction)}pct"] = opening[
                    opening.half_rise_hour <= threshold
                ]
            for name, frame in conditions.items():
                set_id = f"P{promoter}_O{int(overlap*100)}_{name}"
                peaks = sorted(frame.peak.unique())
                genes = sorted(
                    candidates.loc[candidates.peak.isin(peaks), "gene"].unique()
                )
                definitions.append(
                    dict(
                        set_id=set_id,
                        module=name,
                        promoter_kb=promoter,
                        reciprocal_overlap=overlap,
                        n_peaks=len(peaks),
                        n_genes=len(genes),
                        genes=";".join(genes),
                    )
                )
                all_sets.extend(dict(set_id=set_id, peak=peak) for peak in peaks)
    pd.concat(audits, ignore_index=True).to_csv(
        a.out / "external_cell_qc_depth_strata.tsv", sep="\t", index=False
    )
    pd.concat(mapped_rows, ignore_index=True).to_csv(
        a.out / "target_temporal_annotations.tsv.gz", sep="\t", index=False
    )
    pd.DataFrame(definitions).to_csv(
        a.out / "temporal_set_definitions.tsv", sep="\t", index=False
    )
    pd.DataFrame(all_sets).to_csv(
        a.out / "temporal_set_memberships.tsv.gz", sep="\t", index=False
    )
    print(
        pd.DataFrame(definitions)
        .query("promoter_kb==2 and reciprocal_overlap==0.5")[
            ["module", "n_peaks", "n_genes"]
        ]
        .to_string(index=False)
    )


if __name__ == "__main__":
    main()
