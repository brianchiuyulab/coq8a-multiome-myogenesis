"""Test fixed external temporal peak sets in matched GSE208248 nuclei.

Both COQ8A contrasts and both TSS gates use the same frozen memberships.
P values refer to matched nuclei, not four independent human donors. Programme
BH adjustment is within one QC/contrast/stage setting; sensitivity settings are
not pooled into a multiple-testing family. Gene-level screening has its own BH
family. Raw p values, directions and complete outputs are retained.
"""

import argparse
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from scipy import sparse, stats
from statsmodels.stats.multitest import multipletests

from multiome_core import h5_for, read_barcodes, paired_stats


def summarize(frame):
    high = frame.high_open.to_numpy()
    low = frame.low_open.to_numpy()
    s = paired_stats(high - low)
    return dict(
        n_pairs=len(frame),
        high_open_pct=100 * high.mean(),
        low_open_pct=100 * low.mean(),
        fold_open=high.mean() / low.mean() if low.mean() > 0 else np.nan,
        delta_pp=100 * s["difference"],
        delta_ci_low_pp=100 * s["ci_low"],
        delta_ci_high_pp=100 * s["ci_high"],
        p_pair=s["p_pair"],
    )


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for arg in ["h5-root", "tables", "temporal", "out"]:
        p.add_argument("--" + arg, type=Path, required=True)
    a = p.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)
    definitions = pd.read_csv(a.temporal / "temporal_set_definitions.tsv", sep="\t")
    definitions = definitions[definitions.n_peaks > 0].copy()
    definitions["scope"] = "programme"
    members = pd.read_csv(a.temporal / "temporal_set_memberships.tsv.gz", sep="\t")
    candidates = pd.read_csv(a.tables / "candidate_peak_gene.tsv.gz", sep="\t")
    additions = []
    newmembers = []
    for module in [
        "All_opening",
        "Early_by24h",
        "Middle_24to48h",
        "Earliest_25pct",
        "Late_after48h",
    ]:
        peaks = set(members.loc[members.set_id == "P2_O50_" + module, "peak"])
        for gene, frame in candidates[candidates.peak.isin(peaks)].groupby("gene"):
            selected = sorted(frame.peak.unique())
            set_id = f"Gene_{module}_{gene}"
            additions.append(
                dict(
                    set_id=set_id,
                    module=module,
                    promoter_kb=2,
                    reciprocal_overlap=0.5,
                    n_peaks=len(selected),
                    n_genes=1,
                    genes=gene,
                    scope="gene",
                )
            )
            newmembers.extend(dict(set_id=set_id, peak=peak) for peak in selected)
    definitions = pd.concat([definitions, pd.DataFrame(additions)], ignore_index=True)
    members = pd.concat([members, pd.DataFrame(newmembers)], ignore_index=True)
    definitions["column"] = np.arange(len(definitions))
    definitions.to_csv(a.out / "tested_set_definitions.tsv", sep="\t", index=False)
    members.to_csv(a.out / "tested_set_memberships.tsv.gz", sep="\t", index=False)
    colindex = dict(zip(definitions.set_id, definitions.column))
    mapping = pd.read_csv(a.tables / "consensus_peak_map.tsv.gz", sep="\t").set_index(
        "peak"
    )
    pairs = pd.read_csv(a.tables / "matched_pairs.tsv.gz", sep="\t")
    allscores = []
    for gsm, ps in pairs.groupby("gsm"):
        barcodes = sorted(set(ps.high_barcode) | set(ps.low_barcode))
        with h5py.File(h5_for(a.h5_root, gsm)) as h5:
            _, matrix, _, peaks = read_barcodes(h5, barcodes)
        bidx = {b: i for i, b in enumerate(barcodes)}
        pidx = {str(x): i for i, x in enumerate(peaks)}
        local = members.copy()
        target = (
            local.peak
            if gsm == "GSM6339597"
            else local.peak.map(mapping[gsm + "_peak"])
        )
        local["row"] = target.map(pidx)
        local["column"] = local.set_id.map(colindex)
        assert local.row.notna().all()
        local = local.drop_duplicates(["row", "column"])
        denominator = local.groupby("column").size()
        weights = sparse.csr_matrix(
            (
                1 / local.column.map(denominator).to_numpy(),
                (local.row.astype(int), local.column),
            ),
            shape=(len(peaks), len(definitions)),
        )
        values = (matrix @ weights).toarray()
        for r in ps.itertuples():
            score = pd.DataFrame(
                {
                    "set_id": definitions.set_id,
                    "high_open": values[bidx[r.high_barcode]],
                    "low_open": values[bidx[r.low_barcode]],
                }
            )
            for key in ["gsm", "source", "stage", "gate", "contrast", "pair"]:
                score[key] = getattr(r, key)
            allscores.append(score)
        print(gsm, "sets", len(definitions), "matched comparisons", len(ps), flush=True)
    scores = pd.concat(allscores, ignore_index=True)
    scores.to_csv(a.out / "matched_pair_temporal_scores.tsv.gz", sep="\t", index=False)
    libraries = []
    for key, frame in scores.groupby(
        ["set_id", "gate", "contrast", "gsm", "source", "stage"]
    ):
        row = dict(zip(["set_id", "gate", "contrast", "gsm", "source", "stage"], key))
        row.update(summarize(frame))
        libraries.append(row)
    libraries = pd.DataFrame(libraries)
    libraries.loc[libraries.delta_pp.abs() < 1e-10, "delta_pp"] = 0
    libraries.merge(definitions, on="set_id").to_csv(
        a.out / "temporal_effects_by_library.tsv.gz", sep="\t", index=False
    )
    outputs = []
    contrasts = []
    for stage in ["pooled", "stem", "differentiated"]:
        data = scores if stage == "pooled" else scores[scores.stage == stage]
        if data.empty:
            continue
        for key, frame in data.groupby(["set_id", "gate", "contrast"]):
            row = dict(zip(["set_id", "gate", "contrast"], key))
            row["stage_analysis"] = stage
            row.update(summarize(frame))
            ll = libraries[
                (libraries.set_id == key[0])
                & (libraries.gate == key[1])
                & (libraries.contrast == key[2])
            ]
            if stage != "pooled":
                ll = ll[ll.stage == stage]
            row["positive_libraries"] = int((ll.delta_pp > 1e-10).sum())
            row["n_libraries"] = len(ll)
            row["n_source_lines"] = frame.source.nunique()
            outputs.append(row)
        # Compare relative effects in early and late regions on the same pairs.
        for gate, contrast in (
            data[["gate", "contrast"]]
            .drop_duplicates()
            .itertuples(index=False, name=None)
        ):
            for promoter in [2, 1]:
                for overlap in [50, 25]:
                    prefix = f"P{promoter}_O{overlap}_"
                    subset = data[(data.gate == gate) & (data.contrast == contrast)]
                    early = subset[subset.set_id == prefix + "Early_by24h"]
                    late = subset[subset.set_id == prefix + "Late_after48h"]
                    joint = early.merge(
                        late, on=["gsm", "pair", "source"], suffixes=("_early", "_late")
                    )
                    if len(joint) < 2:
                        continue
                    le = joint.low_open_early.mean()
                    ll = joint.low_open_late.mean()
                    if min(le, ll) <= 0:
                        continue
                    difference = (joint.high_open_early - joint.low_open_early) / le - (
                        joint.high_open_late - joint.low_open_late
                    ) / ll
                    result = paired_stats(difference)
                    rng = np.random.default_rng(208248)
                    array = joint[
                        [
                            "high_open_early",
                            "low_open_early",
                            "high_open_late",
                            "low_open_late",
                        ]
                    ].to_numpy()
                    indices = [
                        np.flatnonzero(joint.gsm.to_numpy() == gsm)
                        for gsm in joint.gsm.unique()
                    ]
                    boot = []
                    for _ in range(2000):
                        selected = np.concatenate(
                            [rng.choice(ix, len(ix), replace=True) for ix in indices]
                        )
                        h1, l1, h2, l2 = array[selected].mean(axis=0)
                        if l1 > 0 and l2 > 0:
                            boot.append(h1 / l1 - h2 / l2)
                    contrasts.append(
                        dict(
                            stage_analysis=stage,
                            gate=gate,
                            contrast=contrast,
                            promoter_kb=promoter,
                            reciprocal_overlap=overlap / 100,
                            n_pairs=len(joint),
                            early_minus_late_fold=joint.high_open_early.mean() / le
                            - joint.high_open_late.mean() / ll,
                            bootstrap_ci_low=np.quantile(boot, 0.025),
                            bootstrap_ci_high=np.quantile(boot, 0.975),
                            p_pair_scaled_difference=result["p_pair"],
                        )
                    )
    result = pd.DataFrame(outputs).merge(definitions, on="set_id")
    result["q_within_family"] = np.nan
    grouping = [
        "gate",
        "contrast",
        "stage_analysis",
        "promoter_kb",
        "reciprocal_overlap",
        "scope",
    ]
    # Gene families are distinct temporal-module screens; programmes are a
    # single family per sensitivity setting, without correction across settings.
    result["family_module"] = np.where(
        result.scope == "gene", result.module, "all_programmes"
    )
    for _, ids in result.groupby(grouping + ["family_module"]).groups.items():
        valid = result.loc[ids, "p_pair"].dropna().index
        if len(valid):
            result.loc[valid, "q_within_family"] = multipletests(
                result.loc[valid, "p_pair"], method="fdr_bh"
            )[1]
    result.to_csv(a.out / "temporal_effect_grid.tsv", sep="\t", index=False)
    pd.DataFrame(contrasts).to_csv(
        a.out / "early_vs_late_comparison.tsv", sep="\t", index=False
    )
    print(
        result.query(
            "scope=='programme' and promoter_kb==2 and reciprocal_overlap==0.5 and stage_analysis=='pooled'"
        )[
            [
                "module",
                "gate",
                "contrast",
                "n_peaks",
                "fold_open",
                "p_pair",
                "positive_libraries",
            ]
        ].to_string(
            index=False
        )
    )


if __name__ == "__main__":
    main()
