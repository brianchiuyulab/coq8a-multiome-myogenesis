"""Define enhancer sets from human chromatin states and independent MYOD ChIP."""

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse, stats
from statsmodels.stats.multitest import multipletests


def annotate(mapping, reference):
    output = []
    grouped = {c: g.sort_values("start") for c, g in reference.groupby("chrom")}
    for r in mapping.itertuples():
        g = grouped.get(r.chrom19)
        if g is None:
            output.append((r.peak, 0, 0.0))
            continue
        s = g[(g.start < r.end19) & (g.end > r.start19)]
        covered = []
        for row in s.itertuples():
            covered.append((max(r.start19, row.start), min(r.end19, row.end)))
        union = []
        for start, end in sorted(covered):
            if union and start <= union[-1][1]:
                union[-1][1] = max(union[-1][1], end)
            else:
                union.append([start, end])
        output.append(
            (r.peak, len(s), sum(e - s for s, e in union) / (r.end19 - r.start19))
        )
    return pd.DataFrame(output, columns=["peak", "overlap_regions", "covered_fraction"])


def main(a):
    a.out.mkdir(parents=True, exist_ok=True)
    mapping = pd.read_csv(
        a.root / "results/temporal/candidate_hg19_mapping.tsv", sep="\t"
    )
    assert len(mapping) == 5097
    mapping = mapping[mapping.mapping_pass].copy()
    print("Valid unique hg19 mappings:", len(mapping), "of 5097", flush=True)
    states = pd.read_csv(
        a.hmm,
        sep="\t",
        header=None,
        names=[
            "chrom",
            "start",
            "end",
            "state",
            "score",
            "strand",
            "thick_start",
            "thick_end",
            "rgb",
        ],
    )
    enh = states[states.state.isin(["4_Strong_Enhancer", "5_Strong_Enhancer"])]
    h = annotate(mapping, enh).rename(
        columns={
            "overlap_regions": "strong_enhancer_segments",
            "covered_fraction": "strong_enhancer_fraction",
        }
    )
    chip = []
    manifest = []
    for name in [
        "GSM1218849_MB135GMMD.peak.txt.gz",
        "GSM1218850_MB135DMMD.peak.txt.gz",
    ]:
        file = a.chip / name
        d = pd.read_csv(file, sep="\t", header=None).iloc[:, :3]
        d.columns = ["chrom", "start", "end"]
        chip.append(d)
        manifest.append(
            dict(file=name, sha256=hashlib.sha256(file.read_bytes()).hexdigest())
        )
    b = annotate(mapping, pd.concat(chip)).rename(
        columns={
            "overlap_regions": "MYOD_chip_calls",
            "covered_fraction": "MYOD_bound_fraction",
        }
    )
    ann = mapping.merge(h, on="peak").merge(b, on="peak")
    ann["strong_enhancer"] = ann.strong_enhancer_segments > 0
    ann["MYOD_bound_enhancer"] = ann.strong_enhancer & (ann.MYOD_chip_calls > 0)
    ann.to_csv(a.out / "external_peak_annotations.tsv", sep="\t", index=False)
    original = pd.read_csv(
        a.root / "results/tables/candidate_peak_gene.tsv.gz", sep="\t"
    )
    sets = {
        "HSMM_strong_enhancers": set(ann.loc[ann.strong_enhancer, "peak"]),
        "HSMM_strong_MYOD_bound_enhancers": set(
            ann.loc[ann.MYOD_bound_enhancer, "peak"]
        ),
    }
    for gene in ["MYOD1", "MYF5", "MYF6", "MYOG"]:
        sets[gene + "_strong_MYOD_bound"] = sets[
            "HSMM_strong_MYOD_bound_enhancers"
        ] & set(original.loc[original.gene.eq(gene), "peak"])
    if sets["MYF5_strong_MYOD_bound"] == sets["MYF6_strong_MYOD_bound"]:
        sets["MYF5_MYF6_shared_strong_MYOD_bound"] = sets.pop("MYF5_strong_MYOD_bound")
        sets.pop("MYF6_strong_MYOD_bound")
    core_sets = [
        name for name in sets if name.startswith(("MYOD1_", "MYF5_", "MYF6_", "MYOG_"))
    ]
    sets["MRF_core_union"] = set().union(*(sets[name] for name in core_sets))
    pd.DataFrame(
        [(name, peak) for name, peaks in sets.items() for peak in sorted(peaks)],
        columns=["set_id", "peak"],
    ).to_csv(a.out / "set_membership.tsv", sep="\t", index=False)
    pd.DataFrame(
        [(name, len(peaks)) for name, peaks in sets.items()],
        columns=["set_id", "n_peaks"],
    ).to_csv(a.out / "set_sizes.tsv", sep="\t", index=False)
    effects = pd.read_csv(a.root / "results/unstratified/ATAC_all_5097.tsv", sep="\t")
    pieces = []
    for name, peaks in sets.items():
        if not peaks:
            continue
        sub = effects[effects.peak.isin(peaks)].copy()
        sub["set_id"] = name
        sub["q_within_external_set"] = multipletests(
            sub.p_pair_binomial, method="fdr_bh"
        )[1]
        pieces.append(sub)
    pd.concat(pieces).to_csv(
        a.out / "peak_effects_by_external_set.tsv", sep="\t", index=False, na_rep="NA"
    )

    peaks = (a.cache / "peaks.txt").read_text().splitlines()
    idx = {p: i for i, p in enumerate(peaks)}
    ps = pd.read_csv(a.root / "results/tables/matched_pairs.tsv.gz", sep="\t")
    ps = ps[(ps.gate == "TSS_ge_3") & (ps.contrast == "3plus_vs_1")]
    rows = []
    for gsm, g in ps.groupby("gsm"):
        meta = pd.read_csv(a.cache / (gsm + "_meta.tsv"), sep="\t").set_index("barcode")
        x = sparse.load_npz(a.cache / (gsm + "_atac.npz"))
        hi = meta.index.get_indexer(g.high_barcode)
        lo = meta.index.get_indexer(g.low_barcode)
        for name, members in sets.items():
            if not members:
                continue
            indices = [idx[p] for p in sorted(members)]
            high = np.asarray(x[indices][:, hi].mean(axis=0)).ravel()
            low = np.asarray(x[indices][:, lo].mean(axis=0)).ravel()
            rows.extend(
                dict(set_id=name, gsm=gsm, pair=i, high=h, low=l)
                for i, (h, l) in enumerate(zip(high, low))
            )
    scores = pd.DataFrame(rows)
    scores.to_csv(a.out / "matched_set_scores.tsv", sep="\t", index=False)
    summaries = []
    for name, g in scores.groupby("set_id"):
        delta = g.high - g.low
        t = stats.ttest_1samp(delta, 0)
        signs = g.assign(delta=delta).groupby("gsm").delta.mean()
        summaries.append(
            dict(
                set_id=name,
                n_peaks=len(sets[name]),
                n_pairs=len(g),
                FC=g.high.mean() / g.low.mean(),
                delta_pp=100 * delta.mean(),
                p=t.pvalue,
                positive_libraries=int(sum(signs > 0)),
            )
        )
    summary = pd.DataFrame(summaries)
    summary["q_all_testable_sets"] = multipletests(summary.p, method="fdr_bh")[1]
    summary.to_csv(a.out / "module_effects.tsv", sep="\t", index=False)
    (a.out / "reference_manifest.json").write_text(
        json.dumps(
            {
                "hmm_url": "https://hgdownload.soe.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeBroadHmm/wgEncodeBroadHmmHsmmHMM.bed.gz",
                "hmm_sha256": hashlib.sha256(a.hmm.read_bytes()).hexdigest(),
                "assembly": "hg19",
                "strong_states": ["4_Strong_Enhancer", "5_Strong_Enhancer"],
                "overlap_rule": "Any positive overlap; fractions reported; no target effect filter",
                "MYOD_ChIP_sources": manifest,
            },
            indent=2,
        )
        + "\n"
    )
    print(summary.to_string(index=False))
    focus = pd.concat(pieces)
    print(
        focus[(focus.peak == "chr11:17649919-17650798")][
            ["set_id", "fold_open", "p_pair_binomial", "q_within_external_set"]
        ].to_string(index=False)
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    p.add_argument("--hmm", type=Path, required=True)
    p.add_argument("--chip", type=Path, required=True)
    p.add_argument("--cache", type=Path, required=True)
    p.add_argument("--out", type=Path, required=True)
    main(p.parse_args())
