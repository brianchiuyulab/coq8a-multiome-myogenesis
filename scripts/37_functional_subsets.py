"""Freeze GO functional subsets and test externally dynamic myogenic peaks.

Run freeze before test. Membership uses external GO annotations, the existing
221-gene universe, fixed 100-kb candidate neighborhoods, and external temporal
labels only. COQ8A effect sizes and P values are not membership inputs.
"""

import argparse
import hashlib
import json
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from scipy import sparse
from statsmodels.stats.multitest import multipletests

from multiome_core import h5_for, read_barcodes, paired_stats


SETS = {
    "Differentiation": "GOBP_MYOBLAST_DIFFERENTIATION",
    "Fusion": "GOBP_MYOBLAST_FUSION",
    "Positive_regulation": "GOBP_POSITIVE_REGULATION_OF_MUSCLE_CELL_DIFFERENTIATION",
    "Negative_regulation": "GOBP_NEGATIVE_REGULATION_OF_MUSCLE_CELL_DIFFERENTIATION",
}
PHASES = {
    "Early": "Early_by24h",
    "Middle": "Middle_24to48h",
    "Late": "Late_after48h",
    "Closing": "Closing",
}


def freeze(root, out):
    sources = {}
    for name in ["muscle_subprogrammes.json", "functional_regulation_supplement.json"]:
        sources.update(json.loads((root / "reference" / name).read_text()))
    universe = set(
        pd.read_csv(root / "results/tables/gene_search_space.tsv", sep="\t").gene
    )
    candidates = pd.read_csv(
        root / "results/tables/candidate_peak_gene.tsv.gz", sep="\t"
    )
    candidates = candidates[candidates.n_libraries.eq(4)]
    assert candidates.peak.nunique() == 5097
    temporal = pd.read_csv(
        root / "results/temporal/temporal_set_memberships.tsv.gz", sep="\t"
    )
    phases = pd.concat(
        [
            temporal.loc[temporal.set_id.eq("P2_O50_" + suffix), ["peak"]].assign(
                phase=phase
            )
            for phase, suffix in PHASES.items()
        ],
        ignore_index=True,
    )
    assert len(phases) == 410 and phases.peak.is_unique
    genes, members, provenance = [], [], {}
    for label, name in SETS.items():
        source = sources[name]
        selected = sorted(universe & set(source["geneSymbols"]))
        genes.extend(dict(function=label, gene=gene, source=name) for gene in selected)
        linked = candidates[candidates.gene.isin(selected)].merge(phases, on="peak")
        linked["function"] = label
        linked["direction"] = np.where(linked.phase.eq("Closing"), "Closing", "Opening")
        linked["set_id"] = linked.function + "__" + linked.direction
        members.append(linked)
        provenance[label] = {
            "name": name,
            "genes_in_221": selected,
            "source_url": source["download_url"],
            "source_sha256": source["download_sha256"],
            "exact_source": source.get("exactSource"),
            "retrieved_utc": source.get("retrieved_utc"),
        }
    genes = pd.DataFrame(genes)
    members = pd.concat(members, ignore_index=True)
    genes.to_csv(out / "functional_genes.tsv", sep="\t", index=False)
    members.to_csv(out / "functional_peak_memberships.tsv", sep="\t", index=False)
    counts = []
    for label in SETS:
        gs = set(genes.loc[genes.function.eq(label), "gene"])
        row = dict(
            function=label,
            n_genes=len(gs),
            genes=";".join(sorted(gs)),
            n_common_peaks=candidates[candidates.gene.isin(gs)].peak.nunique(),
        )
        sub = members[members.function.eq(label)]
        for phase in PHASES:
            row["n_" + phase.lower()] = sub[sub.phase.eq(phase)].peak.nunique()
        counts.append(row)
    pd.DataFrame(counts).to_csv(out / "functional_counts.tsv", sep="\t", index=False)
    manifest = dict(
        sources=provenance,
        status="exploratory external-annotation subset analysis",
        primary_gate="TSS_ge_3",
        primary_contrast="3plus_vs_1",
        module_family="all nonempty function-by-opening/closing modules within each QC/contrast",
        peak_family="union of distinct peaks in all four functions, opening and closing together",
        membership_hash_encoding="UTF-8 with LF newlines",
        membership_sha256=hashlib.sha256(
            (out / "functional_peak_memberships.tsv")
            .read_text(encoding="utf-8")
            .encode("utf-8")
        ).hexdigest(),
    )
    (out / "analysis_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(pd.DataFrame(counts).drop(columns="genes").to_string(index=False), flush=True)


def summarize(frame):
    high, low = frame.high_open.to_numpy(), frame.low_open.to_numpy()
    stats = paired_stats(high - low)
    return dict(
        n_pairs=len(frame),
        high_open_pct=high.mean() * 100,
        low_open_pct=low.mean() * 100,
        fold_open=high.mean() / low.mean() if low.mean() else np.nan,
        delta_pp=stats["difference"] * 100,
        ci_low_pp=stats["ci_low"] * 100,
        ci_high_pp=stats["ci_high"] * 100,
        p_pair=stats["p_pair"],
    )


def test(root, out, h5_root):
    manifest = json.loads((out / "analysis_manifest.json").read_text())
    assert (
        hashlib.sha256(
            (out / "functional_peak_memberships.tsv")
            .read_text(encoding="utf-8")
            .encode("utf-8")
        ).hexdigest()
        == manifest["membership_sha256"]
    )
    members = pd.read_csv(out / "functional_peak_memberships.tsv", sep="\t")
    unique = members[["set_id", "peak"]].drop_duplicates()
    definitions = unique.groupby("set_id").size().rename("n_peaks").reset_index()
    ids = definitions.set_id.tolist()
    pairs = pd.read_csv(root / "results/tables/matched_pairs.tsv.gz", sep="\t")
    mapping = pd.read_csv(
        root / "results/tables/consensus_peak_map.tsv.gz", sep="\t"
    ).set_index("peak")
    scoreframes = []
    for gsm, ps in pairs.groupby("gsm"):
        barcodes = sorted(set(ps.high_barcode) | set(ps.low_barcode))
        with h5py.File(h5_for(h5_root, gsm)) as h5:
            _, matrix, _, peaks = read_barcodes(h5, barcodes)
        bidx = {b: i for i, b in enumerate(barcodes)}
        pidx = {str(p): i for i, p in enumerate(peaks)}
        target = (
            unique.peak
            if gsm == "GSM6339597"
            else unique.peak.map(mapping[gsm + "_peak"])
        )
        rows = target.map(pidx)
        cols = unique.set_id.map({k: i for i, k in enumerate(ids)})
        assert rows.notna().all()
        denom = unique.groupby("set_id").size()
        weights = sparse.csr_matrix(
            (1 / unique.set_id.map(denom).to_numpy(), (rows.astype(int), cols)),
            shape=(len(peaks), len(ids)),
        )
        values = (matrix @ weights).toarray()
        for r in ps.itertuples():
            block = pd.DataFrame(
                dict(
                    set_id=ids,
                    high_open=values[bidx[r.high_barcode]],
                    low_open=values[bidx[r.low_barcode]],
                )
            )
            for key in ["gsm", "source", "stage", "gate", "contrast", "pair"]:
                block[key] = getattr(r, key)
            scoreframes.append(block)
        print(gsm, "complete", flush=True)
    scores = pd.concat(scoreframes, ignore_index=True)
    scores.to_csv(out / "matched_functional_scores.tsv.gz", sep="\t", index=False)
    outputs = []
    libraries = []
    for keys, frame in scores.groupby(["gate", "contrast", "set_id", "gsm", "stage"]):
        libraries.append(
            dict(zip(["gate", "contrast", "set_id", "gsm", "stage"], keys))
            | summarize(frame)
        )
    libs = pd.DataFrame(libraries)
    libs.to_csv(out / "functional_effects_by_library.tsv", sep="\t", index=False)
    for keys, frame in scores.groupby(["gate", "contrast", "set_id"]):
        row = dict(zip(["gate", "contrast", "set_id"], keys)) | summarize(frame)
        ll = libs[
            (libs.gate == keys[0])
            & (libs.contrast == keys[1])
            & (libs.set_id == keys[2])
        ]
        row["positive_libraries"] = int((ll.delta_pp > 1e-10).sum())
        row["negative_libraries"] = int((ll.delta_pp < -1e-10).sum())
        outputs.append(row)
    result = pd.DataFrame(outputs).merge(definitions, on="set_id")
    for _, frame in result.groupby(["gate", "contrast"]):
        result.loc[frame.index, "q_modules"] = multipletests(
            frame.p_pair, method="fdr_bh"
        )[1]
    result.to_csv(out / "functional_effects.tsv", sep="\t", index=False)
    effects = pd.read_csv(
        root / "results/temporal/bidirectional_peak_effects.tsv", sep="\t"
    )
    effects = effects[effects.peak.isin(unique.peak)].copy()
    for _, frame in effects.groupby(["gate", "contrast"]):
        assert len(frame) == unique.peak.nunique()
        effects.loc[frame.index, "q_functional_union"] = multipletests(
            frame.p_exact, method="fdr_bh"
        )[1]
    labels = members.groupby("peak").function.agg(lambda x: ";".join(sorted(set(x))))
    effects["functions"] = effects.peak.map(labels)
    effects.to_csv(out / "functional_peak_effects.tsv", sep="\t", index=False)
    # Reconcile each aggregate FC against previously computed single-peak counts.
    for r in result.itertuples():
        selected = set(unique.loc[unique.set_id.eq(r.set_id), "peak"])
        subset = effects[
            (effects.gate == r.gate)
            & (effects.contrast == r.contrast)
            & effects.peak.isin(selected)
        ]
        assert len(subset) == r.n_peaks
        assert np.isclose(r.fold_open, subset.high_open.sum() / subset.low_open.sum())
    main = result[(result.gate == "TSS_ge_3") & (result.contrast == "3plus_vs_1")]
    print(
        main[["set_id", "n_peaks", "fold_open", "p_pair", "q_modules"]].to_string(
            index=False
        )
    )
    print("Distinct functional peaks:", unique.peak.nunique())


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode", choices=["freeze", "test"])
    parser.add_argument(
        "--root", type=Path, default=Path(__file__).resolve().parents[1]
    )
    parser.add_argument("--h5-root", type=Path)
    args = parser.parse_args()
    destination = args.root / "results/functional"
    destination.mkdir(parents=True, exist_ok=True)
    if args.mode == "freeze":
        freeze(args.root, destination)
    else:
        if args.h5_root is None:
            parser.error("test requires --h5-root")
        test(args.root, destination, args.h5_root)
