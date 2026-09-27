"""Rebuild the functional scope and test undifferentiated cultures separately.

The primary comparison is TSS >=3, COQ8A >=2 versus exactly 1 UMI in
undifferentiated libraries. Day 7 and >=3 versus 1 are separate comparisons.
"""

import argparse
import hashlib
import importlib.util
import json
from pathlib import Path
import subprocess
import sys

import numpy as np
import pandas as pd
from scipy import sparse, stats
from statsmodels.stats.multitest import multipletests
from threadpoolctl import threadpool_limits


def save(x, path):
    x.to_csv(
        path,
        sep="\t",
        index=False,
        na_rep="NA",
        compression=(
            {"method": "gzip", "mtime": 0} if str(path).endswith(".gz") else None
        ),
    )


def bh(p):
    return multipletests(p, method="fdr_bh")[1]


def main(a):
    root, work = a.root.resolve(), a.work.resolve()
    out = root / "results/task1_stage_specific"
    out.mkdir(parents=True, exist_ok=True)
    universe = set(
        pd.read_csv(root / "reference/myogenesis_221_gene_sources.tsv", sep="\t").gene
    )
    external = json.loads((root / "reference/muscle_subprogrammes.json").read_text())
    diff = universe & set(external["GOBP_MYOBLAST_DIFFERENTIATION"]["geneSymbols"])
    fusion = universe & set(external["GOBP_MYOBLAST_FUSION"]["geneSymbols"])
    genes = sorted(diff | fusion)
    assert (len(universe), len(diff), len(fusion), len(genes)) == (221, 17, 12, 25)
    save(
        pd.DataFrame(
            {
                "gene": genes,
                "differentiation": [g in diff for g in genes],
                "fusion": [g in fusion for g in genes],
            }
        ),
        out / "genes.tsv",
    )
    near = pd.read_csv(root / "results/tables/candidate_peak_gene.tsv.gz", sep="\t")
    common = near[near.n_libraries.eq(4) & near.nearest_tss_bp.le(100000)]
    member = common[common.gene.isin(genes)].copy()
    names = sorted(member.peak.unique())
    assert common.peak.nunique() == 5097 and len(names) == 565
    save(member, out / "peak_membership.tsv")
    group = {g: pd.Index(names).get_indexer(d.peak) for g, d in member.groupby("gene")}

    # Reapply the archived QC and matching algorithm to all libraries independently.
    qcwork = work / "stage_specific_qc"
    subprocess.run(
        [
            sys.executable,
            str(root / "scripts/02_qc_and_matching.py"),
            "--nuclei",
            str(root / "results/tables/nuclei.tsv.gz"),
            "--fragment-qc",
            str(root / "reference/fragment_qc"),
            "--out",
            str(qcwork),
        ],
        check=True,
    )
    allpairs = pd.read_csv(qcwork / "matched_pairs.tsv.gz", sep="\t")
    pairs = allpairs[allpairs.gate.eq("TSS_ge_3")].copy()
    previous = pd.read_csv(root / "results/tables/matched_pairs.tsv.gz", sep="\t")
    keys = ["gsm", "gate", "contrast", "high_barcode", "low_barcode"]
    assert set(map(tuple, allpairs[keys].to_numpy())) == set(
        map(tuple, previous[keys].to_numpy())
    )
    save(pairs, out / "matched_pairs.tsv.gz")
    save(pd.read_csv(qcwork / "qc_by_library.tsv", sep="\t"), out / "qc_by_library.tsv")
    inventory = pd.read_csv(qcwork / "qc_nuclei.tsv.gz", sep="\t")
    inventory = inventory[inventory.pass_tss3]
    inv = []
    for gsm, d in inventory.groupby("gsm"):
        inv.append(
            dict(
                gsm=gsm,
                stage=d.stage.iloc[0],
                qc_nuclei=len(d),
                zero=int(d.COQ8A_umi.eq(0).sum()),
                low1=int(d.COQ8A_umi.eq(1).sum()),
                high2=int(d.COQ8A_umi.ge(2).sum()),
                high3=int(d.COQ8A_umi.ge(3).sum()),
            )
        )
    save(pd.DataFrame(inv), out / "group_inventory.tsv")

    peak_index = pd.Index((work / "peaks.txt").read_text().splitlines()).get_indexer(
        names
    )
    rna_index = pd.Index((work / "genes.txt").read_text().splitlines()).get_indexer(
        genes
    )
    assert min(peak_index) >= 0 and min(rna_index) >= 0
    data = {}
    for gsm in pairs.gsm.unique():
        meta = pd.read_csv(work / f"{gsm}_meta.tsv", sep="\t").set_index("barcode")
        atac = sparse.load_npz(work / f"{gsm}_atac.npz")[peak_index].astype(np.float32)
        assert np.isin(atac.data, [0, 1]).all()
        rna = sparse.load_npz(work / f"{gsm}_rna.npz")[rna_index]
        data[gsm] = (meta, atac, rna)
    regional = []
    peaks = []
    libraries = []
    rnas = []
    balances = []
    modules = []
    rng = np.random.default_rng(a.seed)
    for stage in ["stem", "differentiated"]:
        for contrast in ["2plus_vs_1", "3plus_vs_1"]:
            ps = pairs[pairs.stage.eq(stage) & pairs.contrast.eq(contrast)]
            high = []
            low = []
            rh = []
            rl = []
            for gsm, d in ps.groupby("gsm"):
                meta, atac, rna = data[gsm]
                hi = meta.index.get_indexer(d.high_barcode)
                lo = meta.index.get_indexer(d.low_barcode)
                assert min(hi) >= 0 and min(lo) >= 0
                h = atac[:, hi].toarray().T
                l = atac[:, lo].toarray().T
                high.append(h)
                low.append(l)
                hr = (
                    rna[:, hi].toarray()
                    / meta.iloc[hi].total_rna_umi.to_numpy()
                    * 10000
                )
                lr = (
                    rna[:, lo].toarray()
                    / meta.iloc[lo].total_rna_umi.to_numpy()
                    * 10000
                )
                rh.append(hr)
                rl.append(lr)
                for j, p in enumerate(names):
                    libraries.append(
                        dict(
                            stage=stage,
                            contrast=contrast,
                            gsm=gsm,
                            peak=p,
                            n_pairs=len(d),
                            high_open=int(h[:, j].sum()),
                            low_open=int(l[:, j].sum()),
                        )
                    )
                for col in [
                    "total_rna_umi",
                    "total_open_peaks",
                    "myogenesis_score",
                    "percent_mito",
                ]:
                    balances.append(
                        dict(
                            stage=stage,
                            contrast=contrast,
                            gsm=gsm,
                            metric=col,
                            n_pairs=len(d),
                            high_mean=meta.iloc[hi][col].mean(),
                            low_mean=meta.iloc[lo][col].mean(),
                        )
                    )
                for j, g in enumerate(genes):
                    rnas.append(
                        dict(
                            stage=stage,
                            contrast=contrast,
                            gsm=gsm,
                            gene=g,
                            n_pairs=len(d),
                            mean_high=hr[j].mean(),
                            mean_low=lr[j].mean(),
                            FC=(
                                hr[j].mean() / lr[j].mean()
                                if lr[j].mean() > 0
                                else np.nan
                            ),
                        )
                    )
            h = np.vstack(high)
            l = np.vstack(low)
            delta = h - l
            hs = h.sum(0)
            ls = l.sum(0)
            hp = ((h > 0) & (l == 0)).sum(0)
            lp = ((l > 0) & (h == 0)).sum(0)
            pv = np.array(
                [
                    stats.binomtest(int(x), int(x + y)).pvalue if x + y else 1.0
                    for x, y in zip(hp, lp)
                ]
            )
            pq = bh(pv)
            for j, p in enumerate(names):
                peaks.append(
                    dict(
                        stage=stage,
                        contrast=contrast,
                        peak=p,
                        n_pairs=len(ps),
                        high_open=int(hs[j]),
                        low_open=int(ls[j]),
                        high_fraction=hs[j] / len(ps),
                        low_fraction=ls[j] / len(ps),
                        FC=hs[j] / ls[j] if ls[j] > 0 else np.nan,
                        high_only=int(hp[j]),
                        low_only=int(lp[j]),
                        p=pv[j],
                        q565=pq[j],
                    )
                )
            mh = h.mean(1)
            ml = l.mean(1)
            modules.append(
                dict(
                    stage=stage,
                    contrast=contrast,
                    n_pairs=len(ps),
                    high=mh.mean(),
                    low=ml.mean(),
                    FC=mh.mean() / ml.mean(),
                    p=stats.ttest_rel(mh, ml).pvalue,
                )
            )
            hr = np.hstack(rh)
            lr = np.hstack(rl)
            rp = np.nan_to_num(
                stats.ttest_rel(np.log1p(hr), np.log1p(lr), axis=1).pvalue, nan=1.0
            )
            rq = bh(rp)
            for j, g in enumerate(genes):
                rnas.append(
                    dict(
                        stage=stage,
                        contrast=contrast,
                        gsm="combined",
                        gene=g,
                        n_pairs=len(ps),
                        mean_high=hr[j].mean(),
                        mean_low=lr[j].mean(),
                        FC=hr[j].mean() / lr[j].mean() if lr[j].mean() > 0 else np.nan,
                        p=rp[j],
                        q25=rq[j],
                    )
                )
            normalized = delta / np.sqrt(np.maximum((delta**2).sum(0), 1))
            observed = normalized.sum(0)
            modes = {"opening": lambda x: x, "two_sided": np.abs}
            obs = {
                m: np.array([max(0, fn(observed[group[g]]).max()) for g in genes])
                for m, fn in modes.items()
            }
            counters = {m: np.zeros(len(genes), dtype=np.int64) for m in modes}
            with threadpool_limits(limits=2):
                for start in range(0, a.permutations, 1000):
                    batch = min(1000, a.permutations - start)
                    signs = (
                        rng.integers(0, 2, size=(batch, len(ps)), dtype=np.int8) * 2 - 1
                    ).astype(np.float32)
                    sim = signs @ normalized
                    for m, fn in modes.items():
                        v = fn(sim)
                        for j, g in enumerate(genes):
                            counters[m][j] += np.count_nonzero(
                                np.maximum(0, v[:, group[g]].max(1)) >= obs[m][j] - 1e-5
                            )
                    if (start + batch) % 250000 == 0:
                        print(
                            stage, contrast, start + batch, "permutations", flush=True
                        )
            for m, fn in modes.items():
                pv = (counters[m] + 1) / (a.permutations + 1)
                qv = bh(pv)
                for j, g in enumerate(genes):
                    ids = group[g]
                    focal = ids[np.argmax(fn(observed[ids]))]
                    regional.append(
                        dict(
                            stage=stage,
                            contrast=contrast,
                            mode=m,
                            gene=g,
                            n_peaks=len(ids),
                            n_pairs=len(ps),
                            statistic=obs[m][j],
                            p=pv[j],
                            q25=qv[j],
                            top_peak=names[focal],
                            top_peak_FC=(
                                hs[focal] / ls[focal] if ls[focal] > 0 else np.nan
                            ),
                            exceedances=int(counters[m][j]),
                            permutations=a.permutations,
                            seed=a.seed,
                        )
                    )
            save(pd.DataFrame(regional), out / "regional_tests.tsv")
            print(
                pd.DataFrame(regional)
                .query('stage==@stage and contrast==@contrast and mode=="opening"')
                .sort_values("p")[["gene", "p", "q25", "top_peak_FC"]]
                .head(5)
                .to_string(index=False),
                flush=True,
            )
    save(pd.DataFrame(peaks), out / "peak_tests.tsv.gz")
    save(pd.DataFrame(libraries), out / "peak_by_library.tsv.gz")
    save(pd.DataFrame(rnas), out / "RNA_supplement.tsv")
    save(pd.DataFrame(balances), out / "depth_and_state_balance.tsv")
    save(pd.DataFrame(modules), out / "whole_scope_summary.tsv")

    # Recombine existing within-library fits by stage, using two eligible sources.
    per = pd.read_csv(
        root / "results/unstratified/correlations_by_library.tsv.gz", sep="\t"
    )
    per = per[
        per.peak.isin(names) & per.model.isin(["technical_coq", "technical_state_coq"])
    ].copy()
    per["stage"] = per.gsm.map(
        {
            "GSM6339597": "stem",
            "GSM6339601": "stem",
            "GSM6339599": "differentiated",
            "GSM6339603": "differentiated",
        }
    )
    links = []
    for detection in [0.05, 0.01]:
        eligible = per[
            (per.peak_detected / per.n > detection)
            & (per.gene_detected / per.n > detection)
            & per.r.notna()
        ]
        for (stage, model, p, g), d in eligible.groupby(
            ["stage", "model", "peak", "gene"], sort=False
        ):
            if len(d) != 2:
                continue
            w = (d.n - d["rank"] - 2).to_numpy()
            z = np.arctanh(np.clip(d.r.to_numpy(), -0.999999, 0.999999))
            zm = np.average(z, weights=w)
            links.append(
                dict(
                    stage=stage,
                    model=model,
                    peak=p,
                    gene=g,
                    detection=detection,
                    n_libraries=2,
                    r=np.tanh(zm),
                    p=2 * stats.norm.sf(abs(zm) * np.sqrt(w.sum())),
                    r_min=d.r.min(),
                    r_max=d.r.max(),
                    gene_type=d.gene_type.iloc[0],
                    distance_bp=d.distance_bp.iloc[0],
                )
            )
    links = pd.DataFrame(links)
    for _, d in links.groupby(["stage", "model", "detection"]):
        links.loc[d.index, "q_links"] = bh(d.p)
        links.loc[d.index, "n_tests"] = len(d)
    save(links, out / "stage_specific_local_RNA_links.tsv.gz")
    audit = {
        "scope_genes": len(genes),
        "distinct_peaks": len(names),
        "primary_stage": "stem",
        "primary_contrast": "2plus_vs_1",
        "TSS": 3,
        "permutations": a.permutations,
        "seed": a.seed,
        "matching_reproduced": True,
        "stage_is_culture_condition_not_cell_annotation": True,
        "candidate_catalog": "four-library common; unchanged from previous scope",
        "input_sha256": {
            str(p.relative_to(root)): hashlib.sha256(p.read_bytes()).hexdigest()
            for p in [
                root / "reference/myogenesis_221_gene_sources.tsv",
                root / "reference/muscle_subprogrammes.json",
                root / "results/tables/candidate_peak_gene.tsv.gz",
                root / "results/tables/matched_pairs.tsv.gz",
            ]
        },
    }
    (out / "provenance.json").write_text(json.dumps(audit, indent=2) + "\n")


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    p.add_argument("--work", type=Path, required=True)
    p.add_argument("--permutations", type=int, default=2000000)
    p.add_argument("--seed", type=int, default=20260928)
    main(p.parse_args())
