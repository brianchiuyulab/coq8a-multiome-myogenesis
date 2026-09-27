"""Fixed Day-0 sensitivity grid for local opening and closing associations.

All settings, including unavailable ones, are retained. Marker-positive subsets
are sensitivity populations, not independently validated cell-type annotations.
"""

import argparse
import hashlib
import importlib.util
import itertools
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse, stats
from sklearn.neighbors import NearestNeighbors
from statsmodels.stats.multitest import multipletests
from threadpoolctl import threadpool_limits

SAMPLES = ["GSM6339597", "GSM6339601"]
CONTRASTS = [
    "2plus_vs_1",
    "3plus_vs_1",
    "4plus_vs_1",
    "detected_vs_zero",
    "2plus_vs_zero",
    "2plus_vs_le1",
    "positive_CP10k_Q25",
]


def save(d, p):
    d.to_csv(
        p,
        sep="\t",
        index=False,
        na_rep="NA",
        compression={"method": "gzip", "mtime": 0} if str(p).endswith(".gz") else None,
    )


def load_module(path):
    spec = importlib.util.spec_from_file_location("extract_module", path)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def prepare(root, h5, work, names, scope):
    work.mkdir(parents=True, exist_ok=True)
    genes = sorted(
        set(
            pd.read_csv(
                root / "reference/myogenesis_221_gene_sources.tsv", sep="\t"
            ).gene
        )
        | {"COQ8A", "PAX7", "MYOD1", "MYOG", "MYF5", "PTPRC", "PECAM1", "PDGFRA"}
    )
    extractor = load_module(root / "scripts/44_expand_middle_links.py")
    mapping = pd.read_csv(
        root / "results/tables/consensus_peak_map.tsv.gz", sep="\t"
    ).set_index("peak")
    qc = pd.read_csv(root / "results/tables/qc_nuclei.tsv.gz", sep="\t")
    scoregenes = sorted(
        set(genes) - set(scope) - {"COQ8A", "PAX7", "PTPRC", "PECAM1", "PDGFRA"}
    )
    gi = {g: i for i, g in enumerate(genes)}
    for gsm in SAMPLES:
        meta = qc[qc.gsm.eq(gsm) & qc.pass_tss2].copy().reset_index(drop=True)
        local = [
            p if gsm == SAMPLES[0] else mapping.loc[p, gsm + "_peak"] for p in names
        ]
        mat, mito = extractor.extract(
            next(h5.glob(gsm + "_*filtered_feature_bc_matrix.h5")),
            genes,
            local,
            meta.barcode.tolist(),
        )
        rna = mat[: len(genes)].astype(float)
        atac = mat[len(genes) :].astype(np.float32)
        atac.data[:] = 1
        coq = rna[gi["COQ8A"]].toarray().ravel()
        assert np.array_equal(coq, meta.COQ8A_umi)
        meta["coq_CP10k"] = coq / meta.total_rna_umi * 10000
        meta["myogenic_detected"] = (
            rna[[gi[g] for g in ["PAX7", "MYF5", "MYOD1", "MYOG"]]].sum(0).A1 > 0
        )
        meta["PAX7_detected"] = rna[gi["PAX7"]].toarray().ravel() > 0
        meta["state_score"] = np.log1p(
            rna[[gi[g] for g in scoregenes]].toarray()
            / meta.total_rna_umi.to_numpy()
            * 10000
        ).mean(0)
        meta["percent_mito"] = mito / meta.total_rna_umi * 100
        save(meta, work / f"{gsm}_meta.tsv")
        sparse.save_npz(work / f"{gsm}_atac.npz", atac)
        sparse.save_npz(work / f"{gsm}_rna.npz", rna)
    (work / "genes.txt").write_text("\n".join(genes) + "\n")
    (work / "peaks.txt").write_text("\n".join(names) + "\n")


def match(meta, c):
    mask = meta.tss_enrichment.ge(c["TSS"])
    if c["population"] != "all_culture":
        mask &= meta[c["population"]]
    d = meta[mask].copy()
    coq = d.COQ8A_umi
    definition = c["contrast"]
    if definition == "positive_CP10k_Q25":
        pos = d[d.COQ8A_umi > 0].coq_CP10k
        if not len(pos):
            return [], dict(n_eligible=len(d), high=0, low=0)
        lowq, highq = pos.quantile([0.25, 0.75])
        hm = coq.gt(0) & d.coq_CP10k.ge(highq)
        lm = coq.gt(0) & d.coq_CP10k.le(lowq)
        if (hm & lm).any():
            return [], dict(n_eligible=len(d), high=0, low=0)
    elif definition == "detected_vs_zero":
        hm = coq.gt(0)
        lm = coq.eq(0)
    else:
        cut = int(definition[0])
        hm = coq.ge(cut)
        lm = (
            coq.eq(1)
            if definition.endswith("_1")
            else coq.eq(0) if definition.endswith("_zero") else coq.le(1)
        )
    hi = d[hm]
    lo = d[lm]
    support = dict(n_eligible=len(d), high=len(hi), low=len(lo))
    if not len(hi) or not len(lo):
        return [], support
    columns = ["total_rna_umi", "total_open_peaks"]
    hx = np.log1p(hi[columns].to_numpy())
    lx = np.log1p(lo[columns].to_numpy())
    limit = np.array([c["caliper"], c["caliper"]])
    if c["matching"] == "depth_state":
        sd = d.state_score.std()
        if not np.isfinite(sd) or sd == 0:
            sd = 1
        hx = np.column_stack([hx, hi.state_score.to_numpy() / sd])
        lx = np.column_stack([lx, lo.state_score.to_numpy() / sd])
        limit = np.r_[limit, 0.25]
        search_h = hx / limit
        search_l = lx / limit
    else:
        search_h = hx
        search_l = lx
    dist, neighbors = (
        NearestNeighbors(n_neighbors=min(100, len(lo)), algorithm="kd_tree")
        .fit(search_l)
        .kneighbors(search_h)
    )
    used = set()
    pairs = []
    for i in np.argsort(-dist[:, 0]):
        for j in neighbors[i]:
            if int(j) not in used and (np.abs(hx[i] - lx[j]) <= limit).all():
                pairs.append((int(hi.index[i]), int(lo.index[j])))
                used.add(int(j))
                break
    return pairs, support


def permutation(delta, groups, genes, B, seed):
    x = delta / np.sqrt(np.maximum((delta**2).sum(0), 1))
    observed = x.sum(0)
    masks = {"opening": lambda z: z, "closing": lambda z: -z}
    obs = {
        m: np.array([max(0, fn(observed[groups[g]]).max()) for g in genes])
        for m, fn in masks.items()
    }
    counters = {m: np.zeros(len(genes), dtype=np.int64) for m in masks}
    rng = np.random.default_rng(seed)
    with threadpool_limits(limits=2):
        for start in range(0, B, 1000):
            k = min(1000, B - start)
            s = (
                rng.integers(0, 2, size=(k, len(delta)), dtype=np.int8) * 2 - 1
            ).astype(np.float32) @ x
            for m, fn in masks.items():
                v = fn(s)
                for j, g in enumerate(genes):
                    counters[m][j] += np.count_nonzero(
                        np.maximum(0, v[:, groups[g]].max(1)) >= obs[m][j] - 1e-5
                    )
    rows = []
    for m, fn in masks.items():
        p = (counters[m] + 1) / (B + 1)
        q = multipletests(p, method="fdr_bh")[1]
        for j, g in enumerate(genes):
            ids = groups[g]
            best = ids[np.argmax(fn(observed[ids]))]
            rows.append(
                dict(
                    mode=m,
                    gene=g,
                    statistic=obs[m][j],
                    p=p[j],
                    q25=q[j],
                    top_index=int(best),
                    exceedances=int(counters[m][j]),
                    permutations=B,
                    seed=seed,
                    MC_lower=stats.beta.ppf(
                        0.025, counters[m][j] + 1, B - counters[m][j] + 1
                    ),
                    MC_upper=stats.beta.ppf(
                        0.975, counters[m][j] + 1, B - counters[m][j] + 1
                    ),
                )
            )
    return rows


def main(a):
    root = a.root.resolve()
    work = a.work.resolve()
    out = root / "results/day0_sensitivity_grid"
    out.mkdir(parents=True, exist_ok=True)
    member = pd.read_csv(
        root / "results/task1_stage_specific/peak_membership.tsv", sep="\t"
    )
    names = sorted(member.peak.unique())
    genes = sorted(member.gene.unique())
    groups = {g: pd.Index(names).get_indexer(d.peak) for g, d in member.groupby("gene")}
    specs = []
    for i, (t, pop, coq, cal, mt) in enumerate(
        itertools.product(
            [2, 3],
            ["all_culture", "myogenic_detected", "PAX7_detected"],
            CONTRASTS,
            [0.05, 0.1, 0.2, 0.3],
            ["depth", "depth_state"],
        )
    ):
        specs.append(
            dict(
                config=f"D{i:03}",
                TSS=t,
                population=pop,
                contrast=coq,
                caliper=cal,
                matching=mt,
            )
        )
    # Freeze the complete grid before evaluating any target association.
    save(pd.DataFrame(specs), out / "configurations.tsv")
    (out / "design.json").write_text(
        json.dumps(
            dict(
                configurations=len(specs),
                screen_permutations=a.permutations,
                refinement_permutations=2000000,
                refinement_rule="top 3 configs per direction with concordant top-peak library directions; retain fixed primary as well",
                seed=20260928,
                scope_genes=25,
                peaks=565,
                stage="undifferentiated",
                marker_subsets="sensitivity proxies, not validated MuSC annotation",
                state_score="mean log1p CP10k over broad myogenesis genes excluding all 25 scope genes",
                state_caliper_SD=0.25,
                correction="BH25 within direction/config; BH565 for individual peaks; no across-grid FDR",
            ),
            indent=2,
        )
        + "\n"
    )
    if not a.reuse:
        prepare(root, a.h5.resolve(), work, names, genes)
    data = {
        gsm: (
            pd.read_csv(work / f"{gsm}_meta.tsv", sep="\t"),
            sparse.load_npz(work / f"{gsm}_atac.npz").toarray().T,
        )
        for gsm in SAMPLES
    }
    regional = []
    peaks = []
    pair_rows = []
    support = []
    status = []
    balance = []
    pairs_by_config = {}
    for n, c in enumerate(specs):
        ds = []
        hs = []
        ls = []
        counts = []
        stored = {}
        for gsm, (meta, atac) in data.items():
            pairs, s = match(meta, c)
            support.append(dict(**c, gsm=gsm, n_pairs=len(pairs), **s))
            stored[gsm] = pairs
            if not pairs:
                continue
            hi, lo = np.array(pairs).T
            h = atac[hi]
            l = atac[lo]
            hs.append(h)
            ls.append(l)
            ds.append(h - l)
            counts.append((h.sum(0), l.sum(0)))
            for i, j in pairs:
                pair_rows.append(
                    dict(
                        config=c["config"],
                        gsm=gsm,
                        high_barcode=meta.iloc[i].barcode,
                        low_barcode=meta.iloc[j].barcode,
                    )
                )
            for col in [
                "total_rna_umi",
                "total_open_peaks",
                "state_score",
                "percent_mito",
            ]:
                balance.append(
                    dict(
                        config=c["config"],
                        gsm=gsm,
                        metric=col,
                        high_mean=meta.iloc[hi][col].mean(),
                        low_mean=meta.iloc[lo][col].mean(),
                    )
                )
        if len(ds) != 2:
            status.append(dict(**c, status="no_pairs_in_one_or_both_sources"))
            continue
        delta = np.vstack(ds)
        h = np.vstack(hs)
        l = np.vstack(ls)
        pairs_by_config[c["config"]] = stored
        high = h.sum(0)
        low = l.sum(0)
        ho = ((h > 0) & (l == 0)).sum(0)
        lo = ((l > 0) & (h == 0)).sum(0)
        p = np.array(
            [
                stats.binomtest(int(x), int(x + y)).pvalue if x + y else 1.0
                for x, y in zip(ho, lo)
            ]
        )
        q = multipletests(p, method="fdr_bh")[1]
        concord_up = (counts[0][0] > counts[0][1]) & (counts[1][0] > counts[1][1])
        concord_down = (counts[0][0] < counts[0][1]) & (counts[1][0] < counts[1][1])
        for j, name in enumerate(names):
            peaks.append(
                dict(
                    config=c["config"],
                    peak=name,
                    n_pairs=len(delta),
                    high=int(high[j]),
                    low=int(low[j]),
                    FC=high[j] / low[j] if low[j] > 0 else np.nan,
                    p=p[j],
                    q565=q[j],
                    high_only=int(ho[j]),
                    low_only=int(lo[j]),
                    both_up=bool(concord_up[j]),
                    both_down=bool(concord_down[j]),
                    source1_high=int(counts[0][0][j]),
                    source1_low=int(counts[0][1][j]),
                    source1_pairs=len(ds[0]),
                    source2_high=int(counts[1][0][j]),
                    source2_low=int(counts[1][1][j]),
                    source2_pairs=len(ds[1]),
                )
            )
        rows = permutation(delta, groups, genes, a.permutations, 20260928 + n)
        for row in rows:
            j = row.pop("top_index")
            row.update(c)
            row.update(
                n_pairs=len(delta),
                top_peak=names[j],
                top_peak_FC=high[j] / low[j] if low[j] > 0 else np.nan,
                top_peak_p=p[j],
                top_peak_q565=q[j],
                both_sources_same_direction=bool(
                    concord_up[j] if row["mode"] == "opening" else concord_down[j]
                ),
                n_peaks=len(groups[row["gene"]]),
            )
            regional.append(row)
        status.append(dict(**c, status="evaluated"))
        if n % 12 == 0:
            save(pd.DataFrame(regional), out / "regional_screen.tsv.gz")
            print(
                "completed",
                n + 1,
                "/",
                len(specs),
                c["config"],
                "pairs",
                len(delta),
                "min q",
                min(r["q25"] for r in rows),
                flush=True,
            )
    save(pd.DataFrame(regional), out / "regional_screen.tsv.gz")
    save(pd.DataFrame(peaks), out / "peak_screen.tsv.gz")
    save(pd.DataFrame(pair_rows), out / "matched_pairs.tsv.gz")
    save(pd.DataFrame(support), out / "support.tsv")
    save(pd.DataFrame(status), out / "status.tsv")
    save(pd.DataFrame(balance), out / "balance.tsv")
    reg = pd.DataFrame(regional)
    chosen = []
    for mode in ["opening", "closing"]:
        d = reg[reg["mode"].eq(mode) & reg.both_sources_same_direction].sort_values(
            ["q25", "p", "config", "gene"]
        )
        chosen.extend(d.drop_duplicates("config").head(3).config)
    primary = [
        c["config"]
        for c in specs
        if c["TSS"] == 3
        and c["population"] == "all_culture"
        and c["contrast"] == "2plus_vs_1"
        and c["caliper"] == 0.1
        and c["matching"] == "depth"
    ][0]
    chosen = sorted(set(chosen + [primary]))
    (out / "refined_configurations.json").write_text(
        json.dumps(chosen, indent=2) + "\n"
    )
    refined = []
    for config in chosen:
        stored = pairs_by_config[config]
        ds = []
        for gsm, (meta, atac) in data.items():
            hi, lo = np.array(stored[gsm]).T
            ds.append(atac[hi] - atac[lo])
        rows = permutation(
            np.vstack(ds), groups, genes, 2000000, 20260928 + int(config[1:])
        )
        c = next(c for c in specs if c["config"] == config)
        for row in rows:
            j = row.pop("top_index")
            row.update(c)
            row["top_peak"] = names[j]
            refined.append(row)
        save(pd.DataFrame(refined), out / "regional_refined.tsv")
        print("refined", config, "min q", min(r["q25"] for r in rows), flush=True)
    print("DONE", flush=True)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    p.add_argument("--work", type=Path, required=True)
    p.add_argument("--h5", type=Path, required=True)
    p.add_argument("--reuse", action="store_true")
    p.add_argument("--permutations", type=int, default=20000)
    main(p.parse_args())
