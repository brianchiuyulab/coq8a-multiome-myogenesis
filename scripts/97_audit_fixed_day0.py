"""Audit the fixed Day-0 analysis against deposited matrices and annotations.

Recalculates candidate membership, joint QC, pairing, every peak test, local
RNA models, and (optionally) all 25 regional permutation tests. The audit
does not change COQ8A definitions, eligibility thresholds, or test families.
"""

import argparse
import gzip
import hashlib
import importlib.util
import importlib.metadata
import json
import re
import sys
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from scipy import sparse, stats


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "results/fixed_day0_D206"
CHECKS = []


def read(path):
    return pd.read_csv(path, sep="\t")


def check(name, condition, detail):
    if not condition:
        raise AssertionError(f"{name}: {detail}")
    CHECKS.append({"check": name, "passed": True, "detail": detail})
    print(name, "PASS", flush=True)


def bh(values):
    p = np.asarray(values, dtype=float)
    order = np.argsort(p)
    adjusted = np.minimum.accumulate(
        (p[order] * len(p) / np.arange(1, len(p) + 1))[::-1]
    )[::-1]
    result = np.empty(len(p))
    result[order] = np.minimum(adjusted, 1)
    return result


def load_module(path):
    spec = importlib.util.spec_from_file_location(path.stem, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def parse_peak(peak):
    chrom, bounds = peak.split(":")
    start, end = map(int, bounds.split("-"))
    return chrom, start, end


def raw_extract(path, requested, barcodes):
    """Independent chunked extraction, including complete RNA/ATAC depths."""
    with h5py.File(path) as handle:
        matrix = handle["matrix"]
        names = np.char.decode(matrix["features/name"][:])
        kinds = matrix["features/feature_type"][:]
        bars = np.char.decode(matrix["barcodes"][:])
        frequency = pd.Series(names).value_counts()
        if not frequency.reindex(requested).eq(1).all():
            raise ValueError("A requested feature is absent or has an ambiguous name")
        name_index = {name: i for i, name in enumerate(names)}
        indices = np.array([name_index[name] for name in requested])
        if (indices < 0).any():
            raise ValueError("Requested feature missing from deposited matrix")
        lookup = np.full(len(names), -1, dtype=int)
        lookup[indices] = np.arange(len(requested))
        barcode_index = pd.Index(bars).get_indexer(barcodes)
        assert (barcode_index >= 0).all()
        column_lookup = np.full(len(bars), -1, dtype=int)
        column_lookup[barcode_index] = np.arange(len(barcodes))
        output = np.zeros((len(barcodes), len(requested)), dtype=np.float32)
        totals = np.zeros((len(barcodes), 2), dtype=np.int64)
        ptr = matrix["indptr"][:]
        for start in range(0, len(bars), 1024):
            end = min(start + 1024, len(bars))
            lo, hi = int(ptr[start]), int(ptr[end])
            features = matrix["indices"][lo:hi]
            values = matrix["data"][lo:hi]
            columns = np.repeat(np.arange(start, end), np.diff(ptr[start : end + 1]))
            mapped = column_lookup[columns]
            keep = mapped >= 0
            mapped, features, values = mapped[keep], features[keep], values[keep]
            for k, modality in enumerate([b"Gene Expression", b"Peaks"]):
                mask = kinds[features] == modality
                weights = values[mask] if k == 0 else (values[mask] > 0)
                totals[:, k] += np.bincount(
                    mapped[mask], weights=weights, minlength=len(barcodes)
                ).astype(np.int64)
            selected = lookup[features] >= 0
            np.add.at(
                output, (mapped[selected], lookup[features[selected]]), values[selected]
            )
    return output, totals


def annotations(gtf, universe, selected, valid_rna):
    sites = {}
    local = {}
    for peak in selected:
        chrom, start, end = parse_peak(peak)
        local[peak] = (chrom, (start + end) / 2)
    targets = {}
    with gzip.open(gtf, "rt") as stream:
        for line in stream:
            if line.startswith("#"):
                continue
            fields = line.rstrip().split("\t")
            if fields[2] != "transcript":
                continue
            attrs = dict(re.findall(r'(\w+) "([^"]+)"', fields[8]))
            gene = attrs["gene_name"]
            position = int(fields[3] if fields[6] == "+" else fields[4]) - 1
            if gene in universe and attrs["gene_type"] == "protein_coding":
                sites.setdefault(fields[0], set()).add((position, gene))
            if gene in valid_rna:
                for peak, (chrom, center) in local.items():
                    distance = abs(position - center)
                    if fields[0] == chrom and distance <= 500000:
                        key = (peak, gene)
                        targets[key] = min(targets.get(key, np.inf), distance)
    return sites, targets


def main(args):
    output = SOURCE / "audit"
    output.mkdir(exist_ok=True)
    universe = set(read(ROOT / "reference/myogenesis_221_gene_sources.tsv").gene)
    external = json.loads((ROOT / "reference/muscle_subprogrammes.json").read_text())
    diff = universe & set(external["GOBP_MYOBLAST_DIFFERENTIATION"]["geneSymbols"])
    fusion = universe & set(external["GOBP_MYOBLAST_FUSION"]["geneSymbols"])
    scope = sorted(diff | fusion)
    check(
        "external_gene_scope",
        (len(universe), len(diff), len(fusion), len(scope)) == (221, 17, 12, 25),
        "221 intersect (17 differentiation union 12 fusion) = 25",
    )
    peaks = read(SOURCE / "all_peaks.tsv").sort_values("peak").reset_index(drop=True)
    regions = read(SOURCE / "all_regions.tsv")
    pairs = read(SOURCE / "matched_pairs.tsv")
    candidate = read(SOURCE / "local_RNA_candidates.tsv")
    selected = sorted(candidate.peak.unique())
    genes = sorted(candidate.gene.unique())
    mapping = read(ROOT / "results/tables/consensus_peak_map.tsv.gz").set_index("peak")
    first = next(args.h5.glob("GSM6339597_*filtered_feature_bc_matrix.h5"))
    with h5py.File(first) as h:
        f = h["matrix/features"]
        rn = np.char.decode(f["name"][:][f["feature_type"][:] == b"Gene Expression"])
        frequency = pd.Series(rn).value_counts()
        valid = set(frequency[frequency == 1].index)
    sites, target_map = annotations(args.gtf, universe, selected, valid)
    near = []
    sorted_sites = {chrom: sorted(values) for chrom, values in sites.items()}
    site_positions = {
        chrom: np.array([v[0] for v in values])
        for chrom, values in sorted_sites.items()
    }
    for peak, row in mapping[
        mapping.coordinate_qc & mapping.n_libraries.ge(3)
    ].iterrows():
        if row.chrom not in sorted_sites:
            continue
        center = (row.start + row.end) // 2
        left = np.searchsorted(site_positions[row.chrom], center - 100000, side="left")
        right = np.searchsorted(
            site_positions[row.chrom], center + 100000, side="right"
        )
        for position, gene in sorted_sites[row.chrom][left:right]:
            distance = abs(position - (row.start + row.end) // 2)
            if distance <= 100000:
                near.append((peak, gene, distance, row.n_libraries))
    rebuilt = (
        pd.DataFrame(near, columns=["peak", "gene", "nearest_tss_bp", "n_libraries"])
        .sort_values("nearest_tss_bp")
        .drop_duplicates(["peak", "gene"])
    )
    old_near = read(ROOT / "results/tables/candidate_peak_gene.tsv.gz")
    check(
        "zero_based_TSS_membership",
        set(zip(rebuilt.peak, rebuilt.gene)) == set(zip(old_near.peak, old_near.gene)),
        "Rebuilding with 0-based TSS leaves every candidate peak-gene membership unchanged",
    )
    common = rebuilt[rebuilt.n_libraries.eq(4)]
    member = common[common.gene.isin(scope)].copy()
    check(
        "peak_scope",
        common.peak.nunique() == 5097
        and set(member.peak) == set(peaks.peak)
        and len(peaks) == 565,
        "5097 broad peaks; 565 distinct differentiation/fusion peaks",
    )
    member.to_csv(SOURCE / "peak_membership.tsv", sep="\t", index=False)
    if args.correct_tss_distances:
        rebuilt.to_csv(
            ROOT / "results/tables/candidate_peak_gene.tsv.gz",
            sep="\t",
            index=False,
            compression={"method": "gzip", "mtime": 0},
        )
        member.to_csv(
            ROOT / "results/task1_stage_specific/peak_membership.tsv",
            sep="\t",
            index=False,
        )
    check(
        "local_RNA_universe",
        set(target_map) == set(zip(candidate.peak, candidate.gene))
        and all(
            np.isclose(target_map[(r.peak, r.gene)], r.distance_bp)
            for r in candidate.itertuples()
        ),
        f"{len(candidate)} peak-RNA pairs rebuilt from all transcript TSSs within 500 kb",
    )
    canonical = {f"chr{i}" for i in range(1, 23)} | {"chrX", "chrY"}
    m = mapping.loc[peaks.peak]
    check(
        "peak_coordinate_QC",
        m.chrom.isin(canonical).all()
        and m.width.between(200, 2000).all()
        and not m.blacklisted.any()
        and m.n_libraries.eq(4).all(),
        "Canonical, 200-2000 bp, blacklist-free anchor peaks shared across 4 libraries",
    )
    for gsm in ["GSM6339599", "GSM6339601", "GSM6339603"]:
        overlap_ok = []
        for peak, row in m.iterrows():
            chrom, lo, hi = parse_peak(row[gsm + "_peak"])
            overlap = max(0, min(hi, row.end) - max(lo, row.start))
            overlap_ok.append(
                chrom == row.chrom
                and min(overlap / (hi - lo), overlap / row.width) >= 0.5
            )
        check(
            "peak_mapping_" + gsm,
            all(overlap_ok) and m[gsm + "_peak"].is_unique,
            "Reciprocal overlap >=50%; one-to-one",
        )
    inventory = read(ROOT / "results/tables/nuclei.tsv.gz")
    inventory = inventory[
        inventory.total_rna_umi.ge(500) & inventory.total_open_peaks.ge(500)
    ].copy()
    cap = inventory.groupby("gsm")[["total_rna_umi", "total_open_peaks"]].transform(
        lambda x: x.quantile(0.95)
    )
    inventory = inventory[
        (inventory.total_rna_umi <= cap.total_rna_umi)
        & (inventory.total_open_peaks <= cap.total_open_peaks)
    ]
    for kind in ["tss", "nucleosome", "blacklist"]:
        table = pd.concat(
            [
                read(p)
                for p in sorted(
                    (ROOT / "reference/fragment_qc").glob(
                        f"GSM*_fragment_{kind}_qc.tsv"
                    )
                )
            ]
        )
        inventory = inventory.merge(table, on=["gsm", "barcode"], validate="one_to_one")
    inventory["eligible"] = (
        inventory.nucleosome_signal.lt(4)
        & (inventory.atac_peak_region_fragments / inventory.atac_fragments).ge(0.25)
        & (inventory.blacklist_fragments / inventory.atac_fragments).lt(0.05)
        & inventory.tss_enrichment.ge(3)
    )
    grid = load_module(ROOT / "scripts/88_day0_sensitivity_grid.py")
    design = json.loads((SOURCE / "design.json").read_text())
    check(
        "fixed_design",
        design["TSS"] == 3
        and design["contrast"] == "2plus_vs_zero"
        and design["population"] == "all_culture"
        and design["caliper"] == 0.3
        and design["matching"] == "depth",
        "Same fixed definition for all 565 peaks and 25 regions",
    )
    cache_genes = (args.work / "genes.txt").read_text().splitlines()
    rna_genes = sorted(set(cache_genes) | set(genes))
    gi = {g: i for i, g in enumerate(rna_genes)}
    peak_index = pd.Index(peaks.peak)
    source_links = read(SOURCE / "links_by_source.tsv")
    high, low, rna_high, rna_low = [], [], [], []
    for gsm in ["GSM6339597", "GSM6339601"]:
        meta = read(args.work / f"{gsm}_meta.tsv")
        eligible = meta[meta.pass_tss3]
        audited = inventory[inventory.gsm.eq(gsm) & inventory.eligible]
        check(
            "joint_QC_" + gsm,
            set(eligible.barcode) == set(audited.barcode),
            f"{len(eligible)} eligible Day-0 nuclei; thresholds reapplied",
        )
        local = (
            list(peaks.peak)
            if gsm == "GSM6339597"
            else mapping.loc[peaks.peak, gsm + "_peak"].tolist()
        )
        path = next(args.h5.glob(gsm + "_*filtered_feature_bc_matrix.h5"))
        raw, totals = raw_extract(path, rna_genes + local, meta.barcode.tolist())
        check(
            "raw_depth_" + gsm,
            np.array_equal(totals[:, 0], meta.total_rna_umi)
            and np.array_equal(totals[:, 1], meta.total_open_peaks),
            "Complete raw RNA counts and detected-peak depths match metadata",
        )
        atac = (raw[:, len(rna_genes) :] > 0).astype(np.float32)
        check(
            "raw_cache_" + gsm,
            np.array_equal(
                atac, sparse.load_npz(args.work / f"{gsm}_atac.npz").toarray().T
            )
            and np.array_equal(
                raw[:, [gi[g] for g in cache_genes]],
                sparse.load_npz(args.work / f"{gsm}_rna.npz").toarray().T,
            ),
            "All 565 ATAC peaks and cached RNA features independently re-extracted",
        )
        pr = pairs[pairs.gsm.eq(gsm)]
        ix = pd.Index(meta.barcode)
        hi, lo = ix.get_indexer(pr.high_barcode), ix.get_indexer(pr.low_barcode)
        reproduced, _ = grid.match(meta, design)
        check(
            "pairing_" + gsm,
            reproduced == list(zip(hi, lo)) and len(set(hi) | set(lo)) == 2 * len(pr),
            f"All {len(pr)} pairs reproduced; no reused nucleus",
        )
        depth = np.log1p(meta[["total_rna_umi", "total_open_peaks"]].to_numpy())
        check(
            "pair_definition_" + gsm,
            (raw[hi, gi["COQ8A"]] >= 2).all()
            and (raw[lo, gi["COQ8A"]] == 0).all()
            and (np.abs(depth[hi] - depth[lo]) <= 0.3).all(),
            "Raw COQ8A definition and both depth calipers verified",
        )
        high.append(atac[hi])
        low.append(atac[lo])
        cp = raw[:, [gi[g] for g in genes]].astype(float) / totals[:, 0, None] * 10000
        rna_high.append(cp[hi])
        rna_low.append(cp[lo])
        for population, ids in [
            ("eligible_D0", np.flatnonzero(meta.pass_tss3)),
            ("matched", np.r_[hi, lo]),
        ]:
            for state in [False, True]:
                cov = np.column_stack(
                    [
                        np.ones(len(ids)),
                        depth[ids],
                        np.log1p(meta.coq_CP10k.to_numpy()[ids]),
                    ]
                )
                if state:
                    cov = np.column_stack([cov, meta.state_score.to_numpy()[ids]])
                q, _ = np.linalg.qr(cov, mode="reduced")
                x = atac[ids].astype(float)
                y = np.log1p(cp[ids])
                x -= q @ (q.T @ x)
                y -= q @ (q.T @ y)
                model = "depth_COQ_state" if state else "depth_COQ"
                saved = source_links[
                    (source_links.gsm == gsm)
                    & (source_links.population == population)
                    & (source_links.model == model)
                ].set_index(["peak", "gene"])
                values = []
                for row in candidate.itertuples():
                    xv, yv = (
                        x[:, peak_index.get_loc(row.peak)],
                        y[:, genes.index(row.gene)],
                    )
                    den = np.linalg.norm(xv) * np.linalg.norm(yv)
                    corr = np.dot(xv, yv) / den if den > 1e-10 else np.nan
                    values.append(corr)
                expected = saved.loc[
                    list(zip(candidate.peak, candidate.gene)), "r"
                ].to_numpy()
                check(
                    f"RNA_model_{gsm}_{population}_{model}",
                    np.allclose(values, expected, equal_nan=True, atol=1e-10),
                    "All 170 partial correlations independently reproduced by QR projection",
                )
    h, l = np.vstack(high), np.vstack(low)
    ho, lo = ((h > 0) & (l == 0)).sum(0), ((l > 0) & (h == 0)).sum(0)
    p = np.array(
        [
            stats.binomtest(int(a), int(a + b)).pvalue if a + b else 1
            for a, b in zip(ho, lo)
        ]
    )
    fc = np.divide(
        h.sum(0), l.sum(0), out=np.full(h.shape[1], np.nan), where=l.sum(0) > 0
    )
    check(
        "all_peak_counts_FC_p_q",
        np.array_equal(h.sum(0), peaks.high)
        and np.array_equal(l.sum(0), peaks.low)
        and np.array_equal(ho, peaks.high_only)
        and np.array_equal(lo, peaks.low_only)
        and np.allclose(fc, peaks.FC, equal_nan=True)
        and np.allclose(p, peaks.p, rtol=1e-12)
        and np.allclose(bh(p), peaks.q565, rtol=1e-12),
        "All 565 count contrasts, fold changes, exact p and BH565 reproduced",
    )
    delta = h - l
    groups = {g: peak_index.get_indexer(d.peak) for g, d in member.groupby("gene")}
    z = delta.sum(0) / np.sqrt(np.maximum((delta**2).sum(0), 1))
    for mode, sign in [("opening", 1), ("closing", -1)]:
        d = regions[regions["mode"].eq(mode)].set_index("gene").loc[scope]
        observed = [max(0, (sign * z[groups[g]]).max()) for g in scope]
        best = [
            peaks.peak.iloc[groups[g][np.argmax(sign * z[groups[g]])]] for g in scope
        ]
        check(
            "regional_statistics_" + mode,
            np.allclose(observed, d.statistic, atol=1e-5)
            and all(
                np.isclose(
                    sign * z[peak_index.get_loc(a)], sign * z[peak_index.get_loc(b)]
                )
                for a, b in zip(best, d.top_peak)
            )
            and np.allclose((d.exceedances + 1) / (d.permutations + 1), d.p)
            and np.allclose(bh(d.p), d.q25),
            "All 25 maxima, peak identities, permutation p formula and BH25 checked",
        )
    if args.permutations:
        permutation = pd.DataFrame(
            grid.permutation(delta, groups, scope, 2000000, 20261134)
        )
        permutation.to_csv(
            output / "regional_permutation_rerun.tsv", sep="\t", index=False
        )
        joined = permutation.merge(
            regions, on=["mode", "gene"], suffixes=("_new", "_saved")
        )
        check(
            "regional_permutation_rerun",
            np.array_equal(joined.exceedances_new, joined.exceedances_saved),
            "2,000,000 swaps; seed 20261134; all 50 exceedance counts exactly reproduced",
        )
    rh, rl = np.vstack(rna_high), np.vstack(rna_low)
    difference = np.log1p(rh) - np.log1p(rl)
    rp = stats.ttest_1samp(difference, 0, axis=0).pvalue
    rfc = np.divide(
        rh.mean(0), rl.mean(0), out=np.full(len(genes), np.nan), where=rl.mean(0) > 0
    )
    saved = read(SOURCE / "RNA_high_low.tsv").set_index("gene").loc[genes]
    check(
        "RNA_contrasts",
        np.allclose(rp, saved.p, equal_nan=True)
        and np.allclose(rfc, saved.FC, equal_nan=True)
        and np.allclose(
            bh(np.nan_to_num(rp, nan=1))[np.isfinite(rp)],
            saved.q_local_RNAs[np.isfinite(rp)],
        ),
        "All 170 RNA fold changes, paired p and BH170 reproduced",
    )
    combined = read(SOURCE / "links_combined.tsv")
    for key, d in source_links.groupby(["population", "model", "peak", "gene"]):
        row = combined.set_index(["population", "model", "peak", "gene"]).loc[key]
        ok = d.r.notna() & d.RNA_detected.ge(10) & d.peak_detected.ge(10)
        assert row.estimable_sources == ok.sum()
        if ok.all():
            w = d.df - 1
            fisher = np.average(
                np.arctanh(np.clip(d.r, -0.999999, 0.999999)), weights=w
            )
            assert np.isclose(np.tanh(fisher), row.r)
            assert np.isclose(
                2 * stats.norm.sf(abs(fisher) * np.sqrt(w.sum())),
                row.p,
                rtol=1e-10,
                atol=1e-200,
            )
        else:
            assert pd.isna(row.p)
    for _, d in combined.groupby(["population", "model"]):
        assert np.allclose(
            bh(d.p.fillna(1))[d.p.notna()], d.loc[d.p.notna(), "q_all_local_links"]
        )
    check(
        "combined_links",
        True,
        "680 combined estimates: Fisher z, detection eligibility and four BH170 families checked",
    )
    environment = {
        package: importlib.metadata.version(package)
        for package in [
            "numpy",
            "pandas",
            "scipy",
            "h5py",
            "statsmodels",
            "scikit-learn",
            "matplotlib",
            "threadpoolctl",
        ]
    }
    environment["python"] = sys.version
    report = {
        "passed": True,
        "checks": CHECKS,
        "permutation_rerun": args.permutations,
        "sources": 2,
        "pairs": len(h),
        "environment": environment,
    }
    (output / "audit_report.json").write_text(json.dumps(report, indent=2) + "\n")
    files = [p for p in SOURCE.glob("*.tsv*")] + [SOURCE / "design.json"]
    manifest = [
        {
            "file": p.relative_to(ROOT).as_posix(),
            "sha256": hashlib.sha256(p.read_bytes()).hexdigest(),
        }
        for p in sorted(files)
    ]
    pd.DataFrame(manifest).to_csv(
        output / "source_table_sha256.tsv", sep="\t", index=False
    )
    print("Audit complete", len(CHECKS), "checks", flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--work", type=Path, required=True)
    parser.add_argument("--h5", type=Path, required=True)
    parser.add_argument("--gtf", type=Path, required=True)
    parser.add_argument("--permutations", action="store_true")
    parser.add_argument("--correct-tss-distances", action="store_true")
    main(parser.parse_args())
