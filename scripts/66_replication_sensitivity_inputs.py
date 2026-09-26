"""Prepare a declared sensitivity grid and genome-wide donor pseudobulks."""

import argparse
import json
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import sparse, io
from sklearn.neighbors import NearestNeighbors

CELLTYPES = ["Satellite Cells", "Fast", "Slow", "Intermediate"]
RULES = [
    "2plus_vs_1",
    "3plus_vs_1",
    "4plus_vs_1",
    "3plus_vs_1to2",
    "positive_q25",
    "positive_q33",
    "detected_vs_zero",
]


def groups(d, rule):
    x = d.COQ8A_raw_UMI.to_numpy()
    if rule in ["2plus_vs_1", "3plus_vs_1", "4plus_vs_1"]:
        return np.flatnonzero(x >= int(rule[0])), np.flatnonzero(x == 1)
    if rule == "3plus_vs_1to2":
        return np.flatnonzero(x >= 3), np.flatnonzero((x >= 1) & (x <= 2))
    if rule == "detected_vs_zero":
        return np.flatnonzero(x >= 1), np.flatnonzero(x == 0)
    fraction = 0.25 if rule == "positive_q25" else 1 / 3
    v = x / d.raw_RNA_library_size.to_numpy()
    pos = v[x > 0]
    if len(pos) < 2:
        return np.array([], int), np.array([], int)
    lo, hi = np.quantile(pos, [fraction, 1 - fraction])
    if hi <= lo:
        return np.array([], int), np.array([], int)
    return np.flatnonzero((x > 0) & (v >= hi)), np.flatnonzero((x > 0) & (v <= lo))


def match(d, hi, lo, caliper):
    if not len(hi) or not len(lo):
        return np.array([], int), np.array([], int)
    if caliper is None:
        return hi, lo
    x = np.log1p(d[["raw_RNA_library_size", "nFeature_ATAC"]].to_numpy())
    nn = NearestNeighbors(n_neighbors=min(100, len(lo)), algorithm="kd_tree").fit(x[lo])
    dist, indices = nn.kneighbors(x[hi])
    used = set()
    h = []
    l = []
    for i in np.argsort(-dist[:, 0]):
        for j in indices[i]:
            if int(j) not in used and np.max(np.abs(x[hi[i]] - x[lo[j]])) <= caliper:
                used.add(int(j))
                h.append(hi[i])
                l.append(lo[j])
                break
    return np.asarray(h, int), np.asarray(l, int)


def main(root, work, out):
    out.mkdir(parents=True, exist_ok=True)
    meta = pd.read_csv(work / "nucleus_inventory.tsv", sep="\t")
    barcodes = (work / "all_ATAC_barcodes.txt").read_text().splitlines()
    meta = meta.set_index("barcode").loc[barcodes].reset_index()
    meta["cell_index"] = np.arange(len(meta))
    meta["donor"] = meta.Sample.str.split("_").str[0]
    meta["celltype"] = meta["refined_annotations_wknn_0.8"]
    caps = meta.groupby(["Sample", "celltype"])[
        ["raw_RNA_library_size", "nFeature_ATAC"]
    ].transform(lambda z: z.quantile(0.95))
    meta["trim95"] = (meta.raw_RNA_library_size <= caps.raw_RNA_library_size) & (
        meta.nFeature_ATAC <= caps.nFeature_ATAC
    )
    configs = []
    cols = []
    members = []
    full_pairs = []
    cohort_rows = []
    groupid = 0
    for ctype in CELLTYPES:
        for time in ["Pre", "Post"]:
            for qc in ["author", "trim95"]:
                base = meta[meta.celltype.eq(ctype) & meta.Time.eq(time)]
                if qc == "trim95":
                    base = base[base.trim95]
                for rule in RULES:
                    for caliper in [None, 0.1, 0.2, 0.3]:
                        cfg = f"C{len(configs):03d}"
                        config = dict(
                            config=cfg,
                            celltype=ctype,
                            time=time,
                            qc=qc,
                            rule=rule,
                            caliper="none" if caliper is None else str(caliper),
                        )
                        nsource = 0
                        nh = nl = 0
                        for sample, d in base.groupby("Sample", sort=True):
                            d = d.reset_index(drop=True)
                            hi, lo = groups(d, rule)
                            raw_h, raw_l = len(hi), len(lo)
                            h, l = match(d, hi, lo, caliper)
                            if not len(h) or not len(l):
                                cohort_rows.append(
                                    config
                                    | dict(
                                        Sample=sample,
                                        donor=d.donor.iloc[0],
                                        raw_high=raw_h,
                                        raw_low=raw_l,
                                        n_high=0,
                                        n_low=0,
                                    )
                                )
                                continue
                            nsource += 1
                            nh += len(h)
                            nl += len(l)
                            means = np.log1p(
                                d[["raw_RNA_library_size", "nFeature_ATAC"]].to_numpy()
                            )
                            imbalance = means[h].mean(0) - means[l].mean(0)
                            cohort_rows.append(
                                config
                                | dict(
                                    Sample=sample,
                                    donor=d.donor.iloc[0],
                                    raw_high=raw_h,
                                    raw_low=raw_l,
                                    n_high=len(h),
                                    n_low=len(l),
                                    log_RNA_imbalance=imbalance[0],
                                    log_ATAC_imbalance=imbalance[1],
                                )
                            )
                            for label, indices in [("low", l), ("high", h)]:
                                g = d.iloc[indices]
                                cols.append(
                                    config
                                    | dict(
                                        column=groupid,
                                        Sample=sample,
                                        donor=d.donor.iloc[0],
                                        group=label,
                                        n_nuclei=len(g),
                                        raw_ATAC_size=g.raw_ATAC_library_size.sum(),
                                        raw_RNA_size=g.raw_RNA_library_size.sum(),
                                        mean_log_RNA=np.log1p(
                                            g.raw_RNA_library_size
                                        ).mean(),
                                        mean_log_ATAC=np.log1p(g.nFeature_ATAC).mean(),
                                    )
                                )
                                members.extend((int(i), groupid) for i in g.cell_index)
                                groupid += 1
                            if caliper is not None:
                                full_pairs.extend(
                                    dict(
                                        config=cfg,
                                        donor=d.donor.iloc[0],
                                        high_index=int(d.iloc[i].cell_index),
                                        low_index=int(d.iloc[j].cell_index),
                                    )
                                    for i, j in zip(h, l)
                                )
                        configs.append(
                            config | dict(n_donors=nsource, n_high=nh, n_low=nl)
                        )
    conf = pd.DataFrame(configs)
    columns = pd.DataFrame(cols)
    conf.to_csv(out / "configurations.tsv", sep="\t", index=False)
    columns.to_csv(out / "pseudobulk_columns.tsv", sep="\t", index=False)
    pd.DataFrame(cohort_rows).to_csv(
        out / "sample_balance.tsv", sep="\t", index=False, na_rep="NA"
    )
    pd.DataFrame(full_pairs).to_csv(
        work / "sensitivity_pairs.tsv.gz", sep="\t", index=False
    )
    arr = np.asarray(members, dtype=np.int32)
    membership = sparse.csc_matrix(
        (np.ones(len(arr), np.int32), (arr[:, 0], arr[:, 1])),
        shape=(len(meta), groupid),
    )
    sparse.save_npz(work / "sensitivity_membership.npz", membership)
    shape = tuple(map(int, (work / "all_ATAC_shape.txt").read_text().splitlines()))
    xp = np.memmap(work / "all_ATAC_x.bin", dtype="<i4", mode="r")
    ix = np.memmap(work / "all_ATAC_i.bin", dtype="<i4", mode="r")
    ptr = np.memmap(work / "all_ATAC_p.bin", dtype="<i4", mode="r")
    atac = sparse.csc_matrix((xp, ix, ptr), shape=shape, copy=False)
    assert np.array_equal(
        np.asarray(atac.sum(axis=0)).ravel(), meta.raw_ATAC_library_size
    )
    output = np.memmap(
        work / "sensitivity_ATAC_pseudobulk.bin",
        dtype="<i4",
        mode="w+",
        shape=(shape[0], groupid),
        order="F",
    )
    for start in range(0, groupid, 96):
        stop = min(start + 96, groupid)
        block = (atac @ membership[:, start:stop]).toarray()
        assert block.min() >= 0
        output[:, start:stop] = block
        output.flush()
        print("ATAC pseudobulk", stop, "/", groupid, flush=True)
    candidate = io.mmread(work / "candidate_ATAC_counts.mtx").tocsc()
    candidate_names = (work / "candidate_ATAC_peaks.txt").read_text().splitlines()
    nuc_barcodes = (work / "barcodes.txt").read_text().splitlines()
    reorder = pd.Index(nuc_barcodes).get_indexer(barcodes)
    candidate = candidate[:, reorder]
    opened = candidate.copy()
    opened.data[:] = 1
    np.save(
        work / "sensitivity_open_counts.npy",
        (opened @ membership).toarray().astype(np.int32),
    )
    rna = io.mmread(work / "candidate_RNA_counts.mtx").tocsc()[:, reorder]
    io.mmwrite(work / "sensitivity_RNA_pseudobulk.mtx", rna @ membership)
    (out / "matrix_shapes.json").write_text(
        json.dumps(
            dict(
                ATAC=list(output.shape),
                RNA=[rna.shape[0], groupid],
                candidate_peaks=len(candidate_names),
            ),
            indent=2,
        )
        + "\n"
    )
    manifest = dict(
        celltypes=CELLTYPES,
        times=["Pre", "Post"],
        QC=["author retained", "within sample/celltype 95th-percentile depth caps"],
        rules=RULES,
        calipers=["none", 0.1, 0.2, 0.3],
        n_configurations=len(conf),
        TSS_sensitivity="not available in author object; not claimed",
        feature_basis="all author ATAC peaks for normalization; reciprocal-50%-overlap candidate scopes for reporting",
        primary_anchor="author QC; 3plus_vs_1; caliper .10; Pre, each celltype separately",
        statistics="donor-blocked edgeR QL for >=2 usable donors; no donor dropping based on outcomes",
    )
    (out / "grid_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print("Complete:", len(conf), "settings,", groupid, "pseudobulks", flush=True)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    p.add_argument("--work", type=Path, required=True)
    p.add_argument("--out", type=Path, required=True)
    a = p.parse_args()
    main(a.root, a.work, a.out)
