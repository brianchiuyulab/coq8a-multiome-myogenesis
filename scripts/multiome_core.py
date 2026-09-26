"""Shared 10x H5 reader and paired descriptive statistics."""

from pathlib import Path

import h5py
import numpy as np
from scipy import sparse, stats


def h5_for(root: Path, gsm: str) -> Path:
    return next(root.glob(f"{gsm}_*filtered_feature_bc_matrix.h5"))


def read_barcodes(h5: h5py.File, barcodes):
    matrix = h5["matrix"]
    features = matrix["features"]
    names = np.char.decode(features["name"][:], "utf-8")
    types = features["feature_type"][:]
    ng = int((types == b"Gene Expression").sum())
    assert np.all(types[:ng] == b"Gene Expression") and np.all(types[ng:] == b"Peaks")
    all_barcodes = np.char.decode(matrix["barcodes"][:], "utf-8")
    index = {b: i for i, b in enumerate(all_barcodes)}
    r_rows, r_cols, r_data, p_rows, p_cols = [], [], [], [], []
    ptr = matrix["indptr"][:]
    for row, barcode in enumerate(barcodes):
        col = index[barcode]
        left, right = int(ptr[col]), int(ptr[col + 1])
        ids, values = matrix["indices"][left:right], matrix["data"][left:right]
        is_rna = ids < ng
        if is_rna.any():
            r_rows.extend([row] * int(is_rna.sum()))
            r_cols.extend(ids[is_rna].tolist())
            r_data.extend(values[is_rna].tolist())
        if (~is_rna).any():
            p_rows.extend([row] * int((~is_rna).sum()))
            p_cols.extend((ids[~is_rna] - ng).tolist())
    rna = sparse.csr_matrix(
        (np.asarray(r_data, dtype=np.float32), (r_rows, r_cols)), shape=(len(barcodes), ng)
    )
    atac = sparse.csr_matrix(
        (np.ones(len(p_cols), dtype=np.int8), (p_rows, p_cols)),
        shape=(len(barcodes), len(names) - ng),
    )
    atac.sum_duplicates()
    atac.data[:] = 1
    return rna, atac, names[:ng], names[ng:]


def paired_stats(difference):
    d = np.asarray(difference, dtype=float)
    n = len(d)
    if n < 2:
        return {
            "n_pairs": n,
            "difference": float(np.mean(d)) if n else np.nan,
            "ci_low": np.nan,
            "ci_high": np.nan,
            "p_pair": np.nan,
        }
    mean = float(d.mean())
    se = float(stats.sem(d))
    t = float(stats.t.ppf(0.975, n - 1))
    p = float(stats.ttest_1samp(d, 0).pvalue) if se > 0 else 1.0
    return {
        "n_pairs": n,
        "difference": mean,
        "ci_low": mean - t * se,
        "ci_high": mean + t * se,
        "p_pair": p,
    }
