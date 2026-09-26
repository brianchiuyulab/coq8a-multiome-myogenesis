"""Extract RNA/ATAC depth and COQ8A UMIs from four public 10x multiome H5 files.

Example: python scripts/01_extract_nuclei.py --h5-root PATH --out results/tables/nuclei.tsv.gz
"""

import argparse
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

SAMPLES = {
    "GSM6339597": ("line1", "stem"),
    "GSM6339599": ("line1", "differentiated"),
    "GSM6339601": ("line2", "stem"),
    "GSM6339603": ("line2", "differentiated"),
}


def extract(path: Path, metrics_path: Path, gsm: str) -> pd.DataFrame:
    with h5py.File(path) as h:
        m = h["matrix"]
        f = m["features"]
        feature_type = f["feature_type"][:]
        names = f["name"][:]
        coq = np.flatnonzero((feature_type == b"Gene Expression") & (names == b"COQ8A"))
        assert len(coq) == 1, f"{gsm}: COQ8A feature is ambiguous"
        ng = int((feature_type == b"Gene Expression").sum())
        assert np.all(feature_type[:ng] == b"Gene Expression")
        barcodes = np.char.decode(m["barcodes"][:], "utf-8")
        ptr = m["indptr"][:]
        rna = np.zeros(len(barcodes), dtype=np.int32)
        atac = np.zeros(len(barcodes), dtype=np.int32)
        coq_umi = np.zeros(len(barcodes), dtype=np.int32)
        for left in range(0, len(barcodes), 250):
            right = min(left + 250, len(barcodes))
            start, end = int(ptr[left]), int(ptr[right])
            idx, val = m["indices"][start:end], m["data"][start:end]
            rel = ptr[left : right + 1] - start
            for j in range(right - left):
                lo, hi = int(rel[j]), int(rel[j + 1])
                item, count = idx[lo:hi], val[lo:hi]
                rna[left + j] = int(count[item < ng].sum())
                atac[left + j] = int((item >= ng).sum())
                hit = np.flatnonzero(item == coq[0])
                if len(hit):
                    coq_umi[left + j] = int(count[hit[0]])
    metric = pd.read_csv(
        metrics_path,
        usecols=[
            "barcode",
            "is_cell",
            "excluded_reason",
            "gex_umis_count",
            "atac_fragments",
            "atac_peak_region_fragments",
        ],
    )
    source, stage = SAMPLES[gsm]
    out = pd.DataFrame(
        {
            "gsm": gsm,
            "barcode": barcodes,
            "source": source,
            "stage": stage,
            "total_rna_umi": rna,
            "total_open_peaks": atac,
            "COQ8A_umi": coq_umi,
        }
    )
    out = out.merge(metric, on="barcode", how="left", validate="one_to_one")
    assert out.is_cell.eq(1).all(), f"{gsm}: uncalled barcode in filtered H5"
    assert out.excluded_reason.fillna(0).eq(0).all(), f"{gsm}: excluded barcode in filtered H5"
    assert out.gex_umis_count.notna().all()
    return out


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--h5-root", type=Path, required=True)
    p.add_argument("--out", type=Path, required=True)
    a = p.parse_args()
    rows = []
    for gsm in SAMPLES:
        matrix = next(a.h5_root.glob(f"{gsm}_*filtered_feature_bc_matrix.h5"))
        metrics = next(a.h5_root.glob(f"{gsm}_*per_barcode_metrics.csv.gz"))
        d = extract(matrix, metrics, gsm)
        print(gsm, len(d), "filtered barcodes", flush=True)
        rows.append(d)
    a.out.parent.mkdir(parents=True, exist_ok=True)
    pd.concat(rows, ignore_index=True).to_csv(
        a.out, sep="\t", index=False, compression={"method": "gzip", "mtime": 0}
    )


if __name__ == "__main__":
    main()
