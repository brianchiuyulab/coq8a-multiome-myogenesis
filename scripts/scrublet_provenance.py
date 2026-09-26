"""Independent RNA-based multiplet sensitivity screen for the multiome nuclei.

Scrublet is run separately on each 10x library. This complements, but does not
replace, ATAC-specific doublet detection in a future publication workflow.
"""

from pathlib import Path
import argparse
import h5py
import numpy as np
import pandas as pd
from scipy import sparse
import scrublet as scr


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--h5-root", type=Path, required=True)
    p.add_argument("--barcodes", type=Path, required=True)
    p.add_argument("--out", type=Path, required=True)
    a = p.parse_args()
    bc = pd.read_csv(a.barcodes, sep="\t")
    out = []
    for gsm in sorted(bc.gsm.unique()):
        wanted = bc.loc[bc.gsm == gsm, "barcode"].tolist()
        path = next(a.h5_root.glob(gsm + "*filtered_feature_bc_matrix.h5"))
        with h5py.File(path) as h:
            m = h["matrix"]
            f = m["features"]
            ng = int(np.sum(f["feature_type"][:] == b"Gene Expression"))
            saved = np.char.decode(m["barcodes"][:], "utf8")
            lookup = {x: i for i, x in enumerate(saved)}
            ptr = m["indptr"][:]
            row = []
            col = []
            val = []
            for i, b in enumerate(wanted):
                c = lookup[b]
                start = int(ptr[c])
                stop = int(ptr[c + 1])
                ids = m["indices"][start:stop]
                use = ids < ng
                row.extend([i] * int(use.sum()))
                col.extend(ids[use])
                val.extend(m["data"][start:stop][use])
        X = sparse.csr_matrix(
            (np.asarray(val, dtype=np.float32), (row, col)), shape=(len(wanted), ng)
        )
        print(gsm, "Scrublet input", X.shape, X.nnz, flush=True)
        detector = scr.Scrublet(X, expected_doublet_rate=0.06, random_state=20260926)
        score, pred = detector.scrub_doublets(
            min_counts=2,
            min_cells=3,
            min_gene_variability_pctl=85,
            n_prin_comps=30,
            use_approx_neighbors=False,
        )
        q = pd.DataFrame(
            {
                "gsm": gsm,
                "barcode": wanted,
                "scrublet_score": score,
                "scrublet_predicted_doublet": pred,
            }
        )
        out.append(q)
        print(
            gsm,
            "doublets",
            int(np.sum(pred)),
            "/",
            len(pred),
            "threshold",
            detector.threshold_,
            flush=True,
        )
    a.out.parent.mkdir(parents=True, exist_ok=True)
    pd.concat(out, ignore_index=True).to_csv(
        a.out, sep="\t", index=False, compression={"method": "gzip", "mtime": 0}
    )


if __name__ == "__main__":
    main()
