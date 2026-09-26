"""Stream GSE109828 counts, retain cell QC and externally overlapping regions.

Coordinates are harmonized before any COQ8A effects are accessed. The external
cell filter follows Pliner et al.: >=1000 accessible sites and accessibility at
one of seven myogenic promoters. Promoters here span +/-2 kb of annotated
transcript TSSs; +/-1 kb is retained for a filter sensitivity analysis.
"""

import argparse
import gzip
import re
from pathlib import Path

import numpy as np
import pandas as pd
from pyliftover import LiftOver
from scipy import sparse

MARKERS = {"MYOG", "MYOD1", "DMD", "TNNT1", "MYH1", "MYH3", "TPM2"}


def lift_interval(lo, chrom, start, end):
    hits = [
        lo.convert_coordinate(chrom, int(x))
        for x in [start, (start + end) // 2, end - 1]
    ]
    if any(len(x) != 1 for x in hits):
        return None
    a, m, b = [x[0] for x in hits]
    if not (a[0] == m[0] == b[0] and a[2] == m[2] == b[2]):
        return None
    left, right = sorted([int(a[1]), int(b[1])])
    if abs((right - left + 1) / (end - start) - 1) > 0.05:
        return None
    return a[0], left, right + 1


def reference_regions(a):
    lo = LiftOver(str(a.raw / "hg38ToHg19.over.chain.gz"))
    candidates = pd.read_csv(a.tables / "candidate_peak_gene.tsv.gz", sep="\t")
    records = []
    for peak in sorted(candidates.loc[candidates.n_libraries == 4, "peak"].unique()):
        chrom, start, end = re.split("[:-]", peak)
        hit = lift_interval(lo, chrom, int(start), int(end))
        records.append(
            dict(
                peak=peak,
                chrom19=hit[0] if hit else "",
                start19=hit[1] if hit else -1,
                end19=hit[2] if hit else -1,
                mapping_pass=hit is not None,
            )
        )
    lifted = pd.DataFrame(records)
    lifted.to_csv(a.out / "candidate_hg19_mapping.tsv", sep="\t", index=False)
    promoters = []
    with gzip.open(a.gtf, "rt", encoding="utf-8") as stream:
        for line in stream:
            if line.startswith("#") or not any(
                f'gene_name "{g}"' in line for g in MARKERS
            ):
                continue
            fields = line.rstrip().split("\t")
            if fields[2] != "transcript":
                continue
            name = re.search(r'gene_name "([^"]+)"', fields[8]).group(1)
            tss = int(fields[3]) - 1 if fields[6] == "+" else int(fields[4]) - 1
            for width in [1000, 2000]:
                hit = lift_interval(lo, fields[0], max(0, tss - width), tss + width + 1)
                if hit:
                    promoters.append(
                        dict(
                            gene=name,
                            width=width,
                            chrom19=hit[0],
                            start19=hit[1],
                            end19=hit[2],
                        )
                    )
    promoters = pd.DataFrame(promoters).drop_duplicates()
    assert set(promoters.gene) == MARKERS
    promoters.to_csv(
        a.out / "myogenic_marker_promoters_hg19.tsv", sep="\t", index=False
    )
    return lifted[lifted.mapping_pass], promoters


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for arg in ["raw", "tables", "gtf", "out"]:
        p.add_argument("--" + arg, type=Path, required=True)
    p.add_argument("--experiment", choices=["1", "2"], required=True)
    a = p.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)
    if (a.out / "candidate_hg19_mapping.tsv").exists() and (
        a.out / "myogenic_marker_promoters_hg19.tsv"
    ).exists():
        lifted = pd.read_csv(a.out / "candidate_hg19_mapping.tsv", sep="\t")
        lifted = lifted[lifted.mapping_pass]
        promoters = pd.read_csv(a.out / "myogenic_marker_promoters_hg19.tsv", sep="\t")
    else:
        lifted, promoters = reference_regions(a)
    targets = {
        c: (d.start19.to_numpy(), d.end19.to_numpy(), d.peak.to_numpy())
        for c, d in lifted.groupby("chrom19")
    }
    prom = {
        c: (d.start19.to_numpy(), d.end19.to_numpy(), d.width.to_numpy())
        for c, d in promoters.groupby("chrom19")
    }
    exp = "HSMM" + a.experiment
    meta = pd.read_csv(next(a.raw.glob(f"*{exp}_indextable.txt.gz")), sep="\t")
    assert not meta.barcode.duplicated().any()
    meta["hour"] = meta.timepoint.str.extract(r"(\d+)").astype(int)
    barcode_index = {b: i for i, b in enumerate(meta.barcode)}
    qc_sites = np.zeros(len(meta), dtype=np.int64)
    qc_reads = np.zeros(len(meta), dtype=np.int64)
    marker1 = np.zeros(len(meta), dtype=bool)
    marker2 = marker1.copy()
    cache = {}
    selected = []
    overlaps = []
    rows = []
    cols = []
    seen_rows = 0
    unknown_rows = 0
    counts = next(a.raw.glob(f"*{exp}_counts.txt.gz"))
    for k, chunk in enumerate(
        pd.read_csv(
            counts,
            sep="\t",
            header=None,
            names=["region", "barcode", "count"],
            dtype={"count": "int32"},
            chunksize=1_000_000,
        )
    ):
        for region in chunk.region.unique():
            if region in cache:
                continue
            chrom, start, end = region.split("_")
            start, end = int(start), int(end)
            links = []
            m1 = m2 = False
            if chrom in targets:
                starts, ends, names = targets[chrom]
                mask = (starts < end) & (ends > start)
                for s, e, name in zip(starts[mask], ends[mask], names[mask]):
                    ov = min(end, e) - max(start, s)
                    reciprocal = min(ov / (end - start), ov / (e - s))
                    links.append(
                        dict(
                            external_peak=region,
                            peak=name,
                            overlap_bp=ov,
                            reciprocal_overlap=reciprocal,
                        )
                    )
            if chrom in prom:
                starts, ends, widths = prom[chrom]
                hit = widths[(starts < end) & (ends > start)]
                m1 = bool((hit == 1000).any())
                m2 = bool((hit == 2000).any())
            col = -1
            if links:
                col = len(selected)
                selected.append(region)
                overlaps.extend(links)
            cache[region] = (col, m1, m2)
        cell = chunk.barcode.map(barcode_index)
        valid = cell.notna().to_numpy()
        unknown_rows += int((~valid).sum())
        ci = cell[valid].to_numpy(dtype=np.int32)
        qc_sites += np.bincount(ci, minlength=len(meta))
        qc_reads += np.bincount(
            ci, weights=chunk.loc[valid, "count"].to_numpy(), minlength=len(meta)
        ).astype(np.int64)
        codes = np.array(chunk.loc[valid, "region"].map(cache).tolist(), dtype=np.int32)
        marker1[ci[codes[:, 1] > 0]] = True
        marker2[ci[codes[:, 2] > 0]] = True
        keep = codes[:, 0] >= 0
        rows.append(ci[keep])
        cols.append(codes[keep, 0])
        seen_rows += len(chunk)
        if k % 10 == 0:
            print(
                exp,
                "rows",
                seen_rows,
                "unique regions",
                len(cache),
                "retained regions",
                len(selected),
                flush=True,
            )
    matrix = sparse.coo_matrix(
        (
            np.ones(sum(map(len, rows)), dtype=np.int8),
            (np.concatenate(rows), np.concatenate(cols)),
        ),
        shape=(len(meta), len(selected)),
    ).tocsr()
    duplicate_selected = int((matrix.data > 1).sum())
    matrix.data[:] = 1
    meta["open_sites"] = qc_sites
    meta["peak_reads"] = qc_reads
    meta["myogenic_promoter_1kb"] = marker1
    meta["myogenic_promoter_2kb"] = marker2
    meta["pass_1kb"] = (qc_sites >= 1000) & marker1
    meta["pass_2kb"] = (qc_sites >= 1000) & marker2
    meta.to_csv(a.out / f"{exp}_cells.tsv.gz", sep="\t", index=False)
    sparse.save_npz(a.out / f"{exp}_selected_binary.npz", matrix)
    pd.DataFrame({"external_peak": selected}).to_csv(
        a.out / f"{exp}_selected_regions.tsv", sep="\t", index=False
    )
    pd.DataFrame(overlaps).to_csv(
        a.out / f"{exp}_overlaps.tsv.gz", sep="\t", index=False
    )
    pd.DataFrame(
        [
            dict(
                experiment=exp,
                rows=seen_rows,
                unknown_barcode_rows=unknown_rows,
                total_external_regions=len(cache),
                selected_regions=len(selected),
                duplicate_selected_cell_peak=duplicate_selected,
            )
        ]
    ).to_csv(a.out / f"{exp}_ingest_audit.tsv", sep="\t", index=False)
    print(meta.groupby("hour")[["pass_1kb", "pass_2kb"]].sum().to_string(), flush=True)


if __name__ == "__main__":
    main()
