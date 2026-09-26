"""Map all shared human candidate peaks to independent C2C12 peak calls."""

import argparse
from bisect import bisect_left
import gzip
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd


def map_positions(chain_path, peaks):
    """Stream UCSC chains and retain only the requested source coordinates."""
    requested = {}
    for peak in peaks:
        chrom, interval = peak.split(":")
        start, end = map(int, interval.split("-"))
        requested.setdefault(chrom, set()).update([start, (start + end) // 2, end - 1])
    requested = {chrom: sorted(coords) for chrom, coords in requested.items()}
    mapped = {(chrom, pos): [] for chrom, coords in requested.items() for pos in coords}
    active = False
    with gzip.open(chain_path, "rt") as handle:
        for line in handle:
            if line.startswith("chain "):
                fields = line.split()
                chrom = fields[2]
                assert fields[4] == "+", "Unexpected reverse source chain"
                tp, te = int(fields[5]), int(fields[6])
                qchrom, qsize, strand = fields[7], int(fields[8]), fields[9]
                qp, chain_id = int(fields[10]), int(fields[12])
                coords = requested.get(chrom, [])
                left, right = bisect_left(coords, tp), bisect_left(coords, te)
                active = left < right
            elif active and line.strip():
                fields = list(map(int, line.split()))
                size = fields[0]
                left, right = bisect_left(coords, tp), bisect_left(coords, tp + size)
                for pos in coords[left:right]:
                    qpos = qp + pos - tp
                    if strand == "-":
                        qpos = qsize - qpos - 1
                    mapped[chrom, pos].append((qchrom, qpos, strand, chain_id))
                if len(fields) == 3:
                    tp += size + fields[1]
                    qp += size + fields[2]
                else:
                    active = False
    return mapped


def main(a):
    a.out.mkdir(parents=True, exist_ok=True)
    effects = pd.read_csv(a.root / "results/unstratified/ATAC_all_5097.tsv", sep="\t")
    assert len(effects) == 5097
    position_map = map_positions(a.reference / "hg38ToMm10.over.chain.gz", effects.peak)
    print("Coordinate mapping complete", flush=True)
    references = {}
    for stage in ["GM", "DM60h"]:
        d = pd.read_csv(
            a.reference / f"GSE224489_C2C12_{stage}_ATAC_peaks.broadPeak.gz",
            sep="\t",
            header=None,
        )
        d = d.iloc[:, :3]
        d.columns = ["chrom", "start", "end"]
        references[stage] = dict(tuple(d.groupby("chrom")))
    rows = []
    for peak in effects.peak:
        chrom, interval = peak.split(":")
        start, end = map(int, interval.split("-"))
        mapped = [
            position_map[chrom, pos] for pos in [start, (start + end) // 2, end - 1]
        ]
        r = dict(peak=peak, mapping_pass=False, reason="nonunique_or_unmapped")
        if all(m is not None and len(m) == 1 for m in mapped):
            m = [x[0] for x in mapped]
            if len({(v[0], v[2], v[3]) for v in m}) == 1:
                s, e = min(m[0][1], m[2][1]), max(m[0][1], m[2][1]) + 1
                ratio = (e - s) / (end - start)
                inside = s <= m[1][1] < e
                if inside and 0.5 <= ratio <= 2:
                    r.update(
                        mapping_pass=True,
                        reason="pass",
                        chrom_mm10=m[0][0],
                        start_mm10=s,
                        end_mm10=e,
                        strand=m[0][2],
                        length_ratio=ratio,
                    )
                    for stage, refs in references.items():
                        g = refs.get(m[0][0])
                        if g is None:
                            r[stage + "_overlap"] = False
                            r[stage + "_covered_fraction"] = 0
                            continue
                        g = g[(g.start < e) & (g.end > s)]
                        intervals = sorted(
                            (max(s, row.start), min(e, row.end))
                            for row in g.itertuples()
                        )
                        union = []
                        for left, right in intervals:
                            if union and left <= union[-1][1]:
                                union[-1][1] = max(right, union[-1][1])
                            else:
                                union.append([left, right])
                        r[stage + "_overlap"] = len(g) > 0
                        r[stage + "_covered_fraction"] = sum(
                            right - left for left, right in union
                        ) / (e - s)
                else:
                    r["reason"] = "length_or_center_inconsistency"
            else:
                r["reason"] = "inconsistent_chain"
        rows.append(r)
    out = pd.DataFrame(rows)
    out.to_csv(
        a.out / "human_to_mouse_ATAC_annotation.tsv", sep="\t", index=False, na_rep="NA"
    )
    paths = [a.reference / "hg38ToMm10.over.chain.gz"] + list(
        a.reference.glob("*.broadPeak.gz")
    )
    manifest = dict(
        assembly="mm10, verified in GSM7025063 metadata",
        hg38_inputs=5097,
        unique_consistent_mappings=int(out.mapping_pass.sum()),
        conditions="GM and DM60h; author merged-replicate broad peak calls",
        mapping_rule="Unique start, midpoint and end-1 on same chain/strand; span ratio 0.5 to 2",
        interpretation="Peak-call overlap, not a quantitative accessibility fold change or a paired-nucleus assay",
        sources=[
            dict(file=p.name, sha256=hashlib.sha256(p.read_bytes()).hexdigest())
            for p in paths
        ],
    )
    (a.out / "mouse_ATAC_reference_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n"
    )
    print(out.groupby("reason").size().to_string())
    focal = [
        "chr11:17649919-17650798",
        "chr3:8733438-8733968",
        "chrX:31220862-31221781",
        "chr11:19201752-19202603",
    ]
    print(out[out.peak.isin(focal)].to_string(index=False))


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    p.add_argument("--reference", type=Path, required=True)
    p.add_argument("--out", type=Path, required=True)
    main(p.parse_args())
