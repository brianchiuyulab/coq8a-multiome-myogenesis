"""Count paired-nucleus ATAC insertions in fixed genomic windows with pysam."""

import argparse
import csv
import gzip
import math
import hashlib
import json
from pathlib import Path

import pysam

FRAGMENTS = {
    "GSM6339597": "GSM6339598",
    "GSM6339599": "GSM6339600",
    "GSM6339601": "GSM6339602",
    "GSM6339603": "GSM6339604",
}


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    p.add_argument("--fragments", type=Path, required=True)
    a = p.parse_args()
    out = a.root / "results/figure_source"
    with (out / "fragment_windows.tsv").open() as f:
        windows = list(csv.DictReader(f, delimiter="\t"))
    with gzip.open(a.root / "results/tables/matched_pairs.tsv.gz", "rt") as f:
        pairs = [
            r
            for r in csv.DictReader(f, delimiter="\t")
            if r["gate"] == "TSS_ge_3" and r["contrast"] == "3plus_vs_1"
        ]
    rows = []
    for gsm, fgsm in FRAGMENTS.items():
        ps = [r for r in pairs if r["gsm"] == gsm]
        group = {r[k + "_barcode"]: k for r in ps for k in ["high", "low"]}
        assert len(group) == 2 * len(ps)
        path = next(a.fragments.glob(fgsm + "*fragments.tsv.gz"))
        with pysam.TabixFile(str(path), index=str(path) + ".tbi") as tbx:
            for w in windows:
                lo, hi, bp = [int(w[k]) for k in ["start", "end", "bin_bp"]]
                n = math.ceil((hi - lo) / bp)
                counts = {g: [0] * n for g in ["high", "low"]}
                for line in tbx.fetch(w["chrom"], max(0, lo), hi):
                    f = line.split("\t")
                    g = group.get(f[3])
                    if g is None:
                        continue
                    # Cell Ranger fragment endpoints already include the Tn5 shift.
                    for cut in [int(f[1]), int(f[2])]:
                        if lo <= cut < hi:
                            counts[g][(cut - lo) // bp] += 1
                for g in ["low", "high"]:
                    for i, count in enumerate(counts[g]):
                        width = min(bp, hi - (lo + i * bp))
                        rows.append(
                            dict(
                                window_id=w["window_id"],
                                kind=w["kind"],
                                gsm=gsm,
                                group=g,
                                n_nuclei=len(ps),
                                bin_start=lo + i * bp,
                                bin_end=min(hi, lo + (i + 1) * bp),
                                relative_bp=lo + i * bp - int(w["center"]),
                                insertions=count,
                                insertions_per_100_nuclei_100bp=count
                                / len(ps)
                                * 100
                                * 100
                                / width,
                            )
                        )
        print(gsm, "completed", len(windows), "windows", flush=True)
    with gzip.open(out / "fragment_profiles.tsv.gz", "wt", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    identities = "\n".join(
        sorted(
            "|".join(r[k] for k in ["gsm", "high_barcode", "low_barcode"])
            for r in pairs
        )
    )
    manifest = dict(
        paired_barcode_sha256=hashlib.sha256(identities.encode()).hexdigest(),
        windows_sha256=hashlib.sha256(
            (out / "fragment_windows.tsv").read_text().encode()
        ).hexdigest(),
        endpoint_convention="Cell Ranger adjusted start and end; no additional Tn5 shift",
        fragment_counting="Each deduplicated fragment row contributes one insertion per endpoint",
    )
    (out / "fragment_profile_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n"
    )


if __name__ == "__main__":
    main()
