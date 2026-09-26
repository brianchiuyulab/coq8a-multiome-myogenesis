"""Annotate fixed middle peaks with independent human MYOD1 ChIP peak calls."""

import argparse
import hashlib
import json
import urllib.request
from pathlib import Path
import pandas as pd

p = argparse.ArgumentParser(description=__doc__)
p.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
p.add_argument("--out", type=Path, required=True)
a = p.parse_args()
a.out.mkdir(parents=True, exist_ok=True)
members = pd.read_csv(
    a.root / "results/temporal/temporal_set_memberships.tsv.gz", sep="\t"
)
peaks = set(members.loc[members.set_id.eq("P2_O50_Middle_24to48h"), "peak"])
mapping = pd.read_csv(a.root / "results/temporal/candidate_hg19_mapping.tsv", sep="\t")
mapping = mapping[mapping.peak.isin(peaks)].copy()
assert len(mapping) == 34 and mapping.mapping_pass.all()
sources = [
    ("myoblast", "GSM1218849", "GSM1218849_MB135GMMD.peak.txt.gz"),
    ("myotube_72h", "GSM1218850", "GSM1218850_MB135DMMD.peak.txt.gz"),
]
rows = []
manifest = []
for stage, gsm, name in sources:
    url = f"https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM1218nnn/{gsm}/suppl/{name}"
    file = a.out / name
    if not file.exists():
        file.write_bytes(urllib.request.urlopen(url, timeout=30).read())
    chip = pd.read_csv(
        file,
        sep="\t",
        header=None,
        names=[
            "chrom",
            "start",
            "end",
            "name",
            "score",
            "strand",
            "signal",
            "minus_log10_p",
            "q",
        ],
    )
    manifest.append(
        dict(
            gsm=gsm,
            stage=stage,
            url=url,
            assembly="hg19",
            sha256=hashlib.sha256(file.read_bytes()).hexdigest(),
            n_chip_peaks=len(chip),
        )
    )
    for r in mapping.itertuples():
        sub = chip[
            (chip.chrom == r.chrom19) & (chip.start < r.end19) & (chip.end > r.start19)
        ]
        overlap = (sub.end.clip(upper=r.end19) - sub.start.clip(lower=r.start19)).sum()
        rows.append(
            dict(
                peak=r.peak,
                stage=stage,
                hg19_region=f"{r.chrom19}:{r.start19}-{r.end19}",
                overlap_chip_peaks=len(sub),
                overlap_bp_sum=int(overlap),
                max_target_overlap_fraction=(
                    float(
                        (
                            (
                                sub.end.clip(upper=r.end19)
                                - sub.start.clip(lower=r.start19)
                            )
                            / (r.end19 - r.start19)
                        ).max()
                    )
                    if len(sub)
                    else 0
                ),
            )
        )
pd.DataFrame(rows).to_csv(a.out / "middle_MYOD1_overlap.tsv", sep="\t", index=False)
(a.out / "MYOD1_source_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
print(
    pd.DataFrame(rows)
    .groupby("stage")
    .overlap_chip_peaks.apply(lambda x: int((x > 0).sum()))
)
