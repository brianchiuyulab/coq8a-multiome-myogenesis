"""Reproduce the external-time-defined COQ8A accessibility analysis."""

import argparse
import hashlib
import json
import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent


def file_digest(path):
    data = path.read_bytes()
    if path.suffix == ".tsv":
        data = data.replace(b"\r\n", b"\n")
    return hashlib.sha256(data).hexdigest()


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument(
        "--raw", type=Path, required=True, help="GSE109828 download directory"
    )
    p.add_argument(
        "--gtf", type=Path, required=True, help="GENCODE v48 hg38 annotation.gtf.gz"
    )
    p.add_argument(
        "--h5-root",
        type=Path,
        required=True,
        help="GSE208248 four processed multiome H5 files",
    )
    p.add_argument(
        "--work", type=Path, required=True, help="External sparse matrices and cell QC"
    )
    p.add_argument("--tables", type=Path, default=ROOT / "results/tables")
    p.add_argument("--out", type=Path, default=ROOT / "results/temporal")
    p.add_argument("--figures", type=Path, default=ROOT / "figures/temporal")
    p.add_argument(
        "--reuse-prepared",
        action="store_true",
        help="Reuse completed external matrices; rebuild time sets and effects",
    )
    a = p.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)

    def run(script, *args):
        subprocess.run(
            [sys.executable, str(ROOT / "scripts" / script), *map(str, args)],
            check=True,
        )

    run("19_fetch_temporal_reference.py", "--out", a.raw)
    if not a.reuse_prepared:
        for exp in [1, 2]:
            run(
                "20_prepare_temporal_reference.py",
                "--raw",
                a.raw,
                "--tables",
                a.tables,
                "--gtf",
                a.gtf,
                "--out",
                a.work,
                "--experiment",
                exp,
            )
    run(
        "21_define_temporal_regions.py",
        "--prepared",
        a.work,
        "--tables",
        a.tables,
        "--out",
        a.out,
    )
    frozen = {
        name: file_digest(a.out / name)
        for name in ["temporal_set_definitions.tsv", "temporal_set_memberships.tsv.gz"]
    }
    run(
        "22_test_temporal_accessibility.py",
        "--h5-root",
        a.h5_root,
        "--tables",
        a.tables,
        "--temporal",
        a.out,
        "--out",
        a.out,
    )
    for name, checksum in frozen.items():
        assert file_digest(a.out / name) == checksum
    (a.out / "frozen_membership_sha256.json").write_text(
        json.dumps(frozen, indent=2) + "\n", encoding="utf-8"
    )
    for pattern in [
        "*ingest_audit.tsv",
        "*cells.tsv.gz",
        "candidate_hg19_mapping.tsv",
        "myogenic_marker_promoters_hg19.tsv",
    ]:
        for path in a.work.glob(pattern):
            shutil.copy2(path, a.out / path.name)
    shutil.copy2(a.raw / "download_manifest.json", a.out / "download_manifest.json")
    run("24_temporal_specificity.py", "--temporal", a.out)
    run(
        "25_temporal_rna_qc.py",
        "--h5-root",
        a.h5_root,
        "--tables",
        a.tables,
        "--temporal",
        a.out,
        "--scrublet",
        ROOT / "reference/rna_scrublet_qc.tsv.gz",
    )
    run("23_plot_temporal_accessibility.py", "--temporal", a.out, "--out", a.figures)
    run(
        "26_validate_temporal.py",
        "--h5-root",
        a.h5_root,
        "--tables",
        a.tables,
        "--temporal",
        a.out,
    )


if __name__ == "__main__":
    main()
