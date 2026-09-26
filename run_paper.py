"""Reproduce the manuscript analysis and figures from processed public inputs."""

import argparse
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys

ROOT = Path(__file__).resolve().parent


def run(script, *args):
    subprocess.run(
        [sys.executable, str(ROOT / "scripts" / script), *map(str, args)],
        cwd=ROOT,
        check=True,
    )


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--mode", choices=["figures", "analysis"], default="figures")
    p.add_argument("--h5-root", type=Path)
    p.add_argument("--gtf", type=Path)
    p.add_argument("--external-raw", type=Path)
    p.add_argument("--external-work", type=Path)
    p.add_argument("--reuse-external", action="store_true")
    p.add_argument(
        "--fragments",
        type=Path,
        help="Regenerate figure fragment profiles; requires pysam",
    )
    a = p.parse_args()
    if a.fragments and a.gtf is None:
        p.error("--fragments requires --gtf for TSS profiles")
    t = ROOT / "results/tables"
    e = ROOT / "results/temporal"
    r = ROOT / "reference"
    if a.mode == "analysis":
        for field in ["h5_root", "gtf", "external_raw", "external_work"]:
            if getattr(a, field) is None:
                p.error("analysis requires --" + field.replace("_", "-"))
        h, g = a.h5_root.resolve(), a.gtf.resolve()
        run("01_extract_nuclei.py", "--h5-root", h, "--out", t / "nuclei.tsv.gz")
        run(
            "02_qc_and_matching.py",
            "--nuclei",
            t / "nuclei.tsv.gz",
            "--fragment-qc",
            r / "fragment_qc",
            "--out",
            t,
        )
        run(
            "03_define_regions.py",
            "--h5-root",
            h,
            "--gtf",
            g,
            "--gene-list",
            r / "myogenesis_221_gene_sources.tsv",
            "--blacklist",
            r / "hg38-blacklist.v2.bed.gz",
            "--out",
            t,
        )
        run(
            "04_gene_programme_effects.py",
            "--h5-root",
            h,
            "--tables",
            t,
            "--genes",
            r / "myogenesis_221_gene_sources.tsv",
        )
        run("05_peak_effects.py", "--h5-root", h, "--tables", t)
        run("19_fetch_temporal_reference.py", "--out", a.external_raw)
        if not a.reuse_external:
            for exp in [1, 2]:
                run(
                    "20_prepare_temporal_reference.py",
                    "--raw",
                    a.external_raw,
                    "--tables",
                    t,
                    "--gtf",
                    g,
                    "--out",
                    a.external_work,
                    "--experiment",
                    exp,
                )
        run(
            "21_define_temporal_regions.py",
            "--prepared",
            a.external_work,
            "--tables",
            t,
            "--out",
            e,
        )
        membership_paths = [
            e / "temporal_set_memberships.tsv.gz",
            e / "temporal_set_definitions.tsv",
        ]
        membership_hashes = {
            path.name: hashlib.sha256(path.read_bytes()).hexdigest()
            for path in membership_paths
        }
        (e / "frozen_membership_sha256.json").write_text(
            json.dumps(membership_hashes, indent=2) + "\n", encoding="utf-8"
        )
        run(
            "22_test_temporal_accessibility.py",
            "--h5-root",
            h,
            "--tables",
            t,
            "--temporal",
            e,
            "--out",
            e,
        )
        for path in membership_paths:
            assert (
                hashlib.sha256(path.read_bytes()).hexdigest()
                == membership_hashes[path.name]
            )
        for pattern in [
            "*cells.tsv.gz",
            "*ingest_audit.tsv",
            "candidate_hg19_mapping.tsv",
            "myogenic_marker_promoters_hg19.tsv",
        ]:
            for path in a.external_work.glob(pattern):
                shutil.copy2(path, e / path.name)
        shutil.copy2(
            a.external_raw / "download_manifest.json", e / "download_manifest.json"
        )
        run("27_temporal_single_peaks.py", "--tables", t, "--temporal", e)
        run("32_bidirectional_peak_screen.py", "--tables", t, "--temporal", e)
        run("37_functional_subsets.py", "freeze")
        run("37_functional_subsets.py", "test", "--h5-root", h)
        run(
            "07_peak_gene_links.py",
            "--h5-root",
            h,
            "--tables",
            t,
            "--gate",
            "TSS_ge_3",
            "--contrast",
            "3plus_vs_1",
        )
        run(
            "30_middle_cis_links.py",
            "--h5-root",
            h,
            "--tables",
            t,
            "--temporal",
            e,
            "--gtf",
            g,
            "--region-family",
            "opening-closing",
        )
        run("33_dynamic_evidence.py", "--tables", t, "--temporal", e)
        run("38_design_inventory.py")
        run("39_prepare_figure_sources.py", "--gtf", g)
    if a.fragments:
        run("40_extract_fragment_windows.py", "--fragments", a.fragments)
        run(
            "14_tss_fragment_profiles.py",
            "--fragments-dir",
            a.fragments,
            "--gtf",
            a.gtf,
            "--genes",
            r / "myogenesis_221_gene_sources.tsv",
            "--pairs",
            t / "matched_pairs.tsv.gz",
            "--out",
            t / "tss_fragment_profile_221.tsv.gz",
        )
    run("41_make_submission_figures.py")
    run("36_plot_tss_heatmap.py", "--tables", t, "--out", ROOT / "figures/supplement")
    run("42_validate_submission.py")


if __name__ == "__main__":
    main()
