"""Reproduce the primary search, exploratory locus analysis and figures.

Inputs are the four GEO 10x Multiome H5 matrices, four matching barcode metric
files, GENCODE v48 hg38 GTF and versioned reference tables. The complete
primary ATAC test families are written before any candidate ranking.
"""

import argparse
from pathlib import Path
import shutil
import subprocess
import sys


ROOT = Path(__file__).resolve().parent


SAMPLES = ("GSM6339597", "GSM6339599", "GSM6339601", "GSM6339603")


def require_input_files(h5_root: Path, gtf: Path) -> None:
    if not gtf.is_file():
        raise FileNotFoundError(f"GENCODE GTF not found: {gtf}")
    for sample in SAMPLES:
        for suffix in ("filtered_feature_bc_matrix.h5", "per_barcode_metrics.csv.gz"):
            matches = list(h5_root.glob(f"{sample}_*{suffix}"))
            if len(matches) != 1:
                raise FileNotFoundError(
                    f"Expected one {sample} {suffix} under {h5_root}; found {len(matches)}"
                )


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument(
        "--h5-root",
        type=Path,
        required=True,
        help="Directory with four GEO filtered_feature_bc_matrix.h5 files and barcode metrics",
    )
    p.add_argument(
        "--gtf", type=Path, required=True, help="GENCODE v48 hg38 annotation.gtf.gz"
    )
    p.add_argument(
        "--out-root",
        type=Path,
        default=ROOT,
        help="Output repository root; default is this directory",
    )
    p.add_argument(
        "--fragments-dir",
        type=Path,
        help="Optional indexed ATAC fragments; rebuild the 221-gene TSS profile instead of copying the versioned profile",
    )
    args = p.parse_args()
    t = args.out_root / "results" / "tables"
    f = args.out_root / "figures"
    t.mkdir(parents=True, exist_ok=True)
    f.mkdir(parents=True, exist_ok=True)
    r = ROOT / "reference"
    h = args.h5_root.resolve()
    g = args.gtf.resolve()
    require_input_files(h, g)
    primary_steps = [
        ("01_extract_nuclei.py", "--h5-root", h, "--out", t / "nuclei.tsv.gz"),
        (
            "02_qc_and_matching.py",
            "--nuclei",
            t / "nuclei.tsv.gz",
            "--fragment-qc",
            r / "fragment_qc",
            "--out",
            t,
        ),
        (
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
        ),
        (
            "04_gene_programme_effects.py",
            "--h5-root",
            h,
            "--tables",
            t,
            "--genes",
            r / "myogenesis_221_gene_sources.tsv",
        ),
        ("05_peak_effects.py", "--h5-root", h, "--tables", t),
        ("06_primary_decision.py", "--tables", t),
    ]
    exploratory_steps = [
        (
            "07_peak_gene_links.py",
            "--h5-root",
            h,
            "--tables",
            t,
            "--motif-scores",
            r / "jaspar2024_peak_scores.tsv.gz",
        ),
        (
            "07_peak_gene_links.py",
            "--h5-root",
            h,
            "--tables",
            t,
            "--motif-scores",
            r / "jaspar2024_peak_scores.tsv.gz",
            "--gate",
            "TSS_ge_2",
        ),
        ("08_candidate_ranking.py", "--tables", t),
        ("09_link_gate_reconciliation.py", "--tables", t),
        ("10_locus_sensitivity.py", "--h5-root", h, "--tables", t),
        (
            "11_doublet_sensitivity.py",
            "--tables",
            t,
            "--scrublet",
            r / "rna_scrublet_qc.tsv.gz",
        ),
        ("12_make_figures.py", "--tables", t, "--figures", f),
        ("13_validate_release.py", "--tables", t, "--figures", f),
    ]
    for section, steps in [
        ("PRIMARY", primary_steps),
        ("EXPLORATORY AND FIGURES", exploratory_steps),
    ]:
        print(f"\n{section}", flush=True)
        for name, *arguments in steps:
            if name == "12_make_figures.py":
                profile = t / "tss_fragment_profile_221.tsv.gz"
                if args.fragments_dir:
                    cmd = [
                        sys.executable,
                        str(ROOT / "scripts" / "14_tss_fragment_profiles.py"),
                        "--fragments-dir",
                        str(args.fragments_dir.resolve()),
                        "--gtf",
                        str(g),
                        "--genes",
                        str(r / "myogenesis_221_gene_sources.tsv"),
                        "--pairs",
                        str(t / "matched_pairs.tsv.gz"),
                        "--out",
                        str(profile),
                    ]
                    print("RUN 14_tss_fragment_profiles.py", flush=True)
                    subprocess.run(cmd, check=True)
                elif not profile.is_file():
                    source = ROOT / "results" / "tables" / profile.name
                    if not source.is_file():
                        raise FileNotFoundError(
                            f"Versioned TSS profile not found: {source}"
                        )
                    shutil.copyfile(source, profile)
            cmd = [sys.executable, str(ROOT / "scripts" / name), *map(str, arguments)]
            print("RUN", name, flush=True)
            subprocess.run(cmd, check=True)


if __name__ == "__main__":
    main()
