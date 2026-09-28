"""Entry point for the current Day-0 analysis, audit, and publication figures."""

import argparse
from pathlib import Path
import subprocess
import sys


ROOT = Path(__file__).resolve().parent


def run(script, *arguments):
    subprocess.run(
        [sys.executable, str(ROOT / "scripts" / script), *map(str, arguments)],
        cwd=ROOT,
        check=True,
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--mode", choices=["figures", "analysis", "audit"], default="figures"
    )
    parser.add_argument("--h5", type=Path)
    parser.add_argument("--work", type=Path)
    parser.add_argument("--gtf", type=Path)
    parser.add_argument("--fragments", type=Path)
    parser.add_argument(
        "--rebuild-upstream",
        action="store_true",
        help="Re-extract depths, reapply archived fragment QC and rebuild four-library peak mapping",
    )
    parser.add_argument(
        "--permutations",
        action="store_true",
        help="Also rerun 2,000,000 permutations in audit mode",
    )
    args = parser.parse_args()
    if args.mode != "figures" and (args.h5 is None or args.work is None):
        parser.error("analysis/audit require --h5 and --work")
    if (
        args.mode == "audit" or args.rebuild_upstream or args.fragments
    ) and args.gtf is None:
        parser.error("audit, --rebuild-upstream and --fragments require --gtf")
    if args.rebuild_upstream:
        tables = ROOT / "results/tables"
        run(
            "01_extract_nuclei.py",
            "--h5-root",
            args.h5,
            "--out",
            tables / "nuclei.tsv.gz",
        )
        run(
            "02_qc_and_matching.py",
            "--nuclei",
            tables / "nuclei.tsv.gz",
            "--fragment-qc",
            ROOT / "reference/fragment_qc",
            "--out",
            tables,
        )
        run(
            "03_define_regions.py",
            "--h5-root",
            args.h5,
            "--gtf",
            args.gtf,
            "--gene-list",
            ROOT / "reference/myogenesis_221_gene_sources.tsv",
            "--blacklist",
            ROOT / "reference/hg38-blacklist.v2.bed.gz",
            "--out",
            tables,
        )
    if args.mode == "analysis":
        run("98_recompute_fixed_statistics.py", "--h5", args.h5, "--work", args.work)
        run(
            "92_fixed_day0_workflow.py",
            "--h5",
            args.h5,
            "--work",
            args.work,
            "--statistics",
            ROOT / "results/fixed_day0_D206",
        )
    elif args.mode == "audit":
        options = ["--work", args.work, "--h5", args.h5, "--gtf", args.gtf]
        if args.permutations:
            options.append("--permutations")
        run("97_audit_fixed_day0.py", *options)
        return
    if args.fragments:
        run(
            "94_fixed_fragment_profiles.py",
            "--root",
            ROOT,
            "--fragments",
            args.fragments,
            "--gtf",
            args.gtf,
        )
    run("93_plot_fixed_day0.py")
    run("96_plot_fragment_panels.py")
    run("99_validate_fixed_tables.py")


if __name__ == "__main__":
    main()
