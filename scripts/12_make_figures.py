"""Build the three main figures and QC supplement from released tables.

All plot values are read from committed analysis tables. The complete primary
test families precede the exploratory MYOD1 ranking in the analysis pipeline.
"""

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

from figure1_tss import make_figure as figure_one
from figure2_screen import make_figure as figure_two
from figure3_myod1 import make_figure as figure_three
from supplement_qc_figure import make_figure as qc_supplement


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tables", type=Path, required=True)
    parser.add_argument("--figures", type=Path, required=True)
    args = parser.parse_args()

    main_dir = args.figures / "main"
    supplement_dir = args.figures / "supplement"
    figure_one(args.tables, main_dir / "Figure_1_Multiome_and_TSS")
    figure_two(args.tables, main_dir / "Figure_2_Global_Screen_and_Ranking")
    figure_three(args.tables, main_dir / "Figure_3_MYOD1_Locus_and_Sensitivity")
    qc_supplement(args.tables, supplement_dir / "Supplementary_Figure_QC")
    print("Wrote three main figures and one QC supplement", flush=True)


if __name__ == "__main__":
    main()
