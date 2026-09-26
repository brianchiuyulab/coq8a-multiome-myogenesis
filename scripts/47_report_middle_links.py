"""Package the middle-region association follow-up and render its source-based figure."""

import argparse
import hashlib
import json
import shutil
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import sparse
from statsmodels.stats.multitest import multipletests


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--root", type=Path, default=Path(__file__).resolve().parents[1]
    )
    parser.add_argument("--work", type=Path)
    args = parser.parse_args()
    out = args.root / "results/link_followup"
    out.mkdir(parents=True, exist_ok=True)
    if args.work:
        files = [
            "candidate_pairs.tsv",
            "programme_genes.txt",
            "ambiguous_rna_symbols_excluded.txt",
            "correlations_by_library.tsv.gz",
            "correlations_combined.tsv",
            "count_models_by_library.tsv",
            "count_models_combined.tsv",
            "scent_style_bootstrap_followup.tsv",
            "scent_source.json",
            "validation.json",
            "external/middle_MYOD1_overlap.tsv",
            "external/MYOD1_source_manifest.json",
        ]
        for file in files:
            shutil.copy2(args.work / file, out / Path(file).name)
        source = args.work / "SCENTfunctions.R"
        manifest = json.loads((out / "scent_source.json").read_text())
        manifest["source_sha256"] = hashlib.sha256(source.read_bytes()).hexdigest()
        manifest["source_url"] = (
            "https://raw.githubusercontent.com/immunogenomics/SCENT/"
            + manifest["commit"]
            + "/R/SCENTfunctions.R"
        )
        (out / "scent_source.json").write_text(json.dumps(manifest, indent=2) + "\n")
        diagnostics = []
        genes = (args.work / "genes.txt").read_text().splitlines()
        for meta_file in sorted(args.work.glob("GSM*_meta.tsv")):
            gsm = meta_file.name.split("_")[0]
            meta = pd.read_csv(meta_file, sep="\t")
            rna = sparse.load_npz(args.work / (gsm + "_rna.npz"))
            for gene in ["COQ8A", "LDB3", "TNNC1", "MRAS", "SVIL", "MYOD1", "MYOG"]:
                y = rna[genes.index(gene)].toarray().ravel()
                diagnostics.append(
                    dict(
                        gsm=gsm,
                        line=meta.source.iloc[0],
                        stage=meta.stage.iloc[0],
                        gene=gene,
                        n_nuclei=len(meta),
                        detected_fraction=(y > 0).mean(),
                        mean_log1p_CP10k=np.log1p(
                            y / meta.total_rna_umi * 10000
                        ).mean(),
                    )
                )
        pd.DataFrame(diagnostics).to_csv(
            out / "RNA_by_library_descriptive.tsv", sep="\t", index=False
        )

    counts = pd.read_csv(out / "count_models_combined.tsv", sep="\t")
    for model, group in counts.groupby("model"):
        assert len(group) == 107, (model, len(group))
        assert np.allclose(group.q, multipletests(group.p, method="fdr_bh")[1])
    effects = pd.read_csv(
        args.root / "results/temporal/bidirectional_peak_effects.tsv", sep="\t"
    )
    effects = effects.query(
        "gate=='TSS_ge_3' and contrast=='3plus_vs_1' and phase=='Middle'"
    )
    assert len(effects) == 34
    links = pd.read_csv(out / "correlations_combined.tsv", sep="\t").query(
        "population=='all_qc' and detection_fraction==0.05 and model=='technical_state_coq'"
    )
    focal = counts.query(
        "model=='technical_state_coq' and gene in ['LDB3','TNNC1','MRAS','SVIL']"
    ).copy()
    focal = focal.rename(
        columns={"q": "q_link_107", "p": "p_link", "RNA_ratio": "RNA_open_closed_ratio"}
    )
    focal = focal.merge(
        links[["peak", "gene", "r", "loo_r_min", "loo_r_max"]], on=["peak", "gene"]
    )
    focal = focal.merge(
        effects[["peak", "fold_open", "p_exact", "q_within_phase_peaks"]], on="peak"
    )
    focal = focal.rename(
        columns={
            "fold_open": "ATAC_COQ8A_high_low_FC",
            "p_exact": "p_ATAC",
            "q_within_phase_peaks": "q_ATAC_34",
        }
    )
    focal.to_csv(out / "focal_evidence.tsv", sep="\t", index=False)

    models = ["technical", "technical_coq", "technical_state", "technical_state_coq"]
    labels = ["Technical covariates", "+ COQ8A", "+ Myogenesis score", "+ Both"]
    colors = ["#64859E", "#3C6B89", "#E4A35B", "#B96723"]
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 10,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 3.9), sharey=True)
    figure_rows = []
    for ax, gene, letter in zip(axes, ["LDB3", "TNNC1"], ["A", "B"]):
        sub = counts[counts.gene.eq(gene)].set_index("model").loc[models]
        assert len(sub) == 4
        for y, (model, row), color in zip(range(4), sub.iterrows(), colors):
            ax.errorbar(
                row.RNA_ratio,
                y,
                xerr=[[row.RNA_ratio - row.CI_low], [row.CI_high - row.RNA_ratio]],
                fmt="o",
                color=color,
                capsize=3,
                markersize=6,
                linewidth=1.7,
            )
            ax.text(
                2.03, y, f"{row.RNA_ratio:.2f}   {row.q:.2g}", va="center", fontsize=9
            )
            figure_rows.append(dict(model=model, **row.to_dict()))
        ax.axvline(1, color="#777777", linestyle="--", linewidth=1)
        ax.set(
            xlim=(0.85, 2.38),
            ylim=(3.65, -0.85),
            xticks=[1, 1.25, 1.5, 1.75],
            xlabel="RNA count ratio: peak open / closed",
            yticks=range(4),
            yticklabels=labels,
        )
        ax.text(2.03, -0.58, "Ratio    q", fontsize=9, weight="bold")
        ax.set_title(gene, weight="bold", pad=12)
        ax.text(-0.12, 1.06, letter, transform=ax.transAxes, weight="bold", fontsize=14)
        ax.spines[["top", "right", "left"]].set_visible(False)
        ax.tick_params(axis="y", length=0)
    fig.subplots_adjust(left=0.20, right=0.985, bottom=0.22, top=0.80, wspace=0.12)
    folder = args.root / "figures/link_followup"
    folder.mkdir(parents=True, exist_ok=True)
    for suffix in ["png", "pdf"]:
        fig.savefig(folder / ("Middle_peak_RNA_associations." + suffix), dpi=220)
    pd.DataFrame(figure_rows).to_csv(out / "figure_source.tsv", sep="\t", index=False)
    print(
        focal[
            [
                "gene",
                "RNA_open_closed_ratio",
                "q_link_107",
                "r",
                "ATAC_COQ8A_high_low_FC",
                "q_ATAC_34",
            ]
        ].to_string(index=False)
    )


if __name__ == "__main__":
    main()
