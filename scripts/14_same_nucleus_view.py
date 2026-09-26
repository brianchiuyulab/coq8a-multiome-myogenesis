"""Display matched nuclei and local RNA/ATAC measurements without new testing.

The UMAP uses the primary matched nuclei, label-blind variable RNA genes
(excluding COQ8A), and all-four-library common peaks in the fixed 221-gene
myogenesis search space. It is a descriptive embedding, not a trajectory or
cell-type annotation. The other panels show the same nuclei and the primary
Hallmark programme definitions.
"""

import argparse
from pathlib import Path

import h5py
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import sparse
from sklearn.decomposition import TruncatedSVD
from sklearn.feature_extraction.text import TfidfTransformer
from sklearn.preprocessing import StandardScaler
from umap import UMAP

from multiome_core import h5_for, read_barcodes


LIBRARIES = ("GSM6339597", "GSM6339599", "GSM6339601", "GSM6339603")
LIB_LABELS = ("L1 stem", "L1 diff", "L2 stem", "L2 diff")
LIB_COLORS = ("#6b83af", "#e9a34a", "#6b9f75", "#a578a7")
LOW_COLOR = "#3976a6"
HIGH_COLOR = "#c34354"
GENES_TO_DISPLAY = ("PAX7", "MYF5", "MYOD1", "MYOG", "MYF6", "MEF2C", "MEF2D")


def load_primary_data(h5_root: Path, tables: Path):
    pairs = pd.read_csv(tables / "matched_pairs.tsv.gz", sep="\t")
    pairs = pairs[(pairs.gate == "TSS_ge_3") & (pairs.contrast == "2plus_vs_1")]
    nuclei = pd.read_csv(tables / "qc_nuclei.tsv.gz", sep="\t")
    candidate = pd.read_csv(tables / "candidate_peak_gene.tsv.gz", sep="\t")
    candidate = candidate[candidate.n_libraries == 4]
    anchors = sorted(candidate.peak.unique())
    consensus = pd.read_csv(tables / "consensus_peak_map.tsv.gz", sep="\t")
    consensus = consensus.set_index("peak")
    rna_parts, atac_parts, metadata = [], [], []
    reference_rna_names = None
    for gsm in LIBRARIES:
        subset = pairs[pairs.gsm == gsm].sort_values("pair")
        high = subset[["high_barcode", "pair"]].rename(
            columns={"high_barcode": "barcode"}
        )
        low = subset[["low_barcode", "pair"]].rename(columns={"low_barcode": "barcode"})
        high = high.assign(coq_group="High")
        low = low.assign(coq_group="Low")
        cells = pd.concat([high, low], ignore_index=True).sort_values("barcode")
        cells = cells.merge(
            nuclei[nuclei.gsm == gsm][["barcode", "total_rna_umi", "COQ8A_umi"]],
            on="barcode",
            validate="one_to_one",
        )
        cells.insert(0, "gsm", gsm)
        cells["library"] = LIB_LABELS[LIBRARIES.index(gsm)]
        cells["source"] = subset.source.iloc[0]
        cells["stage"] = subset.stage.iloc[0]
        with h5py.File(h5_for(h5_root, gsm)) as h5:
            rna, atac, rna_names, peak_names = read_barcodes(h5, cells.barcode.tolist())
        if reference_rna_names is None:
            reference_rna_names = rna_names
        elif not np.array_equal(reference_rna_names, rna_names):
            raise ValueError(f"RNA feature order differs in {gsm}")
        local = {str(peak): index for index, peak in enumerate(peak_names)}
        mapped = [
            anchor if gsm == LIBRARIES[0] else consensus.at[anchor, gsm + "_peak"]
            for anchor in anchors
        ]
        missing = [peak for peak in mapped if peak not in local]
        if missing:
            raise ValueError(f"{gsm}: {len(missing)} common peaks absent from matrix")
        rna_parts.append(rna)
        atac_parts.append(atac[:, [local[peak] for peak in mapped]])
        metadata.append(cells)
    meta = pd.concat(metadata, ignore_index=True)
    rna = sparse.vstack(rna_parts, format="csr")
    atac = sparse.vstack(atac_parts, format="csr")
    if len(meta) != 2 * len(pairs):
        raise ValueError("Each primary matched pair must contribute two nuclei")
    return meta, rna, atac, reference_rna_names, anchors, candidate


def programme_scores(meta, rna, atac, rna_names, anchors, candidate, genes_path):
    genes = pd.read_csv(genes_path, sep="\t")
    hallmark = genes.loc[genes.Hallmark_Myogenesis, "gene"].tolist()
    rna_index = {str(name): index for index, name in enumerate(rna_names)}
    if not set(hallmark).issubset(rna_index):
        raise ValueError("Hallmark genes missing from RNA matrix")
    normalized = rna.multiply((10000 / meta.total_rna_umi.to_numpy())[:, None]).tocsr()
    normalized.data = np.log1p(normalized.data)
    meta["rna_hallmark"] = np.asarray(
        normalized[:, [rna_index[gene] for gene in hallmark]].mean(axis=1)
    ).ravel()
    anchor_index = {peak: index for index, peak in enumerate(anchors)}
    row, col, weight = [], [], []
    for column, gene in enumerate(hallmark):
        assigned = candidate.loc[candidate.gene == gene, "peak"].unique()
        if not len(assigned):
            raise ValueError(f"{gene}: no common ATAC region")
        row.extend(anchor_index[peak] for peak in assigned)
        col.extend([column] * len(assigned))
        weight.extend([1 / len(assigned)] * len(assigned))
    gene_peak_mean = sparse.csr_matrix(
        (weight, (row, col)), shape=(len(anchors), len(hallmark))
    )
    meta["atac_hallmark"] = np.asarray((atac @ gene_peak_mean).mean(axis=1)).ravel()
    for gene in GENES_TO_DISPLAY:
        if gene in rna_index:
            meta[gene + "_raw"] = rna[:, rna_index[gene]].toarray().ravel()
            meta[gene + "_lognorm"] = normalized[:, rna_index[gene]].toarray().ravel()
    return normalized


def joint_embedding(meta, normalized, atac, rna_names):
    means = np.asarray(normalized.mean(axis=0)).ravel()
    second = np.asarray(normalized.power(2).mean(axis=0)).ravel()
    variance = np.maximum(second - means**2, 0)
    eligible = np.array(
        [name != "COQ8A" and not name.startswith("MT-") for name in rna_names]
    ) & (means > 0.01)
    selected = np.flatnonzero(eligible)
    selected = selected[np.argsort(variance[selected])[-2000:]]
    rna_latent = TruncatedSVD(n_components=20, random_state=42).fit_transform(
        normalized[:, selected]
    )
    atac_tfidf = TfidfTransformer(norm="l2").fit_transform(atac)
    atac_latent = TruncatedSVD(n_components=20, random_state=42).fit_transform(
        atac_tfidf
    )
    rna_latent = StandardScaler().fit_transform(rna_latent)
    atac_latent = StandardScaler().fit_transform(atac_latent)
    joint = np.column_stack([rna_latent, atac_latent]) / np.sqrt(20)
    coordinates = UMAP(
        n_neighbors=35, min_dist=0.35, random_state=42, n_jobs=1
    ).fit_transform(joint)
    meta["UMAP1"] = coordinates[:, 0]
    meta["UMAP2"] = coordinates[:, 1]
    return len(selected)


def myod1_peak_display(meta, atac, anchors, candidate):
    peaks = sorted(candidate.loc[candidate.gene == "MYOD1", "peak"].unique())
    columns = [anchors.index(peak) for peak in peaks]
    high = meta.coq_group.eq("High").to_numpy()
    low = ~high
    high_fraction = np.asarray(atac[high][:, columns].mean(axis=0)).ravel()
    low_fraction = np.asarray(atac[low][:, columns].mean(axis=0)).ravel()
    coordinates = pd.Series(peaks).str.extract(r":(\d+)-(\d+)").astype(int)
    return pd.DataFrame(
        {
            "peak": peaks,
            "midpoint_kb_hg38": (coordinates[0] + coordinates[1]) / 2000,
            "high_open_fraction": high_fraction,
            "low_open_fraction": low_fraction,
        }
    ).sort_values("midpoint_kb_hg38")


def draw_violin(ax, meta, field, title):
    positions, values, colors = [], [], []
    for library_index, library in enumerate(LIB_LABELS):
        for group_index, group in enumerate(("Low", "High")):
            positions.append(library_index * 2.5 + group_index)
            values.append(
                meta.loc[
                    (meta.library == library) & (meta.coq_group == group), field
                ].to_numpy()
            )
            colors.append((LOW_COLOR, HIGH_COLOR)[group_index])
    violins = ax.violinplot(values, positions=positions, widths=0.78, showextrema=False)
    for body, color in zip(violins["bodies"], colors):
        body.set_facecolor(color)
        body.set_edgecolor(color)
        body.set_alpha(0.60)
    ax.scatter(positions, [np.median(value) for value in values], s=14, color="#172333")
    ax.set_xticks([index * 2.5 + 0.5 for index in range(4)], LIB_LABELS, fontsize=8)
    ax.set_title(title, loc="left", fontsize=10, fontweight="bold")
    ax.grid(axis="y", color="#e5e9ed", lw=0.6)
    ax.set_axisbelow(True)


def plot(meta, peaks, output: Path, n_variable_genes: int, n_common_peaks: int):
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 8,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "pdf.fonttype": 42,
        }
    )
    fig = plt.figure(figsize=(13.6, 12.0))
    layout = fig.add_gridspec(
        3, 2, height_ratios=[1.0, 0.95, 1.0], hspace=0.44, wspace=0.32
    )
    axes = [
        fig.add_subplot(layout[row, column]) for row in range(3) for column in range(2)
    ]
    for index, library in enumerate(LIB_LABELS):
        cells = meta[meta.library == library]
        axes[0].scatter(
            cells.UMAP1,
            cells.UMAP2,
            s=7,
            alpha=0.60,
            color=LIB_COLORS[index],
            label=library,
        )
    axes[0].legend(frameon=False, markerscale=2, fontsize=7, ncol=2, loc="best")
    axes[0].set_title(
        "A  Joint RNA/ATAC UMAP: source and stage", loc="left", fontweight="bold"
    )
    for group, color in (("Low", LOW_COLOR), ("High", HIGH_COLOR)):
        cells = meta[meta.coq_group == group]
        axes[1].scatter(
            cells.UMAP1,
            cells.UMAP2,
            s=7,
            alpha=0.52,
            color=color,
            label=f"COQ8A {group.lower()}",
        )
    axes[1].legend(frameon=False, markerscale=2, fontsize=8)
    axes[1].set_title(
        "B  Same UMAP: matched COQ8A groups", loc="left", fontweight="bold"
    )
    for ax in axes[:2]:
        ax.set_xlabel("UMAP 1")
        ax.set_ylabel("UMAP 2")
        ax.set_xticks([])
        ax.set_yticks([])
    draw_violin(axes[2], meta, "rna_hallmark", "C  Per-nucleus Hallmark myogenesis RNA")
    axes[2].set_ylabel("Mean log-normalized RNA")
    draw_violin(
        axes[3], meta, "atac_hallmark", "D  Per-nucleus nearby ATAC accessibility"
    )
    axes[3].set_ylabel("Mean open-peak fraction per gene")
    axes[2].scatter([], [], color=LOW_COLOR, label="COQ8A low: 1 UMI")
    axes[2].scatter([], [], color=HIGH_COLOR, label="COQ8A high: ≥2 UMI")
    axes[2].legend(frameon=False, fontsize=7, ncol=2, loc="upper left")
    ax = axes[4]
    rows = []
    for library in LIB_LABELS:
        for group in ("Low", "High"):
            cells = meta[(meta.library == library) & (meta.coq_group == group)]
            for gene in GENES_TO_DISPLAY:
                rows.append(
                    {
                        "library": library,
                        "group": group,
                        "gene": gene,
                        "mean": cells[gene + "_lognorm"].mean(),
                        "detected_pct": 100 * cells[gene + "_raw"].gt(0).mean(),
                    }
                )
    dots = pd.DataFrame(rows)
    dots["z"] = dots.groupby("gene")["mean"].transform(
        lambda values: (values - values.mean()) / max(values.std(ddof=0), 1e-8)
    )
    for row in dots.itertuples():
        y = LIB_LABELS.index(row.library) * 2 + (row.group == "High")
        x = GENES_TO_DISPLAY.index(row.gene)
        ax.scatter(
            x,
            y,
            s=12 + 3.2 * row.detected_pct,
            c=row.z,
            vmin=-1.8,
            vmax=1.8,
            cmap="RdBu_r",
            edgecolor="white",
            linewidth=0.35,
        )
    ax.set_xticks(
        range(len(GENES_TO_DISPLAY)), GENES_TO_DISPLAY, rotation=40, ha="right"
    )
    ax.set_yticks(
        range(8),
        [
            f"{library} · {group.lower()}"
            for library in LIB_LABELS
            for group in ("Low", "High")
        ],
        fontsize=7,
    )
    ax.invert_yaxis()
    ax.set_xlim(-0.5, len(GENES_TO_DISPLAY) - 0.5)
    ax.set_title(
        "E  Myogenic regulators: RNA detection and level",
        loc="left",
        fontweight="bold",
        pad=28,
    )
    ax.text(
        0,
        1.015,
        "Dot area: % nuclei detected; colour: gene-wise mean RNA z-score",
        transform=ax.transAxes,
        fontsize=7,
    )
    ax = axes[5]
    for peak in peaks.itertuples():
        ax.plot(
            [peak.midpoint_kb_hg38 - 0.45, peak.midpoint_kb_hg38 + 0.45],
            [100 * peak.low_open_fraction, 100 * peak.high_open_fraction],
            color="#aab3bd",
            lw=1.1,
        )
    ax.scatter(
        peaks.midpoint_kb_hg38 - 0.45,
        100 * peaks.low_open_fraction,
        s=20,
        color=LOW_COLOR,
        label="COQ8A low",
        zorder=3,
    )
    ax.scatter(
        peaks.midpoint_kb_hg38 + 0.45,
        100 * peaks.high_open_fraction,
        s=20,
        color=HIGH_COLOR,
        label="COQ8A high",
        zorder=3,
    )
    ax.axvline(17719.565, color="#32855c", lw=1, ls="--", label="MYOD1 TSS")
    ax.set_xlabel("chr11 coordinate (kb, hg38)")
    ax.set_ylabel("Nuclei with open peak (%)")
    ax.set_title(
        "F  MYOD1 nearby peaks: observed open fractions", loc="left", fontweight="bold"
    )
    ax.set_ylim(0, 29)
    ax.legend(frameon=False, fontsize=7, ncol=3, loc="upper right")
    ax.grid(axis="y", color="#e5e9ed", lw=0.6)
    fig.suptitle(
        "Same-nucleus view of COQ8A and myogenesis",
        y=0.995,
        fontsize=16,
        fontweight="bold",
    )
    fig.text(
        0.5,
        0.015,
        f"1,916 nuclei from 958 matched pairs · UMAP: {n_variable_genes:,} variable RNA genes + "
        f"{n_common_peaks:,} common local ATAC peaks · Descriptive display, not trajectory",
        ha="center",
        fontsize=8,
        color="#687582",
    )
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(
        output.with_suffix(".pdf"),
        bbox_inches="tight",
        metadata={"CreationDate": None, "ModDate": None},
    )
    fig.savefig(output.with_suffix(".png"), dpi=240, bbox_inches="tight")
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--h5-root", type=Path, required=True)
    parser.add_argument("--tables", type=Path, required=True)
    parser.add_argument("--genes", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    meta, rna, atac, rna_names, anchors, candidate = load_primary_data(
        args.h5_root, args.tables
    )
    normalized = programme_scores(
        meta, rna, atac, rna_names, anchors, candidate, args.genes
    )
    n_variable_genes = joint_embedding(meta, normalized, atac, rna_names)
    peaks = myod1_peak_display(meta, atac, anchors, candidate)
    meta.to_csv(
        args.tables / "same_nucleus_display.tsv.gz",
        sep="\t",
        index=False,
        compression={"method": "gzip", "mtime": 0},
    )
    peaks.to_csv(args.tables / "myod1_peak_display.tsv", sep="\t", index=False)
    plot(meta, peaks, args.out, n_variable_genes, len(anchors))
    print(
        f"Wrote same-nucleus view for {len(meta)} nuclei and {len(anchors)} common peaks"
    )


if __name__ == "__main__":
    main()
