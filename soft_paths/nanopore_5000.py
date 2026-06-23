import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

# Shared data directory - unchanged
samples_dir = Path("/gpfs01/share/BioinfMSc/Matt_Projects/samples")

# Output directory - derived relative to this script's location
# Assumes this script lives at: <project_root>/nanopore_profile/nanopore_5000.py
# Output goes to:               <project_root>/scatter_plot_nanopore/5000/
SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_DIR = SCRIPT_DIR.parent
outdir = PROJECT_DIR / "scatter_plot_nanopore" / "5000"
outdir.mkdir(parents=True, exist_ok=True)

# autosomes only
chromosomes = [f"chr{i}" for i in range(1, 23)]

def bin_by_estimate_reads(values, bin_width, target_reads):
    x_list = []
    y_list = []

    group = []
    reads_sum = 0
    start = 0

    for i, cnv_value in enumerate(values):
        estimate_reads = cnv_value / 2 * 100
        group.append(cnv_value)
        reads_sum += estimate_reads
        if reads_sum >= target_reads:
            y_list.append(np.median(group))
            x_list.append(((start + i + 1) / 2) * bin_width)
            group = []
            reads_sum = 0
            start = i + 1
    if group:
        x_list.append(((start + len(values)) / 2) * bin_width)
        y_list.append(np.median(group))
    return np.array(x_list), np.array(y_list)


#### loop through all sample folders
for sample_dir in sorted(samples_dir.iterdir()):

    if not sample_dir.is_dir():
        continue

    sample = sample_dir.name

    cnv_file = sample_dir / "nanopore" / "CNV.npy"
    cnv_dict_file = sample_dir / "nanopore" / "CNV_dict.npy"

    output = outdir / f"{sample}_nanopore_merged_reads.png"

    cnv = np.load(cnv_file, allow_pickle=True).item()
    cnv_dict = np.load(cnv_dict_file, allow_pickle=True).item()

    bin_width = cnv_dict["bin_width"]

    target_reads = 5000

    merged_bin_width = target_reads / 100 * bin_width
    print(
        sample,
        "original_bin_width:", bin_width,
        "target_reads:", target_reads,
        "merged_bin_width:", round(merged_bin_width)
    )

    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(18, 4))

    total = 0
    xticks = []
    xticklabels = []

    for contig in chromosomes:
        values = np.array(cnv[contig], dtype=float)
        chromo_length = len(values) * bin_width

        x, y = bin_by_estimate_reads(values, bin_width, target_reads)

        ax.scatter(
            x=x + total,
            y=y,
            s=0.1
        )

        xticks.append(total + chromo_length / 2)
        xticklabels.append(contig.replace("chr", ""))

        total += chromo_length

    ax.set_xticks(xticks)
    ax.set_xticklabels(xticklabels)

    ax.set_ylim((0, 8))
    ax.set_xlim((0, total))

    ax.set_xlabel("Chromosome")
    ax.set_ylabel("Copy number")
    ax.set_title(f"Nanopore CNV plot: {sample}")

    fig.savefig(output, dpi=300)
    plt.close(fig)
