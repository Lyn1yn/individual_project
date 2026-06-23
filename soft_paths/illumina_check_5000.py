import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path

# Shared data directory - unchanged
samples_dir = Path("/gpfs01/share/BioinfMSc/Matt_Projects/samples")

# Input bin width file and output directory - derived relative to this script's location
# Assumes this script lives at: <project_root>/illumina_profile/illumina_check_5000.py
# bin_width_file is expected at: <project_root>/scatter_plot_nanopore/all_bin_width.txt
# Output goes to:                <project_root>/scatter_plot_illumina/5000/
SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_DIR = SCRIPT_DIR.parent
bin_width_file = PROJECT_DIR / "scatter_plot_nanopore" / "all_bin_width.txt"
outdir = PROJECT_DIR / "scatter_plot_illumina" / "5000"
outdir.mkdir(parents=True, exist_ok=True)

# use the same target reads as nanopore script
target_np_reads = 5000

bin_width_df = pd.read_csv(
    bin_width_file,
    sep=r"\s+",
    header=None,
    names=["sample", "bin_width"],
)


def merge_by_reads(sub, reads, target_reads):
    x_list = []
    y_list = []

    reads_sum = 0
    x_group = []
    cn_group = []

    for index, row in sub.iterrows():
        reads_sum += row[reads]
        x_group.append(row["start"])
        cn_group.append(row["copy_number"])

        if reads_sum >= target_reads:
            x_list.append(pd.Series(x_group).median())
            y_list.append(pd.Series(cn_group).median())

            reads_sum = 0
            x_group = []
            cn_group = []

    if x_group:
        x_list.append(pd.Series(x_group).median())
        y_list.append(pd.Series(cn_group).median())
    return x_list, y_list


for sample, bin_width in zip(bin_width_df["sample"], bin_width_df["bin_width"]):

    file = samples_dir / sample / "illumina" / f"{sample}.tumor.target.counts.gz"

    if not Path(file).exists():
        print(f"{sample}: illumina file not found, skipping")
        continue

    df = pd.read_csv(file, sep="\t", comment="#")

    reads = df.columns[4]

    autosomes = [f"chr{i}" for i in range(1, 23)]
    df = df[df["contig"].isin(autosomes)].copy()

    median_reads = df[reads].median()
    df["copy_number"] = 2 * df[reads] / median_reads

    df["chrom_num"] = df["contig"].str.replace("chr", "").astype(int)
    df = df.sort_values(["chrom_num", "start"]).copy()

    illumina_bin_size = (df["stop"] - df["start"]).median()

    merged_bin_width = target_np_reads / 100 * bin_width
    target_reads = merged_bin_width / illumina_bin_size * median_reads
    target_reads = round(target_reads)

    approx_illumina_bins_per_point = target_reads / median_reads
    illumina_approx_merged_bin_width = approx_illumina_bins_per_point * illumina_bin_size

    print(
        f"{sample} "
        f"nanopore_original_bin_width: {bin_width:.0f} "
        f"nanopore_target_reads: {target_np_reads} "
        f"target_merged_bin_width: {merged_bin_width:.0f} "
        f"illumina_original_bin_size: {illumina_bin_size:.0f} "
        f"illumina_median_reads_per_bin: {median_reads:.2f} "
        f"illumina_target_reads: {target_reads} "
        f"approx_illumina_bins_per_point: {approx_illumina_bins_per_point:.1f} "
        f"illumina_approx_merged_bin_width: {illumina_approx_merged_bin_width:.0f}"
    )

    output = outdir / f"{sample}_illumina_merged_reads.png"

    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(18, 4))

    total = 0
    xticks = []
    xticklabels = []

    for contig in autosomes:
        sub = df[df["contig"] == contig].copy()
        chromo_length = sub["stop"].max()

        x, y = merge_by_reads(sub, reads, target_reads)
        x = [i + total for i in x]

        ax.scatter(x=x, y=y, s=0.1)

        xticks.append(total + chromo_length / 2)
        xticklabels.append(contig.replace("chr", ""))

        total += chromo_length

    ax.set_xticks(xticks)
    ax.set_xticklabels(xticklabels)

    ax.set_ylim((0, 8))
    ax.set_xlim((0, total))
    plt.xlabel("Chromosome")
    plt.ylabel("Estimated copy number")
    plt.title(f"Illumina CNV plot: {sample}")

    plt.savefig(output, dpi=300)
    plt.close(fig)
