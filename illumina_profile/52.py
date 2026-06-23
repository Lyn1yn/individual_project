import pyBigWig
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path


sample = "STG05-52_c-E03"

bw_file = Path(
    f"/gpfs01/share/BioinfMSc/Matt_Projects/samples/{sample}/illumina/{sample}.tumor.target.counts.bw"
)

bin_width_file = Path(
    "/gpfs01/home/mbxll1/CNS_cancer_project/scatter_plot_illumina/all_bin_width.txt"
)

outdir = Path(
    "/gpfs01/home/mbxll1/CNS_cancer_project/scatter_plot_illumina/4_fixed_target_reads/5000"
)
outdir.mkdir(parents=True, exist_ok=True)

target_np_reads = 5000


def merge_by_reads(sub, reads_col, target_reads):
    x_list = []
    y_list = []

    reads_sum = 0
    x_group = []
    cn_group = []

    for _, row in sub.iterrows():
        reads_sum += row[reads_col]
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


# 读取这个 sample 对应的 nanopore bin width
bin_width_df = pd.read_csv(
    bin_width_file,
    sep=r"\s+",
    header=None,
    names=["sample", "bin_width"],
)

bin_width_row = bin_width_df[bin_width_df["sample"] == sample]

if bin_width_row.empty:
    raise ValueError(f"{sample} not found in all_bin_width.txt")

bin_width = float(bin_width_row["bin_width"].iloc[0])


# 读取 bigWig
if not bw_file.exists():
    raise FileNotFoundError(f"BigWig file not found: {bw_file}")

bw = pyBigWig.open(str(bw_file))

records = []

# 判断 chromosome 名字是 chr1 还是 1
chroms = bw.chroms()

if "chr1" in chroms:
    autosomes = [f"chr{i}" for i in range(1, 23)]
else:
    autosomes = [str(i) for i in range(1, 23)]

for contig in autosomes:
    if contig not in chroms:
        print(f"{contig} not found in bw, skipping")
        continue

    intervals = bw.intervals(contig)

    if intervals is None:
        print(f"{contig}: no intervals")
        continue

    for start, stop, value in intervals:
        if value is None:
            continue

        records.append(
            {
                "contig": contig,
                "start": start,
                "stop": stop,
                "reads": value,
            }
        )

bw.close()

df = pd.DataFrame(records)

if df.empty:
    raise ValueError("No data extracted from BigWig file.")


# 只保留 autosomes
df = df[df["contig"].isin(autosomes)].copy()

# 计算 copy number
reads = "reads"
median_reads = df[reads].median()

df["copy_number"] = 2 * df[reads] / median_reads

# 染色体排序
df["chrom_num"] = df["contig"].str.replace("chr", "", regex=False).astype(int)
df = df.sort_values(["chrom_num", "start"]).copy()


# 计算 illumina target reads
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


# 输出文件名
output = outdir / f"{sample}_illumina_merged_reads.png"

fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(18, 4))

total = 0
xticks = []
xticklabels = []

for contig in autosomes:
    sub = df[df["contig"] == contig].copy()

    if sub.empty:
        continue

    chromo_length = sub["stop"].max()

    x, y = merge_by_reads(sub, reads, target_reads)

    x = [i + total for i in x]

    ax.scatter(
        x=x,
        y=y,
        s=0.1,
    )

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

print(f"Saved: {output}")
