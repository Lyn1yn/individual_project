import pyBigWig
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path

sample = "STG05-52_c-E03"

# input: hared data directory
bw_file = Path(
    f"/gpfs01/share/BioinfMSc/Matt_Projects/samples/{sample}/illumina/{sample}.tumor.target.counts.bw"
)

# Input bin width file and output directory - derived relative to this script's location
# Assumes this script lives at: <project_root>/illumina_profile/52.py
# bin_width_file is expected at: <project_root>/scatter_plot_nanopore/all_bin_width.txt
# Output goes to:                <project_root>/scatter_plot_illumina/5000/
SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_DIR = SCRIPT_DIR.parent
bin_width_file = PROJECT_DIR / "scatter_plot_nanopore" / "all_bin_width.txt"
outdir = PROJECT_DIR / "scatter_plot_illumina" / "5000"
outdir.mkdir(parents=True, exist_ok=True)

target_np_reads = 5000 #nanopore each point contains 5000reads, to infer illumina's expected bin width


def merge_by_reads(sub, reads_col, target_reads):
    x_list = [] #genome position
    y_list = [] #copy number

    reads_sum = 0 #record the accumulation of reads
    x_group = [] #start position for the current group
    cn_group = [] #copy number for current group

    for _, row in sub.iterrows(): #iterate row by row
        reads_sum += row[reads_col] #add this bin's read count
        x_group.append(row["start"]) #record bin's start position
        cn_group.append(row["copy_number"]) #record bin's copy number

        if reads_sum >= target_reads: #when threshold reached, close the group and emit one point
            x_list.append(pd.Series(x_group).median()) #median of the group's start positions becomes the point's x
            y_list.append(pd.Series(cn_group).median()) #median copy number become y

            reads_sum = 0 #reset all 3 accumulators and start the next group
            x_group = []
            cn_group = []

    if x_group: #leftover tail after the loop that never reached the target
        x_list.append(pd.Series(x_group).median())
        y_list.append(pd.Series(cn_group).median())

    return x_list, y_list #2 lists


# Read nanopore bin width for this sample
bin_width_df = pd.read_csv(
    bin_width_file,
    sep=r"\s+",
    header=None,
    names=["sample", "bin_width"],
)

bin_width_row = bin_width_df[bin_width_df["sample"] == sample] #boolean index to select the row of current sample

if bin_width_row.empty:
    raise ValueError(f"{sample} not found in all_bin_width.txt") #if nothing found, the sample is missing

bin_width = float(bin_width_row["bin_width"].iloc[0]) #take out bin width number


# Read bigWig
if not bw_file.exists(): #check the file is present or not
    raise FileNotFoundError(f"BigWig file not found: {bw_file}")

bw = pyBigWig.open(str(bw_file)) #open the bigwig file

records = [] #build dataframe

chroms = bw.chroms() #get the chromosomes

if "chr1" in chroms:
    autosomes = [f"chr{i}" for i in range(1, 23)]
else:
    autosomes = [str(i) for i in range(1, 23)]

for contig in autosomes: #process one chromosome
    if contig not in chroms:
        print(f"{contig} not found in bw, skipping")
        continue

    intervals = bw.intervals(contig) #get the position and reads number

    if intervals is None:
        print(f"{contig}: no intervals")
        continue

    for start, stop, value in intervals: #skip missing value
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

df = pd.DataFrame(records) #change list into a table

if df.empty:
    raise ValueError("No data extracted from BigWig file.")

df = df[df["contig"].isin(autosomes)].copy() #filter autosome again

reads = "reads" #store the column name
median_reads = df[reads].median() #get the median reads per bin

df["copy_number"] = 2 * df[reads] / median_reads #calculate copy number

df["chrom_num"] = df["contig"].str.replace("chr", "", regex=False).astype(int) #order chromosome
df = df.sort_values(["chrom_num", "start"]).copy()

illumina_bin_size = (df["stop"] - df["start"]).median() #calculate bin size

merged_bin_width = target_np_reads / 100 * bin_width #merged bin width for each sample

target_reads = merged_bin_width / illumina_bin_size * median_reads #calculate how many reads need to be merged
target_reads = round(target_reads) #round to an integer for use

approx_illumina_bins_per_point = target_reads / median_reads #sanity check, how many bins required for each point
illumina_approx_merged_bin_width = approx_illumina_bins_per_point * illumina_bin_size #sanity check, what size is the bin width, should be close to merged_bin_width

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






#plot
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

print(f"Saved: {output}")
