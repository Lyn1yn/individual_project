import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

# input sample folder
samples_dir = Path("/gpfs01/share/BioinfMSc/Matt_Projects/samples")

# output folder
outdir = Path("/gpfs01/home/mbxll1/CNS_cancer_project/scatter_plot_nanopore/5000")
outdir.mkdir(parents=True, exist_ok=True)

# autosomes only
chromosomes = [f"chr{i}" for i in range(1, 23)]

def bin_by_estimate_reads(values, bin_width, target_reads): #define a function contains 3 inputs: 

    x_list = [] #create a list contains the position of x-axis, when creat a new plot bin, put a new coordinate in this list
    y_list = [] #create a list contains copy number

    group = [] #create a list to store the CNV copy numbers
    reads_sum = 0 #to record the estimated reads number accumulated for this group so far. Don't have any bins at first, so is 0
    start = 0 #to record this group started with which bin, to start with bin 0

    for i, cnv_value in enumerate(values): #read the CNV copy number of chromosomes one by one. The i here is the bins' number.
        estimate_reads = cnv_value / 2 *100 #to calculate estimate reads
        group.append(cnv_value) #put CNV copy numbers into the list created before
        reads_sum += estimate_reads #put the bin reads into the merged reads
        if reads_sum >= target_reads: #to check whether the merged reads is enough
            y_list.append(np.median(group)) #if the merged reads is enough, calculate the median copy numbers
            x_list.append(((start +i + 1) / 2) * bin_width) #calculate the x positon for the points
            group = [] #clear the group for the next group of merged bins
            reads_sum = 0 #clear the reads sum, then calculate for a new merge group
            start = i + 1 #the new merge group start from the next bin
    if group: ## If group is not empty after the loop ends, it means there are still unprocessed bins remaining
        x_list.append(((start + len(values)) / 2) * bin_width) #get the x axis midpoint coordinate for the last group
        y_list.append(np.median(group)) #get the median CNV value
    return np.array(x_list), np.array(y_list) #get the x axis and copy number values for y axis position


#### loop through all sample folders
for sample_dir in sorted(samples_dir.iterdir()): #Iterate through everything in samples_dir

    if not sample_dir.is_dir():
        continue #Skip this if it's not a folder (exclude readme file)

    sample = sample_dir.name #give sample name, to name the images later

    cnv_file = sample_dir / "nanopore" / "CNV.npy" #automatically creat input file path for each sample
    cnv_dict_file = sample_dir / "nanopore" / "CNV_dict.npy" #automatically creat input file path for each sample

    output = outdir / f"{sample}_nanopore_merged_reads.png" #set name for each picture

    # import result data from cnv_from_bam
    cnv = np.load(cnv_file, allow_pickle=True).item()
    cnv_dict = np.load(cnv_dict_file, allow_pickle=True).item()

    # get bin width
    bin_width = cnv_dict["bin_width"]

    # target reads 5000 as reference
    target_reads = 5000 #reads of sample 5000


    #calculate merged bin width
    merged_bin_width = target_reads / 100 * bin_width
    print(
        sample,
        "original_bin_width:", bin_width,
        "target_reads:", target_reads,
        "merged_bin_width:", round(merged_bin_width)
    )


    # creat a picture
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(18, 4))

    # to draw chromosome1 from origin point
    total = 0

    # create 2 list for chromosome labels
    xticks = []
    xticklabels = []

    # plot
    for contig in chromosomes:
        values = np.array(cnv[contig], dtype=float) # only select first 22 chromosomes
        chromo_length = len(values) * bin_width

        # merge bins, 500reads per point
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

    fig.savefig(output, dpi=300) #save picture
    plt.close(fig) #close picture to free up memory
