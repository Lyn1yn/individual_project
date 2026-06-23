import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path


#input
samples_dir = Path("/gpfs01/share/BioinfMSc/Matt_Projects/samples")
#bin width get from nanopore data
bin_width_file = "/gpfs01/home/mbxll1/CNS_cancer_project/scatter_plot_illumina/all_bin_width.txt"
#output
outdir = Path("/gpfs01/home/mbxll1/CNS_cancer_project/scatter_plot_illumina/5000")
outdir.mkdir(parents=True, exist_ok=True)


# use the same target reads as nanopore script
target_np_reads = 5000

bin_width_df = pd.read_csv(  #read all_bin_width.txt from nanopore data, to creat dataframe
    bin_width_file,
    sep=r"\s+",  #read file using spaces as delimiters
    header=None, #this txt file do not have a header
    names=["sample", "bin_width",]) #name the 2 columns




def merge_by_reads(sub, reads, target_reads): #define a function contains 3 inputs
    x_list = [] #create a list contains the position of x-axis, when creat a new plot bin, put a new coordinate in this list
    y_list = [] #create a list contains copy number

    reads_sum = 0 #record the total number of reads (No add bin, so start with 0)
    x_group = [] #record the x position of each bin in the current merge group
    cn_group = [] #record the copy number of each bin in the current merge group


    for index, row in sub.iterrows(): #read the data row by row, index is the row number

        reads_sum += row[reads] #add read of current bin into cumulative reads (until greater than target_reads)
        x_group.append(row["start"]) #put the start position of current bin into x_group
        cn_group.append(row["copy_number"]) #put the current bin copy number into cn_group

        if reads_sum >= target_reads: #to determine whether the cumulative reads is enough (if enough, output as a point)
            x_list.append(pd.Series(x_group).median()) #calculate median x-value of the bins as x-coordinate
            y_list.append(pd.Series(cn_group).median()) #calculate median copy number of the bins as y-coordinate.

            reads_sum = 0 #clear cumulative reads, prepare for the next group merge
            x_group = [] #clear x position of current merge group
            cn_group = [] #clear copy number of current merge group


    if x_group: #after for loop, if the last group of reads does not reach target_reads, also keep data from the last group
        x_list.append(pd.Series(x_group).median()) #calculate median x_value of the last group
        y_list.append(pd.Series(cn_group).median()) #calculate median copy number of the last group
    return x_list, y_list




#loop through all samples
for sample, bin_width in zip (bin_width_df["sample"], bin_width_df["bin_width"]):  # bind the sample name and its corresponding bin_width together in a loop.

    #input illumina target counts file
    file = f"{samples_dir}/{sample}/illumina/{sample}.tumor.target.counts.gz"


    ## skip if illumina file does not exist (because there are bw files)
    if not Path(file).exists():
        print(f"{sample}: illumina file not found, skipping")
        continue

    df = pd.read_csv(file, sep="\t", comment="#")


    #read illumina target counts
    df = pd.read_csv(file, sep="\t", comment="#")

    ##the 5th column is read counts
    reads = df.columns[4]

    
    #only select autosomes
    autosomes = [f"chr{i}" for i in range(1, 23)]
    df = df[df["contig"].isin(autosomes)].copy()   #only keep contig column in autosomes


    #prepare y axis
    #calculate median read counts
    median_reads = df[reads].median()
    
    #calculate copy number
    df["copy_number"] = 2 * df[reads] / median_reads

    
    #convert chr1..chr22 to integers for sorting
    df["chrom_num"] = df["contig"].str.replace("chr", "").astype(int)

    #sort by chromosome and position
    df = df.sort_values(["chrom_num", "start"]).copy()

    #prepare x axis


    ##caculate nanopore style target reads for samples
    illumina_bin_size = (df["stop"] - df["start"]).median()

    # calculate nanopore merged bin width, different for each sample
    merged_bin_width = target_np_reads / 100 * bin_width

    # convert nanopore merged bin width to Illumina reads, median reads is reads per bin for illumina
    target_reads = merged_bin_width / illumina_bin_size * median_reads
    target_reads = round(target_reads)



    ### estimate how many Illumina original bins are merged into one point
    approx_illumina_bins_per_point = target_reads / median_reads

    # estimate Illumina merged bin width
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





    #output picture name
    output = outdir / f"{sample}_illumina_merged_reads.png"

    #create plot
    fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(18,4))

    total = 0  # to draw chromosome1 from origin point
    xticks = []  # create a list for chromosome labels
    xticklabels = []  # create a list for chromosome labels

    ##plot each chromosome
    for contig in autosomes:
        sub = df[df["contig"] == contig].copy()
        chromo_length = sub["stop"].max()

        #merge neibouring bins by accumulated reads
        x,y = merge_by_reads(sub, reads, target_reads)
        #add total chromosomes one after another
        x = [i + total for i in x]

        #plot
        ax.scatter(
            x=x,
            y=y,
            s=0.1)

        xticks.append(total + chromo_length / 2)
        xticklabels.append(contig.replace("chr",""))

        total += chromo_length


    #label
    ax.set_xticks(xticks)
    ax.set_xticklabels(xticklabels)

    ax.set_ylim((0, 8))
    ax.set_xlim((0, total))
    plt.xlabel("Chromosome")
    plt.ylabel("Estimated copy number")
    plt.title(f"Illumina CNV plot: {sample}")


    plt.savefig(output, dpi=300)
    plt.close(fig) #close picture to free up memory
