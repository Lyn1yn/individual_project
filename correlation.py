import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pyBigWig
from pathlib import Path
from adjustText import adjust_text
from scipy.ndimage import gaussian_filter1d



#1. Bacic settings
# input sample folder
samples_dir = Path("/gpfs01/share/BioinfMSc/Matt_Projects/samples")

# output folder
PROJECT_DIR = Path(__file__).resolve().parent.parent

outdir = PROJECT_DIR / "correlation"


outdir.mkdir(parents=True, exist_ok=True)

chromosomes = [f"chr{i}" for i in range(1, 23)] #only select autosomes
target_np_reads = 5000 #target reads in nanopore script




#2.merge nanopore bins
def merge_nanopore(values, bin_width, target_reads): #merge nanopore bins
    x_list = [] #create a list contains the position of x-axis, when creat a new plot bin, put a new coordinate in this list
    y_list = [] #create a list contains copy number

    group = [] #create a list to store the CNV copy numbers
    reads_sum = 0 #to record the estimated reads number accumulated for this group so far. Don't have any bins at first, so is 0
    start = 0 #to record this group started with which bin, to start with bin 0

    for i, cnv_value in enumerate(values): #read the CNV copy number of chromosomes one by one. The i here is the bins' number.
        reads_sum += cnv_value / 2 * 100 #to calculate estimate reads
        group.append(cnv_value) #put CNV copy numbers into the list created before

        if reads_sum >= target_reads: #to check whether the merged reads is enough
            x_list.append(((start + i + 1) / 2) * bin_width) #calculate the x positon for the points
            y_list.append(np.median(group)) #if the merged reads is enough, calculate the median copy numbers

            group = [] #clear the group for the next group of merged bins
            reads_sum = 0 #clear the reads sum, then calculate for a new merge group
            start = i + 1 #the new merge group start from the next bin

    if group: #after for loop, if the last group of reads does not reach target_reads, also keep data from the last group
        x_list.append(((start + len(values)) / 2) * bin_width) #get the x axis midpoint coordinate for the last group
        y_list.append(np.median(group)) #get the median CNV value

    return np.array(x_list), np.array(y_list) #get the x axis and copy number values for y axis position





#3.merge illumina bins
def merge_illumina(sub, target_reads):
    reads = sub["reads"].values #extract reads column as numpy array
    starts = sub["start"].values #extract start position column as numpy array
    cn = sub["copy_number"].values #extract copy number column as numpy array

    x_list = [] #create a list contains the position of x-axis, when creat a new plot bin, put a new coordinate in this list
    y_list = [] #create a list contains copy number

    reads_sum = 0 #record the total number of reads (No add bin, so start with 0)
    x_group = [] #record the x position of each bin in the current merge group
    cn_group = [] #record the copy number of each bin in the current merge group

    for i in range(len(reads)):
        reads_sum += reads[i] # add current bin's reads to cumulative sum
        x_group.append(starts[i]) #record current bin's start position
        cn_group.append(cn[i]) #record current bin's copy number

        if reads_sum >= target_reads: #to determine whether the cumulative reads is enough (if enough, output as a point)
            x_list.append(np.median(x_group)) #calculate median x-value of the bins as x-coordinate
            y_list.append(np.median(cn_group)) #calculate median copy number of the bins as y-coordinate

            reads_sum = 0 #clear cumulative reads, prepare for the next group merge
            x_group = [] #clear x position of current merge group
            cn_group = [] #clear copy number of current merge group

    if x_group:  #after for loop, if the last group of reads does not reach target_reads, also keep data from the last group
        x_list.append(np.median(x_group)) #calculate median x_value of the last group
        y_list.append(np.median(cn_group)) #calculate median copy number of the last group

    return np.array(x_list), np.array(y_list) #get the x axis and copy number values for y axis position



#3b. raw (no-merge) illumina points, used only for computing no-merge Illumina CN load colouring
def raw_illumina_points(sub):
    """Use original Illumina target bins directly, without merging."""
    if len(sub) == 0:
        return np.array([], dtype=float), np.array([], dtype=float)

    starts = sub["start"].values.astype(float)
    stops = sub["stop"].values.astype(float)
    x = (starts + stops) / 2
    y = sub["copy_number"].values.astype(float)
    return x, y



#4. read illumina input
def read_illumina_counts(sample_dir, sample): #read illumina input (get the gz and bw file in a same format)
    gz_file = sample_dir / "illumina" / f"{sample}.tumor.target.counts.gz" #input
    bw_file = sample_dir / "illumina" / f"{sample}.tumor.target.counts.bw" #input

    if gz_file.exists(): #for the gz files
        df = pd.read_csv(gz_file, sep="\t", comment="#") #read gz files into pandas dataframe. Seperate dataframe by tab, skip the comment line start with #.
        reads_col = df.columns[4] #take the 5th column as reads

        df = df[["contig", "start", "stop", reads_col]].copy() #only keep these 4 columns
        df = df.rename(columns={reads_col: "reads"}) #change the 5th column name into reads (different files probably have a different name)

        return df, "gz"

    if bw_file.exists(): #for the bw files
        bw = pyBigWig.open(str(bw_file)) #use pybigwig to open bw files
        chrom_sizes = bw.chroms() #read the chromosomes and the length of chromosomes

        records = [] #creat a list to store information from bw files

        for chrom in chromosomes: #read autosomes for the chromosomes
            for start, stop, value in bw.intervals(chrom, 0, chrom_sizes[chrom]): #read all intervals form 0 to the end of the chromosomes from bw files. Each intervals have 3 parts: start, stop and value(the reads)
               records.append({ #store current intervals into the records list
                    "contig": chrom, #the chromosome name
                    "start": start, #the start of interval
                    "stop": stop, #the end of interval
                    "reads": value #the interval value
                }) #close the record

        bw.close() #close bw files

        df = pd.DataFrame(records) #change record into pandas dataframe

        return df, "bw" #"bw" is a lable to tell the sample comes form a bw file




#5. prepare illumina dataframe
def prepare_illumina_df(illumina_df): #prepare illumina dataframe
    illumina_df = illumina_df[illumina_df["contig"].isin(chromosomes)].copy() #keep autosomes only

    median_reads = illumina_df["reads"].median() #calculate median read counts
    illumina_df["copy_number"] = 2 * illumina_df["reads"] / median_reads #calculate copy number
    
    illumina_df["chrom_num"] = illumina_df["contig"].str.replace("chr", "").astype(int) #convert chr1..chr22 to integers for sorting
    illumina_df = illumina_df.sort_values(["chrom_num", "start"]).copy() #sort by chromosome and position
    
    illumina_bin_size = (illumina_df["stop"] - illumina_df["start"]).median() ##caculate nanopore style target reads for samples

    return illumina_df, median_reads, illumina_bin_size




#6. Match illumina and nanopore points (by position)
def match_by_nearest_x(illumina_x, illumina_y, nanopore_x, nanopore_y, tolerance): #match illumina and nanopore points by nearest genomic position
    matched_illumina = [] #save CN value for successfuly matched illumina points
    matched_nanopore = [] #save CN vaalue for successfully matched nanopore points
    matched_distances = [] #store the distance between 2 points

    used_nanopore = set() #record which nanopore points have already been used, each nanopore point can only be matched once

    for i, ix in enumerate(illumina_x): #go through each illumina merged point
        distances = np.abs(nanopore_x - ix) #calculate distance between this illumina point and all nanopore points
        j = np.argmin(distances) #find the index of the nearest nanopore point

        if distances[j] <= tolerance and j not in used_nanopore: #only match if within tolerance and nanopore point not already used
            matched_illumina.append(illumina_y[i]) #record CN value for current illumina point
            matched_nanopore.append(nanopore_y[j]) #record CN value for current nanopore point
            matched_distances.append(distances[j]) #record distance between 2 points
            used_nanopore.add(j) #Mark this Nanopore point as used to avoid duplicate matching.

    return ( 
        np.array(matched_illumina),
        np.array(matched_nanopore),
        len(matched_illumina), #matched points number
        np.array(matched_distances)
    )




#7.main loop
summary = [] #creat a list, to store the correlation result for each samples

##prepare data
for sample_dir in sorted(samples_dir.iterdir()): #iterate through everything under the shared sample directory

    if not sample_dir.is_dir(): #skip non-folder (e.g. README.md)
        continue

    sample = sample_dir.name #get the sample names, can create pathway for the samples later
    nanopore_cnv_file = sample_dir / "nanopore" / "CNV.npy" #input nanopore copy number data
    nanopore_dict_file = sample_dir / "nanopore" / "CNV_dict.npy" #input nanopore original bin width data
    illumina_df, illumina_source = read_illumina_counts(sample_dir, sample) #get illumina input and get dataframe
    
    #prepare nanopore data
    cnv = np.load(nanopore_cnv_file, allow_pickle=True).item() #read cnv data for nanopore
    cnv_dict = np.load(nanopore_dict_file, allow_pickle=True).item() #read cnv_dict data for nanopore
    nanopore_bin_width = cnv_dict["bin_width"] #take bin width out
    

    # prepare Illumina data
    illumina_df, median_reads, illumina_bin_size = prepare_illumina_df(illumina_df) #calls the previously written functions to organize the ilumina table

    # apply Nanopore bin width to Illumina bin width sample to sample
    merged_bin_width = target_np_reads / 100 * nanopore_bin_width #calculate nanopore merged bin width, different for each sample
    illumina_target_reads = round(merged_bin_width / illumina_bin_size * median_reads) #convert nanopore merged bin width to Illumina reads

    all_illumina = [] #build list to store all matched copy number values
    all_nanopore = [] #build list to store all matched copy number values

    nanopore_total_estimated_reads = 0

    illumina_total_points = 0 #how many illumina have after merge
    nanopore_total_points = 0 #how may nanopore have after merge
    matched_total_points = 0 #the matched points number for correlation

    match_distance_values = [] #store genomic distances between matched Illumina and Nanopore points

    # per-chromosome matched CN values (needed for intra-chrom threshold)
    _illumina_per_chrom = {}
    _nanopore_per_chrom  = {}

    # per-chromosome ALL raw (no-merge) Illumina CN values (no Nanopore matching),
    # only used to colour plots by no-merge Illumina CN load
    _illumina_per_chrom_no_merge = {}


##process chromosome-by-chromosome
    #illumina and nanopore data for each chromosome
    for chrom in chromosomes: #from chromosome 1 to chromosome 22
        
        sub = illumina_df[illumina_df["contig"] == chrom].copy() #only take out the information current chromosome have

        # Nanopore reads-based merged points
        nanopore_values = np.array(cnv[chrom], dtype=float) #take out the copy number array from cnv for current chromosome

        nanopore_estimated_reads = nanopore_values / 2 * 100 #calculate reads for each bin
        nanopore_total_estimated_reads += np.sum(nanopore_estimated_reads) #add the reads of current chromosome into total reads of the sample, to get the approximate coverage.

        np_x, np_y = merge_nanopore( #merge bins, 5000reads per points. get the genomic position(x-axis) and copy number(y-axis) for merged points
            nanopore_values,
            nanopore_bin_width,
            target_np_reads
        )

        # Illumina reads-based merged points
        il_x, il_y = merge_illumina( #merge points, get the genomic position(x-axis) and copy number(y-axis) for merged points
            sub,
            illumina_target_reads
        )

        #smoothing
        sigma = 2  #use gaussian filter to smooth
        if len(np_y) > 1: #need at least 2 points to apply smoothing
            np_y = gaussian_filter1d(np_y, sigma=sigma) #smooth nanopore CN value to reduce noise
        if len(il_y) > 1: #need at least 2 points to apply smoothing
            il_y = gaussian_filter1d(il_y, sigma=sigma) #smooth illumina CN values to reduce noise


        #set pairing rule
        tolerance = merged_bin_width / 2 #if illumina and nanopore point are within half a merged bin width, consider them matched

        il_common, np_common, n, distances = match_by_nearest_x(
            il_x,
            il_y,
            np_x,
            np_y,
            tolerance
        )

        match_distance_values.extend(distances) #summarise the distances into a list

        all_illumina.extend(il_common) #add the matched copy number values into the list create previously (all_illumina list)
        all_nanopore.extend(np_common) #add the matched copy number values into the list create previously (all_nanopore list)

        # store per-chrom matched CN for intra-chrom threshold calculation
        _illumina_per_chrom[chrom] = np.array(il_common, dtype=float)
        _nanopore_per_chrom[chrom]  = np.array(np_common,  dtype=float)

        illumina_total_points += len(il_y) #calculate the total number of illumina matched points
        nanopore_total_points += len(np_y) #calculate the total number of nanopore matched points
        matched_total_points += n #calculate the matched points number

        # no-merge Illumina CN load: use ALL raw Illumina bins directly,
        # no Nanopore matching. Only used to colour plots.
        _, il_y_raw = raw_illumina_points(sub)
        _illumina_per_chrom_no_merge[chrom] = np.array(il_y_raw, dtype=float)


    all_illumina = np.array(all_illumina, dtype=float) #convert list to numpy array, and make sure it only contains numbers
    all_nanopore = np.array(all_nanopore, dtype=float) #convert list to numpy array, and make sure it only contains numbers



    #scatter plots

    #distance between matched points
    median_match_distance = np.nan #default to NaN in case no points were matched
    if len(match_distance_values) > 0:
        median_match_distance = np.median(match_distance_values) #median genomic distance between matched illumina and nanopore points

    #copy number standard deviation
    illumina_cn_sd = np.std(all_illumina) #calculate cn standard deviation, to see the fluctuation
    nanopore_cn_sd = np.std(all_nanopore) #calculate cn standard deviation, to see the fluctuation

    # platform-specific normalisation features
    illumina_median_cn = np.median(all_illumina) #calculate the median of CN value for illumina (baseline)
    nanopore_median_cn = np.median(all_nanopore) #calculate the meidan of CN value for nanopore (baseline)
    median_cn_difference = nanopore_median_cn - illumina_median_cn #show the difference of baselines between two platforms



    # Copy number load
    # It first labels each matched point as altered/normal, then only counts altered runs that contain at least min_consecutive_points in a row.

    normal_copy_number   = 2.0
    robust_std_floor     = 0.25
    abs_cn_threshold     = 0.60   # CN must be <=1.40 or >=2.60 to pass absolute CN criterion
    abs_z_threshold      = 2.00   # minimum z-score paired with abs threshold
    z_threshold          = 3.00

    # Require several consecutive altered matched points before counting them.
    min_consecutive_points = 3

    def compute_robust_std(values):
        """mad-based robust std, floored to avoid zero."""
        values = np.asarray(values, dtype=float) #convert input to a float numpy array
        values = values[np.isfinite(values)] #remove any NaN or Inf values
        if len(values) == 0: #if nothing is left after filtering, return the floor value
            return robust_std_floor
        mad = np.median(np.abs(values - np.median(values))) #compute mad (median of absolute deviations from the median)
        return max(1.4826 * mad, robust_std_floor) #scale mad to approxiate std, enforce minimum floor to avoid zero


    def is_point_altered(cn_value, robust_std, baseline=2.0):
        """Return True if one merged CN point passes threshold."""
        if not np.isfinite(cn_value): #skip NaN or Inf CN values
            return False

        normal_shift = abs(cn_value - baseline) #absolute deviation of this point from the baseline copy number (default 2)
        normal_z     = normal_shift / robust_std #how many units away from baseline this point is

        return (
            (normal_shift >= abs_cn_threshold and normal_z >= abs_z_threshold) #large absolute shift and moderately high z-score
            or normal_z >= z_threshold #or extreme z-score alone is enough
        )


    def count_consecutive_altered_points(flags, min_run):
        """
        Keep only altered points that are part of a consecutive altered run
        with length >= min_run. Isolated altered points are treated as noise.
        """
        flags = np.asarray(flags, dtype=bool) #ensure flags is a boolean numpy array
        keep = np.zeros(len(flags), dtype=bool) #create a numpy array contains everything set to false (initially make an assumptiojn that no points are retained, subsequently within the loop, only points meeting the condition are set to True, indicating that they are to be retained.)

        run_start = None #track where the run start
        for i, flag in enumerate(flags): #iterate through each point
            if flag and run_start is None: #run start
                run_start = i
            if (not flag or i == len(flags) - 1) and run_start is not None: #this is the last point so ends run
                run_end = i if not flag else i + 1 #when the last point is false, not include this point into run; if the last point is true, include it into run.
                if run_end - run_start >= min_run: #calculate the length of run, only the run meets minimum length keeps
                    keep[run_start:run_end] = True #change the points into 'true' in the keep numpy array.
                run_start = None #run end

        return keep #keep consecutive runs



    def compute_cn_load(all_cn, per_chrom_cn, baseline=2.0):
        """
        Conservative CN load:
        1. label each matched point using Tom's thresholds;
        2. keep only consecutive altered runs;
        3. return percentage of matched points retained as altered.
        robust_std here is GLOBAL (computed once across all chromosomes' matched
        points), not per-chromosome.
        """
        robust_std = compute_robust_std(all_cn) #single global robust std, used for every chromosome

        kept_flags_all = [] #prepare an empty list to keep final result(after consecutive filter)
        point_flags_all = [] #prepare an empty list to keep every points'condition(before consecutive filter)

        for chrom, cn_arr in per_chrom_cn.items(): #process each chrmosome seperatelty
            cn_arr = np.asarray(cn_arr, dtype=float) #ensure cn values are float numpy array
            cn_arr = cn_arr[np.isfinite(cn_arr)] #exlude any NaN or Inf values
            if len(cn_arr) == 0: #skip chrmosomes with no valid points
                continue

            point_flags = np.array([
                is_point_altered(v, robust_std, baseline=baseline)
                for v in cn_arr
            ], dtype=bool) #for every cn value in the chromosome, use "is_point_altered" to define if it is altered or normal(the results are true and false for each points)

            kept_flags = count_consecutive_altered_points(
                point_flags,
                min_consecutive_points
            ) #give the true/false results to the "count_consecutive_altered_points", only keep runs have ≥ 3 alters. Change the true points but not meet requirments into false.

            point_flags_all.extend(point_flags) #keep true/false result before consecutive filter into list
            kept_flags_all.extend(kept_flags) #keep true/false result after filter into list.




        consecutive_load = np.mean(kept_flags_all) * 100 if kept_flags_all else 0.0 #get the true percentage after filter

        point_load = np.mean(point_flags_all) * 100 if point_flags_all else 0.0 #get the true percentage before filter
        #by comparing these two percentage can see how many points are discarded by filtering

        return consecutive_load, point_load, robust_std



    # Build per-chromosome arrays for matched points.
    # need to re-collect them from the per-chrom loop, so attach them here.
    illumina_cn_load, _, _ = compute_cn_load(all_illumina, _illumina_per_chrom)
    nanopore_cn_load, _, _ = compute_cn_load(all_nanopore, _nanopore_per_chrom)

    # Illumina CN load computed from raw (no-merge) bins, used only for colouring plots below.
    all_illumina_no_merge = np.concatenate(list(_illumina_per_chrom_no_merge.values())) if _illumina_per_chrom_no_merge else np.array([], dtype=float)
    illumina_cn_load_no_merge, _, _ = compute_cn_load(all_illumina_no_merge, _illumina_per_chrom_no_merge)

    # Ploidy-aware CN load for samples estimated as ploidy 3 by ichorCNA.
    # The CN profile is normalised so its median maps to 2, but for a ploidy-3 tumour
    # the true baseline is 3 copies. We rescale relative CN (centred on 2) into an
    # approximate absolute CN (centred on ichorCNA ploidy) and measure alteration
    # relative to that ploidy baseline instead of 2.
    ploidy_override = {
        "STG05R2-13_b-E01": 3.0
    }

    ichor_ploidy = ploidy_override.get(sample, 2.0)
    is_ploidy_adjusted_sample = sample in ploidy_override

    illumina_cn_load_ploidy_adjusted = np.nan
    nanopore_cn_load_ploidy_adjusted = np.nan

    if is_ploidy_adjusted_sample:
        # Convert relative CN scaled around 2 into approximate absolute CN scaled around ichorCNA ploidy
        illumina_abs_like = all_illumina * ichor_ploidy / 2.0
        nanopore_abs_like = all_nanopore * ichor_ploidy / 2.0

        illumina_per_chrom_abs_like = {
            chrom: arr * ichor_ploidy / 2.0
            for chrom, arr in _illumina_per_chrom.items()
        }

        nanopore_per_chrom_abs_like = {
            chrom: arr * ichor_ploidy / 2.0
            for chrom, arr in _nanopore_per_chrom.items()
        }

        illumina_cn_load_ploidy_adjusted, _, _ = compute_cn_load(
            illumina_abs_like,
            illumina_per_chrom_abs_like,
            baseline=ichor_ploidy
        )

        nanopore_cn_load_ploidy_adjusted, _, _ = compute_cn_load(
            nanopore_abs_like,
            nanopore_per_chrom_abs_like,
            baseline=ichor_ploidy
        )
    


    #difference between platforms
    cn_load_difference = nanopore_cn_load - illumina_cn_load
    cn_load_abs_difference = abs(cn_load_difference)


    #calculate pearson
    pearson_r = pd.Series(all_illumina).corr( #change CN value into pandas Series. Calculate correlation.
        pd.Series(all_nanopore), #change CN value into pandas Series
        method="pearson"
    )



    summary.append({
        "sample": sample,
        "illumina_source": illumina_source,

        "illumina_points": illumina_total_points,
        "nanopore_points": nanopore_total_points,
        "matched_points": matched_total_points,

        "unused_illumina_points": illumina_total_points - matched_total_points,
        "unused_nanopore_points": nanopore_total_points - matched_total_points,

        "illumina_cn_sd": illumina_cn_sd,
        "nanopore_cn_sd": nanopore_cn_sd,

        "median_match_distance": median_match_distance,

        # platform-specific normalisation
        "illumina_median_cn": illumina_median_cn,
        "nanopore_median_cn": nanopore_median_cn,
        "median_cn_difference": median_cn_difference,

        #nanopore approximate coverage(total estimated reads number)
        "nanopore_total_estimated_reads": nanopore_total_estimated_reads,
        
        #cn load
        "illumina_cn_load": illumina_cn_load,
        "nanopore_cn_load": nanopore_cn_load,
        "cn_load_difference": cn_load_difference,
        "cn_load_abs_difference": cn_load_abs_difference,

        # ploidy-adjusted CN load (only non-NaN for ichorCNA ploidy-3 samples)
        "ichor_ploidy_for_cn_load": ichor_ploidy,
        "is_ploidy_adjusted_sample": is_ploidy_adjusted_sample,
        "illumina_cn_load_ploidy_adjusted": illumina_cn_load_ploidy_adjusted,
        "nanopore_cn_load_ploidy_adjusted": nanopore_cn_load_ploidy_adjusted,

        # Illumina CN load from raw (no-merge) bins, used for colouring plots
        "illumina_cn_load_no_merge": illumina_cn_load_no_merge,


        "pearson_r": pearson_r,
    })




#8.add sample names next to points
def add_sample_labels(x_values, y_values, sample_names):
    texts = []
    for x, y, sample in zip(x_values, y_values, sample_names):
        texts.append(plt.text(x, y, sample, fontsize=6, alpha=0.8))
    adjust_text(texts, arrowprops=dict(arrowstyle="-", color="gray", lw=0.5))



#9. plot
summary_df = pd.DataFrame(summary) #create a dataframe to store summary list from each loop
summary_df.to_csv(outdir / "correlation_summary.csv", index=False) #save the dataframe into csv file
summary_df = summary_df.sort_values("sample") #sort samples by their samples' name
summary_df["sample_id"] = summary_df["sample"].str.extract(r'-(\d+)_')[0] #get the sample number, and save as sample_id
#sort by sample number
summary_df = summary_df.sort_values(
    "sample_id",
    key=lambda s: pd.to_numeric(s, errors="coerce"),
    kind="stable"
).reset_index(drop=True)



summary_df.to_csv(
    outdir / "reads_based_correlation_summary.csv",
    index=False
)

print("\nPearson summary:")
print(summary_df["pearson_r"].describe())


#cn load and pearson relationship results
print("\nRelationship between copy number load and Pearson:")
print(
    "Illumina CN load vs Pearson:",
    summary_df["illumina_cn_load"].corr(summary_df["pearson_r"], method="pearson")
)
print(
    "Nanopore CN load vs Pearson:",
    summary_df["nanopore_cn_load"].corr(summary_df["pearson_r"], method="pearson")
)
print(
    "Absolute CN load difference vs Pearson:",
    summary_df["cn_load_abs_difference"].corr(summary_df["pearson_r"], method="pearson")
)




x = np.arange(len(summary_df)) #create samples' x-axis position


#pearson plot
plt.figure(figsize=(12, 5))

plt.bar( #draw a bar at each position on x axis, the height of bar is the pearson value of each sample
    x,
    summary_df["pearson_r"]
)

plt.xticks(x, summary_df["sample"], rotation=90) #set x-axis labels
plt.ylim(-1, 1) #set y-axis range

plt.ylabel("Pearson correlation coefficient")
plt.xlabel("Sample")
plt.title("Pearson correlation between Illumina and Nanopore CNV profiles")
plt.axhline(0, linestyle="--", linewidth=0.8)

plt.tight_layout() #automatic layout adjustment
plt.savefig(outdir / "pearson_correlation_by_sample.png", dpi=300)
plt.close()





# Illumina copy number load vs Pearson
plt.figure(figsize=(6, 5))

sc = plt.scatter(
    summary_df["illumina_cn_load"],
    summary_df["pearson_r"],
    c=summary_df["illumina_cn_load_no_merge"],
    cmap="viridis",
    s=45,
    edgecolors="black",
    linewidths=0.3
)
add_sample_labels(summary_df["illumina_cn_load"], summary_df["pearson_r"], summary_df["sample_id"])

special = summary_df[summary_df["is_ploidy_adjusted_sample"]].copy()

plt.scatter(
    special["illumina_cn_load_ploidy_adjusted"],
    special["pearson_r"],
    marker="*",
    s=180,
    facecolors="none",
    edgecolors="red",
    linewidths=1.2,
    label="Ploidy-adjusted CN load (ichorCNA ploidy=3)"
)

for _, row in special.iterrows():
    plt.annotate(
        "",
        xy=(row["illumina_cn_load_ploidy_adjusted"], row["pearson_r"]),
        xytext=(row["illumina_cn_load"], row["pearson_r"]),
        arrowprops=dict(arrowstyle="->", color="red", lw=0.8)
    )

plt.legend(fontsize=7)

cbar = plt.colorbar(sc)
cbar.set_label("Illumina CN load from no-merge data (%)")

plt.xlabel("Illumina copy number load (%)")
plt.ylabel("Pearson correlation coefficient")
plt.title("Illumina copy number load vs Pearson correlation")
plt.axhline(0, linestyle="--", linewidth=0.8)

plt.tight_layout()
plt.savefig(outdir / "illumina_cn_load_vs_pearson.png", dpi=300)
plt.close()



# Nanopore copy number load vs Pearson
plt.figure(figsize=(6, 5))

sc = plt.scatter(
    summary_df["nanopore_cn_load"],
    summary_df["pearson_r"],
    c=summary_df["illumina_cn_load_no_merge"],
    cmap="viridis",
    s=45,
    edgecolors="black",
    linewidths=0.3
)
add_sample_labels(summary_df["nanopore_cn_load"], summary_df["pearson_r"], summary_df["sample_id"])

special = summary_df[summary_df["is_ploidy_adjusted_sample"]].copy()

plt.scatter(
    special["nanopore_cn_load_ploidy_adjusted"],
    special["pearson_r"],
    marker="*",
    s=180,
    facecolors="none",
    edgecolors="red",
    linewidths=1.2,
    label="Ploidy-adjusted CN load (ichorCNA ploidy=3)"
)

for _, row in special.iterrows():
    plt.annotate(
        "",
        xy=(row["nanopore_cn_load_ploidy_adjusted"], row["pearson_r"]),
        xytext=(row["nanopore_cn_load"], row["pearson_r"]),
        arrowprops=dict(arrowstyle="->", color="red", lw=0.8)
    )

plt.legend(fontsize=7)

cbar = plt.colorbar(sc)
cbar.set_label("Illumina CN load from no-merge data (%)")

plt.xlabel("Nanopore copy number load (%)")
plt.ylabel("Pearson correlation coefficient")
plt.title("Nanopore copy number load vs Pearson correlation")
plt.axhline(0, linestyle="--", linewidth=0.8)

plt.tight_layout()
plt.savefig(outdir / "nanopore_cn_load_vs_pearson.png", dpi=300)
plt.close()




