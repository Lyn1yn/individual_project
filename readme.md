# 1. Individual Project

This repository contains four related workflows:

1. Nanopore copy-number profile generation
2. Illumina copy-number profile generation
3. Illumina–Nanopore correlation and CN-load analysis
4. ichorCNA analysis

# 2. Project structure

```text
individual_project/
├── README.md
├── scripts/
│   ├── 52.py
│   ├── check_all_bin_width.py
│   ├── correlation.py
│   ├── correlation.sh
│   ├── illumina_check_5000.py
│   ├── make_gc_hg38_500kb.R
│   ├── nanopore_5000.py
│   └── run_batch_ichorCNA_500kb.sh
├── environment/
│   ├── ichorcna.yml
│   └── ichorcna_R.yml
└── .gitignore
```

```text
The following directories are generated during analysis on the HPC:

individual_project/
├── scatter_plot_nanopore/
├── scatter_plot_illumina/
├── correlation/
└── ichorcna_autosome/
    ├── 1_primary_mapped/
    ├── 2_remove_RG/
    ├── 3_wig/
    ├── 4_ichorCNA_results/
    ├── reference/
    ├── software/
    └── logs/
```

# 3. Path setup

Run the following commands from the project root:

```bash
PROJECT_DIR="$(pwd)"
SCRIPTS_DIR="${PROJECT_DIR}/scripts"

NANOPORE_PLOT_DIR="${PROJECT_DIR}/scatter_plot_nanopore"
ILLUMINA_PLOT_DIR="${PROJECT_DIR}/scatter_plot_illumina"
CORRELATION_DIR="${PROJECT_DIR}/correlation"

ICHOR_DIR="${PROJECT_DIR}/ichorcna_autosome"
ICHOR_SOFTWARE_DIR="${ICHOR_DIR}/software"
ICHOR_REFERENCE_DIR="${ICHOR_DIR}/reference"
ICHOR_LOG_DIR="${ICHOR_DIR}/logs"
```

# 4. External raw data

The raw sequencing data were provided by Professor Matthew Loose and are stored in the shared HPC directory:

```text
/gpfs01/share/BioinfMSc/Matt_Projects
```

# 5. Script overview

| Script | Purpose |
|---|---|
| `check_all_bin_width.py` | Reads each sample's `CNV_dict.npy`, extracts the original Nanopore bin width, and writes the summary to `scatter_plot_nanopore/all_bin_width.txt`. |
| `nanopore_5000.py` | Processes Nanopore copy-number data and produces profiles merged to approximately 5000 estimated reads per point. |
| `illumina_check_5000.py` | Processes Illumina `.tumor.target.counts.gz` files and produces profiles matched to the Nanopore 5000-read setting. |
| `52.py` | Separately processes sample `STG05-52_c-E03`, whose Illumina data is only available as a BigWig (`.bw`) file rather than `.tumor.target.counts.gz`. |
| `correlation.py` | Matches Illumina and Nanopore copy-number points, calculates Pearson correlation and CN-load metrics, and generates summary tables and plots. |
| `correlation.sh` | Submits or runs the correlation analysis on the HPC system. |
| `make_gc_hg38_500kb.R` | Generates the hg38 500 kb GC-content WIG reference file. |
| `run_batch_ichorCNA_500kb.sh` | Processes Nanopore BAM files, generates 500 kb read-count WIG files, and runs ichorCNA. |

# 6. Setup

## 6.1 Conda environments

The batch workflow uses two Conda environments. You can either create them manually with the exact package versions below, or install them directly from the provided `environment/` YAML files.

### BAM and WIG processing environment

**Option A — manual install:**

```bash
conda create -n ichorcna \
  -c conda-forge \
  -c bioconda \
  samtools=1.23.1 \
  hmmcopy=0.1.1 \
  -y
```

**Option B — from YAML:**

```bash
conda env create -f environment/ichorcna.yml
```

Check the installation:

```bash
conda activate ichorcna
which samtools
which readCounter
samtools --version
readCounter --help
```


### R and ichorCNA environment

**Option A — manual install:**

```bash
conda create -n ichorcna_R \
  -c conda-forge \
  -c bioconda \
  r-base=4.4.3 \
  r-ichorcna=0.5.1 \
  bioconductor-hmmcopy \
  bioconductor-bsgenome.hsapiens.ucsc.hg38 \
  bioconductor-biostrings \
  r-optparse \
  r-stringr \
  r-readr \
  r-dplyr \
  -y
```

**Option B — from YAML:**

```bash
conda env create -f environment/ichorcna_R.yml
```

The installed ichorCNA R package v0.5.1 was used to access the bundled hg38 mappability and centromere reference files. The analysis itself was run using the GavinHaLab ichorCNA v0.4.0 source repository.

Check the R dependencies:

```bash
conda activate ichorcna_R

Rscript -e '
library(HMMcopy)
library(GenomicRanges)
library(GenomeInfoDb)
library(BSgenome.Hsapiens.UCSC.hg38)
library(Biostrings)
library(data.table)
library(optparse)
cat("R dependencies loaded successfully\n")
'
```

## 6.2 Install ichorCNA v0.4.0

Run these commands from the project root

```bash
PROJECT_DIR="$(pwd)"
ICHOR_DIR="${PROJECT_DIR}/ichorcna_autosome"
ICHOR_SOFTWARE_DIR="${ICHOR_DIR}/software"

mkdir -p "${ICHOR_SOFTWARE_DIR}"
cd "${ICHOR_SOFTWARE_DIR}"
```

Clone the selected version:

```bash
git clone \
  --branch v0.4.0 \
  --depth 1 \
  https://github.com/GavinHaLab/ichorCNA.git \
  GavinHaLab_ichorCNA_v0.4.0
```

Install it into the `ichorcna_R` environment:

```bash
conda activate ichorcna_R

ICHOR_REPO="${ICHOR_SOFTWARE_DIR}/GavinHaLab_ichorCNA_v0.4.0"

R CMD INSTALL "${ICHOR_REPO}"
```

Check the installation:

```bash
Rscript -e '
library(ichorCNA)
cat("ichorCNA version:", as.character(packageVersion("ichorCNA")), "\n")
cat("extdata:", system.file("extdata", package="ichorCNA"), "\n")
'
```

## 6.3 Generate the GC reference

Run these commands from the project root :

```bash
PROJECT_DIR="$(pwd)"
SCRIPTS_DIR="${PROJECT_DIR}/scripts"

conda activate ichorcna_R
cd "${SCRIPTS_DIR}"

Rscript make_gc_hg38_500kb.R
```

The expected output is:

```text
ichorcna_autosome/reference/g_hg38_500kb.wig
```

Check the WIG header:

```bash
PROJECT_DIR="$(pwd)"
GCWIG="${PROJECT_DIR}/ichorcna_autosome/reference/g_hg38_500kb.wig"

head "${GCWIG}"
```

The header should contain:

```text
fixedStep chrom=chr1 start=1 step=500000 span=500000
```

If the file contains `5e+05`, replace it with:

```bash
sed -i \
  's/step=5e+05/step=500000/g; s/span=5e+05/span=500000/g' \
  "${GCWIG}"
```

## 6.4 Reference files

After activating `ichorcna_R`, locate the package reference directory through the active environment:

```bash
conda activate ichorcna_R

EXTDATA="${CONDA_PREFIX}/lib/R/library/ichorCNA/extdata"
```

Check the required files:

```bash
ls -lh "${EXTDATA}/map_hg38_500kb.wig"
ls -lh "${EXTDATA}/GRCh38.GCA_000001405.2_centromere_acen.txt"
```

The GC reference is stored inside the project:

```bash
PROJECT_DIR="$(pwd)"
GCWIG="${PROJECT_DIR}/ichorcna_autosome/reference/g_hg38_500kb.wig"

ls -lh "${GCWIG}"
```

# 7. Workflows

## 7.1 Workflow 1: Nanopore copy-number profiles

### Scripts

```text
scripts/check_all_bin_width.py
scripts/nanopore_5000.py
```

### Steps

1. Read `CNV_dict.npy` for each sample.
2. Extract the original Nanopore bin width.
3. Save the bin-width summary to:
   
   ```text
   scatter_plot_nanopore/all_bin_width.txt
   ```
4. Read Nanopore `CNV.npy` and `CNV_dict.npy`.
5. Estimate reads per bin.
6. Merge adjacent bins until each output point contains approximately 5000 estimated reads.
7. Generate one merged copy-number profile for each sample.

### Output

```text
scatter_plot_nanopore/
├── all_bin_width.txt
└── 5000/
    └── *_nanopore_merged_reads.png
```

### Example commands

```bash
cd "${SCRIPTS_DIR}"

python check_all_bin_width.py
python nanopore_5000.py
```

## 7.2 Workflow 2: Illumina copy-number profiles

### Scripts

```text
scripts/illumina_check_5000.py
scripts/52.py
```

### Steps

1. Read the standard Illumina `.tumor.target.counts.gz` files with `illumina_check_5000.py`.
2. Read the single BigWig Illumina sample separately with `52.py`.
3. Convert read counts to relative copy number.
4. Merge Illumina bins using the corresponding Nanopore bin-width and target-read settings.
5. Generate one merged Illumina copy-number profile for each sample.

### Output

```text
scatter_plot_illumina/
└── 5000/
    └── *_illumina_merged_reads.png
```

### Example commands

```bash
cd "${SCRIPTS_DIR}"

python illumina_check_5000.py
python 52.py
```

## 7.3 Workflow 3: Correlation and CN-load analysis

### Scripts

```text
scripts/correlation.py
scripts/correlation.sh
```

### Steps

1. Read the Illumina and Nanopore copy-number profiles.
2. Process each autosome separately.
3. Merge bins according to the selected read-based strategy.
4. Apply Gaussian smoothing when enabled.
5. Match Illumina and Nanopore points by genomic position.
6. Calculate Pearson correlation.
7. Calculate copy-number standard deviation, matching distance, and CN-load.
8. Generate summary tables and plots.

### Output

```text
correlation/
├── correlation_summary.csv
├── reads_based_correlation_summary.csv
├── pearson_correlation_by_sample.png
├── illumina_cn_load_vs_pearson.png
└── nanopore_cn_load_vs_pearson.png
```

### Example commands

Run directly:

```bash
cd "${SCRIPTS_DIR}"
python correlation.py
```

Or submit through Slurm:

```bash
cd "${SCRIPTS_DIR}"
sbatch correlation.sh
```

## 7.4 Workflow 4: ichorCNA analysis

### Scripts

```text
scripts/make_gc_hg38_500kb.R
scripts/run_batch_ichorCNA_500kb.sh
```

### Steps

1. Generate the hg38 500 kb GC-content WIG file. (Which has been done in the step 6.3 Generate the GC reference)
2. Keep primary mapped reads from each Nanopore BAM.
3. Remove `@RG` entries from the BAM header.
4. Generate a 500 kb read-count WIG using `readCounter`.
5. Replace scientific notation such as `5e+05` in WIG headers with `500000`.
6. Run GavinHaLab ichorCNA v0.4.0.
7. Save tumour-fraction, ploidy, segmentation, and copy-number results.
8. Remove large intermediate BAM files when the run completes successfully.

### Output

```text
ichorcna_autosome/
├── reference/
│   └── g_hg38_500kb.wig
├── 3_wig/
│   └── 500kb/
├── 4_ichorCNA_results/
│   └── 500kb/
└── logs/
```

## 7.5 HPC job submission

### Array configuration

The current batch script contains 11 BAM entries, so the array range is:

```bash
#SBATCH --array=0-10%1
```

The `%1` limit allows only one array task to run at a time.

### Slurm path handling

The batch script should derive the project directory from the Slurm submission directory:

```bash
SCRIPT_DIR="${SLURM_SUBMIT_DIR}"
PROJECT_DIR="$(dirname "${SCRIPT_DIR}")"
WORKDIR="${PROJECT_DIR}/ichorcna_autosome"
```

Submit the job from the project's `scripts` directory:

```bash
PROJECT_DIR="$(pwd)"
SCRIPTS_DIR="${PROJECT_DIR}/scripts"

cd "${SCRIPTS_DIR}"
sbatch run_batch_ichorCNA_500kb.sh
```

Do not submit the script from an unrelated working directory, because `SLURM_SUBMIT_DIR` would then point to the wrong location.

### Slurm logs

The Slurm output paths can be defined relative to the `scripts` submission directory:

```bash
#SBATCH --output=../ichorcna_autosome/logs/batch_ichorCNA_500kb_%A_%a.out
#SBATCH --error=../ichorcna_autosome/logs/batch_ichorCNA_500kb_%A_%a.err
```

Create the log directory before submitting the job:

```bash
PROJECT_DIR="$(pwd)"
ICHOR_LOG_DIR="${PROJECT_DIR}/ichorcna_autosome/logs"

mkdir -p "${ICHOR_LOG_DIR}"
```

Slurm opens `.out` and `.err` files before the script body runs, so the directory must already exist.

### Submit the ichorCNA batch job

Run from the project root:

```bash
PROJECT_DIR="$(pwd)"
SCRIPTS_DIR="${PROJECT_DIR}/scripts"
ICHOR_LOG_DIR="${PROJECT_DIR}/ichorcna_autosome/logs"

mkdir -p "${ICHOR_LOG_DIR}"
cd "${SCRIPTS_DIR}"

sbatch run_batch_ichorCNA_500kb.sh
```

Check the job:

```bash
squeue -u "${USER}"
```

Cancel a job or array:

```bash
scancel JOB_ID
```

# 8. Important notes

- The analysis only includes `chr1` to `chr22`.
- `chrX`, `chrY`, and mitochondrial DNA are excluded.
- GavinHaLab ichorCNA v0.4.0 was used for the analysis scripts and source code.
- The installed ichorCNA R package v0.5.1 was used to access the bundled hg38 mappability and centromere reference files.

# 9. Tools and versions

## Python environment

- Python 3.10.19
- NumPy 2.2.5
- pandas 2.3.3
- SciPy 1.15.3
- Matplotlib 3.10.6
- pyBigWig 0.3.25
- adjustText 1.3.0

## R environment

- R 4.4.3
- BSgenome.Hsapiens.UCSC.hg38 1.4.5
- Biostrings 2.74.0
- ichorCNA R package 0.5.1

## Command-line tools

- samtools 1.23.1
- HMMcopy command-line tools 0.1.1
- GavinHaLab ichorCNA v0.4.0

# 10. Input and output data

Input data

| Data | Format | Number | Approximate size | Description |
|---|---|---:|---:|---|
| Nanopore CN profiles | `CNV.npy` | 30 | 88 KB–1.2 MB | Per-chromosome copy-number arrays |
| Nanopore metadata | `CNV_dict.npy` | 30 | pproximately 512 bytes | Contains original bin width |
| Illumina target counts | `.tumor.target.counts.gz` | 29 | 28 MB–29 MB | Target-bin read counts |
| Illumina BigWig sample | `.bw` | 1 | 34.5 MB | Alternative format for sample 52 |
| Nanopore alignments | `.bam` | 11 | Approximately 20 GB - 50GB | Input to ichorCNA workflow |

Output data


| Workflow | Output directory | Key files | Purpose |
|---|---|---|---|
| 1. Nanopore CN profiles | `scatter_plot_nanopore/5000/` | `*_nanopore_merged_reads.png` | Per-sample CN scatter plot merged to ~5000 reads per point |
| 2. Illumina CN profiles | `scatter_plot_illumina/5000/` | `*_illumina_merged_reads.png` | Per-sample CN scatter plot merged to match the Nanopore setting |
| 3. Correlation & CN-load | `correlation/` | `correlation_summary.csv`<br>`reads_based_correlation_summary.csv`<br>`pearson_correlation_by_sample.png`<br>`illumina_cn_load_vs_pearson.png`<br>`nanopore_cn_load_vs_pearson.png`| Per-sample Illumina/Nanopore Pearson correlation, CN-load metrics, and associated plots |
| 4. ichorCNA | `ichorcna_autosome/4_ichorCNA_results/500kb/{sample}/` | `*_500kb.params.txt` (tumor fraction/ploidy)<br>`*_500kb.cna.seg` / `*.seg` / `*.seg.txt` (CN segmentation)<br>`*_500kb.correctedDepth.txt` (corrected read depth)<br>`*_500kb/*_genomeWide.pdf` (genome-wide CN plot) | Per-sample tumor fraction and ploidy estimates, CN segmentation results, and diagnostic plots |
