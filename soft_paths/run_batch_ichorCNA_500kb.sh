#!/bin/bash
#SBATCH --job-name=batch_ichorCNA_500kb
#SBATCH --output=batch_ichorCNA_500kb_%A_%a.out
#SBATCH --error=batch_ichorCNA_500kb_%A_%a.err
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --array=0-3%1

set -eo pipefail

source ~/.bashrc

# Derive WORKDIR from the script's own location, so it works for any user.
# Assumes this script lives at: <project_root>/ichorcna/run_batch_ichorCNA_500kb.sh
# WORKDIR resolves to:          <project_root>/ichorcna/
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
WORKDIR="${SCRIPT_DIR}"

# Shared raw BAM directory - unchanged
RAW_BAM_DIR="/gpfs01/share/BioinfMSc/Matt_Projects/Nanopore"

BAMS=(
"sort_Intraop0006_A.bam"
"sort_Intraop0021_c.bam"
"sort_Intraop0034_b.bam"
"sort_ds1305_CNVIntraop0067_a.bam"
)

BAM_NAME="${BAMS[$SLURM_ARRAY_TASK_ID]}"

# Strip "sort_" prefix and ".bam" suffix to get sample name
SAMPLE="${BAM_NAME#sort_}"
SAMPLE="${SAMPLE%.bam}"

RAW_BAM="${RAW_BAM_DIR}/${BAM_NAME}"

FILTER_DIR="${WORKDIR}/2_primary_mapped"
NORG_DIR="${WORKDIR}/3_remove_RG"
WIG_DIR="${WORKDIR}/4_wig/500kb"
OUTDIR="${WORKDIR}/5_ichorCNA_results/500kb/${SAMPLE}"
LOGDIR="${WORKDIR}/logs"

mkdir -p "${FILTER_DIR}" "${NORG_DIR}" "${WIG_DIR}" "${OUTDIR}" "${LOGDIR}"

THREADS=8
WINDOW=500000
MAPQ=20

PRIMARY_BAM="${FILTER_DIR}/${SAMPLE}.primary_mapped.bam"
NORG_BAM="${NORG_DIR}/${SAMPLE}.primary_mapped.noRG.bam"
HEADER="${NORG_DIR}/${SAMPLE}.noRG.header.sam"
OUTWIG="${WIG_DIR}/${SAMPLE}.${WINDOW}.mapq${MAPQ}.wig"

CHRS="chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX"

echo "=========================================="
echo "Sample: ${SAMPLE}"
echo "Raw BAM: ${RAW_BAM}"
echo "=========================================="

echo "Check raw BAM"
ls -lh "${RAW_BAM}"

echo "Step 1: make primary mapped BAM"
conda activate ichorcna

if [ ! -f "${PRIMARY_BAM}" ]; then
    samtools view -@ "${THREADS}" -b -F 2308 "${RAW_BAM}" > "${PRIMARY_BAM}"
    samtools index -@ "${THREADS}" "${PRIMARY_BAM}"
fi

echo "Step 2: remove @RG header"
if [ ! -f "${NORG_BAM}" ]; then
    samtools view -H "${PRIMARY_BAM}" | grep -v "^@RG" > "${HEADER}"
    samtools reheader "${HEADER}" "${PRIMARY_BAM}" > "${NORG_BAM}"
    samtools index -@ "${THREADS}" "${NORG_BAM}"
fi

echo "Step 3: make 500kb WIG"
if [ ! -f "${OUTWIG}" ]; then
    readCounter \
        -w "${WINDOW}" \
        -q "${MAPQ}" \
        -c "${CHRS}" \
        "${NORG_BAM}" > "${OUTWIG}"
fi

# Prevent WIG header from showing 5e+05, which ichorCNA cannot parse
sed -i 's/step=5e+05/step=500000/g; s/span=5e+05/span=500000/g' "${OUTWIG}"

echo "WIG generated:"
ls -lh "${OUTWIG}"
head "${OUTWIG}"


echo "Step 4: run ichorCNA"

conda deactivate
conda activate ichorcna_R

ICHOR_REPO="${WORKDIR}/software/GavinHaLab_ichorCNA_v0.4.0"
RUNICHOR="${ICHOR_REPO}/scripts/runIchorCNA.R"

# Reference files: shared extdata from conda env (user-independent path)
EXTDATA="$(conda info --base)/envs/ichorcna_R/lib/R/library/ichorCNA/extdata"

GCWIG="${WORKDIR}/reference/g_hg38_500kb.wig"
MAPWIG="${EXTDATA}/map_hg38_500kb.wig"
CENTROMERE="${EXTDATA}/GRCh38.GCA_000001405.2_centromere_acen.txt"

echo "Check ichorCNA input files"
ls -lh "${RUNICHOR}"
ls -lh "${OUTWIG}"
ls -lh "${GCWIG}"
ls -lh "${MAPWIG}"
ls -lh "${CENTROMERE}"

Rscript "${RUNICHOR}" \
  --libdir "${ICHOR_REPO}" \
  --id "${SAMPLE}_500kb" \
  --WIG "${OUTWIG}" \
  --gcWig "${GCWIG}" \
  --mapWig "${MAPWIG}" \
  --centromere "${CENTROMERE}" \
  --genomeBuild hg38 \
  --genomeStyle UCSC \
  --ploidy "c(2,3)" \
  --normal "c(0.3,0.4,0.5,0.6,0.7,0.8,0.9)" \
  --maxCN 5 \
  --includeHOMD False \
  --chrs "c(\"chr1\",\"chr2\",\"chr3\",\"chr4\",\"chr5\",\"chr6\",\"chr7\",\"chr8\",\"chr9\",\"chr10\",\"chr11\",\"chr12\",\"chr13\",\"chr14\",\"chr15\",\"chr16\",\"chr17\",\"chr18\",\"chr19\",\"chr20\",\"chr21\",\"chr22\",\"chrX\")" \
  --chrTrain "c(\"chr1\",\"chr2\",\"chr3\",\"chr4\",\"chr5\",\"chr6\",\"chr7\",\"chr8\",\"chr9\",\"chr10\",\"chr11\",\"chr12\",\"chr13\",\"chr14\",\"chr15\",\"chr16\",\"chr17\",\"chr18\",\"chr19\",\"chr20\",\"chr21\",\"chr22\")" \
  --estimateNormal True \
  --estimatePloidy True \
  --estimateScPrevalence True \
  --scStates "c(1,3)" \
  --txnE 0.9999 \
  --txnStrength 10000 \
  --plotFileType pdf \
  --plotYLim "c(-2,4)" \
  --outDir "${OUTDIR}/" \
  > "${LOGDIR}/${SAMPLE}_500kb.ichorCNA.log" 2>&1


echo "Remove large intermediate BAM files to save quota"
rm -f "${PRIMARY_BAM}" "${PRIMARY_BAM}.bai"
rm -f "${NORG_BAM}" "${NORG_BAM}.bai"
rm -f "${HEADER}"


echo "Done sample: ${SAMPLE}"
echo "Output directory:"
ls -lh "${OUTDIR}"
