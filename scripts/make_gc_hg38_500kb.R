library(BSgenome.Hsapiens.UCSC.hg38)
library(Biostrings)
options(scipen = 999)


#output
script_dir <- dirname(normalizePath(
  sub(
    "--file=",
    "",
    grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  )
))

PROJECT_DIR <- dirname(script_dir)

out_dir <- file.path(
  PROJECT_DIR,
  "ichorcna_autosome",
  "reference"
)

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

out_file <- file.path(out_dir, "g_hg38_500kb.wig")




genome <- BSgenome.Hsapiens.UCSC.hg38 #change genome object into a shorter name

window <- 500000 #bin width of 500kb, matching the ichorCNA solution
chromosomes <- paste0("chr", c(1:22)) #only select autosomes

con <- file(out_file, open = "w")

for (chr in chromosomes) {
  seq_chr_full <- genome[[chr]] #get a DNAString
  chr_len <- length(seq_chr_full) #get the length of chromosomes
  starts <- seq(1, chr_len, by = window) #get the position of bin start

  cat(
    paste0(
      "fixedStep chrom=", chr,
      " start=1 step=", window,
      " span=", window,
      "\n"
    ),
    file = con
  ) #make a string, and put the string into a file

  for (start in starts) {
    end <- min(start + window - 1, chr_len) #count the end of the bin
    seq_bin <- subseq(seq_chr_full, start = start, end = end) #take out from start to end sequence in the whole chromosome

    counts <- letterFrequency(
      seq_bin,
      letters = c("G", "C", "A", "T"), #count how many times letters show up in the sequence
      as.prob = FALSE
    )

    total_acgt <- sum(counts) #only ATCG, no N 

    if (total_acgt == 0) {
      gc <- 0
    } else {
      gc <- (counts["G"] + counts["C"]) / total_acgt #calculate GC percentage
    }

    cat(sprintf("%.6f\n", gc), file = con)
  }
}

close(con)

cat("GC wig written to:\n")
cat(out_file, "\n")
