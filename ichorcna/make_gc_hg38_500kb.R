library(BSgenome.Hsapiens.UCSC.hg38)
library(Biostrings)

genome <- BSgenome.Hsapiens.UCSC.hg38

out_file <- "/gpfs01/home/mbxll1/CNS_cancer_project/ichorCNA/resources/g_hg38_500kb.wig"

window <- 500000
chromosomes <- paste0("chr", c(1:22, "X"))

con <- file(out_file, open = "w")

for (chr in chromosomes) {
  seq_chr_full <- genome[[chr]]
  chr_len <- length(seq_chr_full)
  starts <- seq(1, chr_len, by = window)

  cat(
    paste0(
      "fixedStep chrom=", chr,
      " start=1 step=", window,
      " span=", window,
      "\n"
    ),
    file = con
  )

  for (start in starts) {
    end <- min(start + window - 1, chr_len)
    seq_bin <- subseq(seq_chr_full, start = start, end = end)

    counts <- letterFrequency(
      seq_bin,
      letters = c("G", "C", "A", "T"),
      as.prob = FALSE
    )

    total_acgt <- sum(counts)

    if (total_acgt == 0) {
      gc <- 0
    } else {
      gc <- (counts["G"] + counts["C"]) / total_acgt
    }

    cat(sprintf("%.6f\n", gc), file = con)
  }
}

close(con)

cat("GC wig written to:\n")
cat(out_file, "\n")
