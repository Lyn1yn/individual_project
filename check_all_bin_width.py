import numpy as np
from pathlib import Path

# input: shared data directory
samples_dir = Path("/gpfs01/share/BioinfMSc/Matt_Projects/samples")

PROJECT_DIR = Path(__file__).resolve().parent.parent

output_file = PROJECT_DIR / "scatter_plot_nanopore" / "all_bin_width.txt"


with open(output_file, "w") as f:
    for sample_dir in sorted(samples_dir.iterdir()):
        if not sample_dir.is_dir():
            continue

        sample = sample_dir.name
        cnv_dict_file = sample_dir / "nanopore" / "CNV_dict.npy"

        cnv_dict = np.load(cnv_dict_file, allow_pickle=True).item()
        bin_width = cnv_dict["bin_width"]

        print(sample, bin_width)
        f.write(f"{sample} {bin_width}\n")
