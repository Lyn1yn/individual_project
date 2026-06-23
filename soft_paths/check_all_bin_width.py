import numpy as np
from pathlib import Path

# Shared data directory - unchanged
samples_dir = Path("/gpfs01/share/BioinfMSc/Matt_Projects/samples")

for sample_dir in sorted(samples_dir.iterdir()):
    if not sample_dir.is_dir():
        continue

    sample = sample_dir.name
    cnv_dict_file = sample_dir / "nanopore" / "CNV_dict.npy"

    if not cnv_dict_file.exists():
        print(sample, "NO_CNV_DICT")
        continue

    cnv_dict = np.load(cnv_dict_file, allow_pickle=True).item()
    bin_width = cnv_dict["bin_width"]

    print(sample, bin_width)
