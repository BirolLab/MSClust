import joblib
import numpy as np
import os

# constants - paramaters - 76603200 reads
N_ESTIMATION = 10000000
GENOME_FILE = "/projects/ABySS_research/hackathons/singlecellmet2025/read_methylation/hg002/hg002v1.0.1.fasta"
KDE_PATH = '/projects/ABySS_research/hackathons/singlecellmet2025/read_methylation/hg002/nanosim_profiles/human_giab_hg002_sub1M_kitv14_dorado_v3.2.1/training_aligned_region.pkl'# I am going to be using aligned portion kde of the aligned reads as we calculate coverage based on aligned reads

if GENOME_FILE != None:
    GENOME_SIZE = 0
    with open(GENOME_FILE, "r") as f_in:
        for line in f_in:
            line = line.strip()
            if line.startswith(">"):
                continue 
            else:
                GENOME_SIZE += len(line)
        
# read the Kernel Density function of reads (both aligned and unaligned)
np.float = float
read_len = []
kde = joblib.load(KDE_PATH)

# sample read sizes for estimating the mean/expected value in this case of kernel density
samples = kde.sample(N_ESTIMATION)
mean = samples.mean()
print(f"Mean length of the aligned portions: {mean}")

# calculate number of reads based on genome size
read_cnt = int(GENOME_SIZE / mean)
print(f"The number of reads for 1x coverage: {read_cnt}")
read_cnt = read_cnt*20*10 # 20 cell types and 10x each
print(f"The number of reads for 20 cell types (10x for cell type): {read_cnt}")

# run nanosim
os.system(f"bash /projects/ABySS_research/hackathons/singlecellmet2025/read_methylation/hg002/run_nanosim.sh {read_cnt}")  
