import os
import glob
import Bio.SeqIO.FastaIO as FIO
import random

random.seed(123)
INPUT_FILE = ""
OUTPUT_DIR = "/projects/ABySS_research/hackathons/singlecellmet2025/read_methylation/reads"

# read the folder to get the names for cell types:
cell_types=glob.glob("/projects/ABySS_research/grant_prep/R21_single-cell-methylation/data/selected_meth_profiles/*.bed")
cell_types = ["-".join(cell_type.split("/")[-1].split("_")[1].split(".")[0].split("-")[:-1]) for cell_type in cell_types]

handles = []
writers = []
for cell_type in cell_types:
    if not os.path.exists(os.path.join(OUTPUT_DIR, cell_type)):
        os.mkdir(os.path.join(OUTPUT_DIR, cell_type))
    for i in range(20):
        handles.append(open(os.path.join(OUTPUT_DIR, cell_type,i), "w"))
        writers.append(FIO.FastaWriter(handles[i]))
        writer[i].write_header()

with open(INPUT_FILE) as handle:
    for values in FIO.SimpleFastaParser(handle):
        random_int = random.randint(0, 399) # 20*20 cells (20 cells for each 20 cell types)
        writers[random_int].write_record(values)

for i in range(400):
    writers[i].write_footer() 
    handles[i].close()