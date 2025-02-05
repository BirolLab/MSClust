import os 

with open("/projects/btl_scratch/bucar/scMethylationReads/aligned_error_profile", "r") as f_in:
    for line in f_in:
        name = line.split("\t")[0]
        with open(f"/projects/btl_scratch/bucar/scMethylationReads/errors/{name}", "a ") as f_out:
            f_out.write(line)


