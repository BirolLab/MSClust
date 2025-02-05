import os
import csv
import argparse
from Bio import SeqIO
import numpy as np
import time

def parse_args():
    """Reads input arguments from command line"""
    parser = argparse.ArgumentParser()
    parser.add_argument("-r", "--reads", type=str, required=True, help="FASTA file with reads")
    parser.add_argument("-m", "--meth", type=str, required=False, help="Methylation marks (BED format)")
    parser.add_argument("-d", "--error_dir", type=str, required=True, help="Directory containing error TSV files")
    parser.add_argument("-o", "--output", type=str, required=False, default="methylated_reads.fa", help="Output FASTA file")
    return parser.parse_args()

def load_methylation_marks(meth_file, haplotype):
    """Loads methylation marks from BED file into a dictionary (chrom -> sorted list of positions)."""
    print("Loading methylation marks...")
    meth_table = {}
    with open(meth_file, "r", encoding="utf-8") as f:
        reader = csv.reader(f, delimiter="\t")
        for row in reader:
            chrom = row[0]
            if haplotype == "maternal":
                pos = int(row[1])
            else:
                pos = int(row[2])
            prob = float(row[3])
            if chrom not in meth_table:
                meth_table[chrom] = []
            if pos == -1 or prob < 0.5: # if it is not aligned to HG002 or if the frequency of methylation is smaller than 0.5 do not consider
                continue
            meth_table[chrom].append(pos)
    return meth_table

def load_error_profile(error_dir, seq_name):
    """Loads error profile for a specific read from a TSV file."""
    error_file = os.path.join(error_dir, seq_name)  # No file extension
    errors = []
    
    if os.path.exists(error_file):
        with open(error_file, "r", encoding="utf-8") as f:
            reader = csv.reader(f, delimiter="\t")
            for row in reader:
                try:
                    error_pos = int(row[1])
                    error_type = row[2]
                    error_length = int(row[3])
                    errors.append((error_pos, error_type, error_length))
                except ValueError:
                    print(f"Skipping malformed error entry in {seq_name}")
    errors.reverse()
    return errors # Return empty list if no errors exist

def apply_methylation_and_errors(seq_name, seq, meth_table, error_dir):
    """Applies methylation and error shifts to the sequence.
        This function assumes that there is no overlapping errors in a single read."""
    split_seq_name = seq_name.split('_')
    chrom = split_seq_name[0].split('-')[0]  # Extract chromosome name
    start_on_read = int(split_seq_name[1])
    end_on_read = start_on_read + int(split_seq_name[6])
    strand = split_seq_name[4]
    clipped_length = int(split_seq_name[5])

    # Get methylation marks (if available)
    if chrom not in meth_table:
        return seq  # No modifications needed

    # Convert sequence to mutable NumPy array (FASTER than string operations)
    seq_arr = np.array(list(seq), dtype="U1")

    # Load errors dynamically from TSV
    errors = load_error_profile(error_dir, seq_name)
    shift = 0
    error_idx = 0
    len_errors = len(errors)
    for meth_pos in meth_table[chrom]:    
        skip = False
        if meth_pos > end_on_read:
            break
        if start_on_read <= meth_pos:  # Methylation lies in the read region
            #for error_pos, error_type, error_length in errors:
            while error_idx >= len_errors:
                error_pos, error_type, error_length = errors[error_idx]   
                
                if error_pos + start_on_read <= meth_pos:
                    if error_type == 'ins':
                        shift += error_length
                        error_idx += 1
                    elif error_type == 'del':
                        if start_on_read + error_pos + error_length - 1 >= meth_pos:
                            skip = True  # Skip this site due to deletion
                            break
                        else:
                            shift -= error_length
                            error_idx += 1
                else:
                    break


            index_of_meth = clipped_length + meth_pos - start_on_read + shift
            if 0 <= index_of_meth < len(seq_arr) - 1 and not skip:
                if strand == 'F':
                    if seq_arr[index_of_meth].upper() == 'C' and seq_arr[index_of_meth + 1].upper() == 'G':
                        seq_arr[index_of_meth] = '1'
                else:
                    rev_index = len(seq_arr) - 1 - index_of_meth
                    if 0 <= rev_index < len(seq_arr) - 1:
                        if seq_arr[rev_index].upper() == 'G' and seq_arr[rev_index - 1].upper() == 'C':
                            seq_arr[rev_index - 1] = '1'

    return "".join(seq_arr)

def process_reads(args):
    """Main function to process FASTA reads with methylation and error corrections."""
    print("Starting processing...")
    meth_table_pat = load_methylation_marks(args.meth, "paternal") if args.meth else {}
    meth_table_mat = load_methylation_marks(args.meth, "maternal") if args.meth else {}
    start = time.time()

    with open(args.output, "w") as writer_fasta:
        count = 0
        for record in SeqIO.parse(args.reads, "fasta"):
            count += 1 
            if count % 100000 == 0:
                print(f"{count / 191251 * 100}% in {(time.time()-start)/60} mins.")
            seq_name, seq = record.id, str(record.seq)
            if "MATERNAL" in seq_name:
                updated_seq = apply_methylation_and_errors(seq_name, seq, meth_table_mat, args.error_dir)
            else:
                updated_seq = apply_methylation_and_errors(seq_name, seq, meth_table_pat, args.error_dir)
            writer_fasta.write(f">{seq_name}\n{updated_seq.upper()}\n")

    print("Processing complete!")

if __name__ == "__main__":
    args = parse_args()
    process_reads(args)

