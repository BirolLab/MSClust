#!/bin/bash
export PATH="/projects/btl_scratch/pkazemi/miniforge/envs/nanosim/bin:$PATH"

REFERENCE_GENOME='/projects/ABySS_research/hackathons/singlecellmet2025/read_methylation/hg002/hg002v1.0.1.fasta'
ANALYSIS_OUTPUT_DIR='/projects/ABySS_research/hackathons/singlecellmet2025/read_methylation/hg002/nanosim_profiles/human_giab_hg002_sub1M_kitv14_dorado_v3.2.1/training'
NUM_THREADS=96

SIMULATION_OUTPUT_DIR='/projects/btl_scratch/bucar/scMethylationReads/scMeth'
SEED=123
NUMBER_OF_READS=1000000

simulator.py genome -rg $REFERENCE_GENOME -c $ANALYSIS_OUTPUT_DIR -o $SIMULATION_OUTPUT_DIR --seed $SEED -t $NUM_THREADS -n $1