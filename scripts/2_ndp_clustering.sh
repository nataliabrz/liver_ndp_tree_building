#!/bin/bash

# Variables
INPUT_DIR=$1          # Path to beta-binomial analysis directory (top level)
SEX_FILE=$2           # Path to sex file
OUT_DIR=$3            # Path to NDP output directory
PATIENT=$4            # Patient ID
burnin=$5             # Number of burn-in iterations

# Create required directories
mkdir -p $OUT_DIR/logfiles
mkdir -p $OUT_DIR/$PATIENT

# Cluster resources
MEM=264000            # Memory in MB
CORE=4                # Number of cores
q=long                # Queue type

# Submit the job
bsub -G team154-grp -J dirichlet_setup_$PATIENT \
    -M"$MEM" -R"select[mem>$MEM] rusage[mem=$MEM] span[hosts=1]" \
    -n $CORE -q $q \
    -e $OUT_DIR/logfiles/$PATIENT.stderr \
    -o $OUT_DIR/logfiles/$PATIENT.stdout \
    Rscript scripts/2a_run_dirichlet_clustering_posthoc.R \
        scripts/utils/ \
        $INPUT_DIR/${PATIENT}_ndp_alt_bb_flt.csv \
        $INPUT_DIR/${PATIENT}_ndp_depth_bb_flt.csv \
        $INPUT_DIR/${PATIENT}_mut_context_GRCh38.txt \
        $SEX_FILE \
        $OUT_DIR/$PATIENT \
        $burnin
