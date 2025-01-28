#!/bin/bash

# Script: 3_tree_building.sh
# Description: Wrapper script to submit phylogenetic tree-building jobs.
# Usage:
#   bash scripts/3_tree_building.sh \
#     <patientID> <ndp_input_dir> <snv_file> <tree_out_dir> [<min_vaf_threshold>] [<min_mut_count>]

# Input arguments
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"

patientID="$1"                               # e.g., PD51606
ndp_input_dir="$2"                           # Directory containing NDP clustering outputs
snv_file="$3"                                # Path to SNV file
tree_out_dir="$4"                            # Directory for tree outputs
min_vaf_threshold="${5:-0.10}"               # Default: 0.10
min_mut_count="${6:-50}"                     # Default: 50

# Validate input arguments
if [ -z "$patientID" ] || [ -z "$ndp_input_dir" ] || [ -z "$snv_file" ] || [ -z "$tree_out_dir" ]; then
  echo "Error: Missing required arguments."
  echo "Usage: bash $0 <patientID> <ndp_input_dir> <snv_file> <tree_out_dir> [<min_vaf_threshold>] [<min_mut_count>]"
  exit 1
fi

echo "Patient ID: $patientID"
echo "NDP input directory: $ndp_input_dir"
echo "SNV file: $snv_file"
echo "Tree output directory: $tree_out_dir"
echo "Min VAF threshold: $min_vaf_threshold"
echo "Min mutation count: $min_mut_count"
echo "----"

# Create output directories
mkdir -p "$tree_out_dir"
mkdir -p "$tree_out_dir/logfiles"

# Check for bsub (LSF HPC environment)
if command -v bsub &> /dev/null; then
  echo "Submitting job to LSF cluster..."
  bsub -G team154-grp -J "$patientID.ndp_tree" \
      -M8000 -R"select[mem>8000] rusage[mem=8000] span[hosts=1]" \
      -n 1 -q normal \
      -e "$tree_out_dir/logfiles/$patientID.stderr" \
      -o "$tree_out_dir/logfiles/$patientID.stdout" \
      "Rscript $SCRIPT_DIR/3a_ndp_tree_generation.R \
          $patientID \
          $ndp_input_dir \
          $SCRIPT_DIR \
          $tree_out_dir \
          $snv_file \
          $min_vaf_threshold \
          $min_mut_count"
else
  echo "Warning: bsub not found. Running locally..."
  Rscript "$SCRIPT_DIR/3a_ndp_tree_generation.R" \
      "$patientID" \
      "$ndp_input_dir" \
      "$SCRIPT_DIR" \
      "$tree_out_dir" \
      "$snv_file" \
      "$min_vaf_threshold" \
      "$min_mut_count"
fi
