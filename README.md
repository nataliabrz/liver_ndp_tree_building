# liver_ndp_tree_building
This repository contains scripts for running n-dimensional dirichlet process (NDP) to cluster mutations based on VAFs from multiple related samples, and building phylogenetic trees based on SNV clustering.

## Workflow
1. Pull extended sequence context for SNVs and indels.
2. Perform n-dimensional dirichlet process (NDP) clustering on SNVs.
3. Generate phylogenetic trees based on clusters.
4. Assign indels to NDP clusters.

## Requirements
- R (>= 4.3.1)
- Required R libraries: `data.table`, `dplyr`, `BSgenome.Hsapiens.UCSC.hg38`, `GenomicRanges`, `Biostrings`
- HPC with job scheduler (e.g., LSF)

### Installation
1. Clone the repository:
```bash
git clone https://github.com/nataliabrz/tree-building-pipeline.git
cd tree-building-pipeline
```

2. Install R dependencies:
```R
install.packages(c(
    "data.table", "dplyr", "stringr", "RColorBrewer",
    "label.switching", "philentropy", "ggplot2", "GGally", 
    "here", "farver", "tibble", "tidyr", "gridExtra", 
    "igraph", "ape", "ggrepel"
))

if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}
BiocManager::install(c(
    "BSgenome.Hsapiens.UCSC.hg38", "GenomicRanges", 
    "Biostrings"
))
```


## Running the Pipeline

### Example Workflow
#### Step 1: Extract the 10 bp sequence context around each mutation for SNVs and indels. This file will be saved in the same directory as input calls.
```bash
Rscript scripts/1_pull_context.R
```

#### Step 2: Run NDP Clustering
```bash
bash scripts/2_ndp_clustering.sh \
  data/filtered_calls/snv \
  data/metadata/sex_file.txt \
  outputs/ndp_results \
  PD51606 \
  15000
```

#### Step 3: Generate phylogenetic trees from NDP clustering outputs and SNV files.
```bash
bash scripts/3_tree_building.sh \
  PD51606 \
  outputs/ndp_results/PD51606 \
  data/filtered_calls/snv/PD51606_bb_pass_snvs_all.csv \
  outputs/trees/PD51606 \
  0.10 \
  50
```

#### Step 4: Assign indels to NDP clusters and output a complete mutation file.
```bash
Rscript scripts/4_allocate_indels_to_ndp_clusters.R

```
