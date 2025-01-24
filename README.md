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
    "RColorBrewer", "label.switching", "philentropy", 
    "ggplot2", "GGally", "dplyr", "data.table", 
    "here", "farver", "gplots"
))

if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}
BiocManager::install(c("BSgenome.Hsapiens.UCSC.hg38", "GenomicRanges", "Biostrings"))

```


## Running the Pipeline

### Example Workflow
#### Step 1: Extract the 10 bp sequence context around each mutation for SNVs and indels.
```bash
Rscript scripts/1_pull_context.R
```

#### Step 2: Run NDP Clustering
```bash
patient=PD51606
bash scripts/2_ndp_clustering.sh \
  data/filtered_calls/snv \
  data/metadata/sex_file.txt \
  outputs/ndp_results \
  $patient \
  15000
```

#### Step 3: Generate Trees
```bash
bash scripts/3_generate_trees.sh --input_clusters outputs/ndp_out/clusters.txt \
  --output_tree outputs/trees/tree.nwk
```

#### Step 4: Assign SNVs and Indels to Clusters
```bash
Rscript scripts/4_assign_clusters_snvs.R --input_snv data/example_snv_data/snv_data.csv \
  --input_tree outputs/trees/tree.nwk --output assignments_snv.csv

Rscript scripts/5_assign_clusters_indels.R --input_indel data/example_indel_data/indel_data.csv \
  --input_tree outputs/trees/tree.nwk --output assignments_indel.csv
```
