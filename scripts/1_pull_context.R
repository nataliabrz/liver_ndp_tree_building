
library(data.table)
library(dplyr)
library(BSgenome.Hsapiens.UCSC.hg38)
library(GenomicRanges)
library(Biostrings)
library(stringr)

# Define types of mutations
types <- c("snv", "indel")

# Define base paths relative to repository structure
input_base <- "data/filtered_calls"
output_base <- "outputs/ndp_input"

for (type in types) {
  # Define type-specific input and output directories
  pull_dir <- file.path(input_base, type)
  out_dir <- file.path(output_base, type)

  # Create output directory if it doesn't exist
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)
  }

  # Get list of input files
  pull_files <- list.files(pull_dir,
                           full.names = TRUE,
                           recursive = TRUE,
                           pattern = "_ndp_alt_bb_flt.csv")

  print(type)

  ### Loop through pull files, generate and write outputs
  for (i in seq_along(pull_files)) {
    # Read the variant file
    this_file <- pull_files[i]
    patient <- gsub(".*?(PD\\d{5})\\w{1}.*", "\\1", pull_files[i])
    print(paste0("Pulling context for patient ", patient))

    this_pos_list <- fread(this_file, sep = ",")
    this_pos_list <- this_pos_list[, c("chrom", "pos", "ref", "alt")]

    # Keep only unique mutations
    this_pos_list <- this_pos_list[!duplicated(this_pos_list), ]

    # Define the context (10 bp each side)
    context_list <- this_pos_list %>%
      dplyr::mutate(start = pos - 10,
                    end = pos + 10) %>%
      dplyr::select(chrom, start, end)

    # Get the sequences
    this_range <- GenomicRanges::makeGRangesFromDataFrame(context_list)
    this_seq <- getSeq(BSgenome.Hsapiens.UCSC.hg38, this_range)
    out_seqs <- GenomicRanges::as.data.frame(this_seq)
    out_table  <- bind_cols(this_pos_list, out_seqs)
    colnames(out_table ) <- c("chrom", "pos", "ref", "alt", "CONTEXT")

    # Write the table to file
    write.table(
      out_table,
      file = file.path(out_dir, paste0(patient, "_mut_context_GRCh38.txt")),
      quote = FALSE,
      row.names = FALSE,
      sep = "\t"
    )
  }
}
