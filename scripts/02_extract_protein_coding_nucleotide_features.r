library(ftrCOOL)
library(seqinr)
library(stringr)

# Source custom functions
source("calculate_gc_percentage.r")
source("calculate_kozak_strength.r")
source("calculate_translation_initiation_efficiency.r")
source("compute_hexamerScore.r")

# Load reference frequency data
reference_frequency <- readRDS(file = "../data/input/referFreq.rds")

# Set filtering parameters
MINIMUM_UTR5_LENGTH <- 10
MINIMUM_CDS_LENGTH <- 30
MINIMUM_UTR3_LENGTH <- 10

# Load and filter protein-coding sequences
process_protein_coding_sequences <- function(file_path) {
  protein_coding_sequences <- read.table(file_path, sep = "\t", quote = "")
  filtered_sequences <- protein_coding_sequences[
    protein_coding_sequences$V4 == "protein_coding" &
    nchar(protein_coding_sequences$V7) >= MINIMUM_UTR5_LENGTH &
    nchar(protein_coding_sequences$V8) >= MINIMUM_CDS_LENGTH &
    nchar(protein_coding_sequences$V9) >= MINIMUM_UTR3_LENGTH &
    grepl("^ATG", protein_coding_sequences$V8),
  ]
  return(filtered_sequences)
}

# Extract sequence features
extract_sequence_features <- function(protein_sequences) {
  # Separate sequences
  utr5_sequences <- protein_sequences$V7
  cds_sequences <- protein_sequences$V8
  utr3_sequences <- protein_sequences$V9
  
  
    # Calculate hexamer scores
  hexamer_scores <- compute_hexamerScore(
    sapply(cds_sequences, function(x) s2c(tolower(x))), 
    referFreq = reference_frequency,
    k = 6,
    step = 1,
    alphabet = c("a", "c", "g", "t"),
    on.ORF = FALSE,
    auto.full = TRUE,
    parallel.cores = 60
  )

  
  
  # Calculate GC content
  utr5_gc_percentages <- unname(sapply(utr5_sequences, calculate_gc_percentage))
  cds_gc_percentages <- unname(sapply(cds_sequences, calculate_gc_percentage))
  utr3_gc_percentages <- unname(sapply(utr3_sequences, calculate_gc_percentage))
  
  # Calculate sequence lengths and percentages
  utr5_lengths <- unname(sapply(utr5_sequences, nchar))
  cds_lengths <- unname(sapply(cds_sequences, nchar))
  utr3_lengths <- unname(sapply(utr3_sequences, nchar))
  
  total_lengths <- utr5_lengths + cds_lengths + utr3_lengths
  utr5_length_percentages <- utr5_lengths / total_lengths
  cds_length_percentages <- cds_lengths / total_lengths
  utr3_length_percentages <- utr3_lengths / total_lengths
  
  
  # Calculate sequence scores
  kozak_strengths <- calculate_kozak_strength(utr5_sequences, cds_sequences)
  translation_initiation_efficiencies <- calculate_translation_initiation_efficiency(utr5_sequences, cds_sequences)
  fickett_scores <- fickettScore(seqs = cds_sequences)
  
  # Calculate regulatory elements
  are_presence <- as.integer(grepl("ATTTA", utr3_sequences))
  pas_presence <- as.integer(grepl("AATAAA", utr3_sequences))
  
  
  # Combine features
  sequence_features <- data.frame(
    row.names = paste(protein_sequences$V1, protein_sequences$V2, protein_sequences$V3,
                     protein_sequences$V4, protein_sequences$V6, sep = "_"),
    hexamer_score = hexamer_scores,
    fickett_score = fickett_scores,
    translation_initiation_efficiency = translation_initiation_efficiencies,
    kozak_strength = kozak_strengths,
    utr5_gc_percentage = utr5_gc_percentages,
    cds_gc_percentage = cds_gc_percentages,
    utr3_gc_percentage = utr3_gc_percentages,
    are_element = are_presence,
    pas_element = pas_presence,
    utr5_length_percentage = utr5_length_percentages,
    cds_length_percentage = cds_length_percentages,
    utr3_length_percentage = utr3_length_percentages
  )
  
  return(sequence_features)
}

# Main execution
protein_sequences <- process_protein_coding_sequences("../data/input/PC.txt")
protein_features <- extract_sequence_features(protein_sequences)

# Save results
write.table(protein_features, "../data/intermediate/protein_coding_nucleotide_features.txt", sep = "\t", quote = FALSE)
