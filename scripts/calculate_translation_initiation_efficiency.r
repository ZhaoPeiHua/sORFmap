calculate_translation_initiation_efficiency <- function(utr5_sequences, cds_sequences) {
   # Load TIS efficiency data
   tis_efficiency_table <- read.table("../data/input/TIS_efficiency.txt", header = TRUE)
   rownames(tis_efficiency_table) <- tis_efficiency_table$sequence
   tis_efficiency_table <- tis_efficiency_table[, -1]
   
   # Extract context sequences
   context_minus6 <- substr(utr5_sequences, nchar(utr5_sequences) - 5, nchar(utr5_sequences))
   context_plus5 <- substr(cds_sequences, 1, 5)
   
   # Combine and convert to RNA notation
   tis_contexts <- paste0(context_minus6, context_plus5)
   tis_contexts <- gsub("T", "U", tis_contexts)
   
   # Retrieve efficiency scores
   efficiency_scores <- tis_efficiency_table[tis_contexts, "efficiency"]
   
   return(efficiency_scores)
}