calculate_kozak_strength <- function(utr5_sequences, cds_sequences) {
   # Extract context nucleotides for Kozak sequence analysis
   context_nucleotide_minus3 <- substr(utr5_sequences, nchar(utr5_sequences) - 2, nchar(utr5_sequences) - 2)
   context_nucleotide_plus4 <- substr(cds_sequences, 4, 4)
   
   # Initialize Kozak strength scores
   kozak_strengths <- rep(1, length(utr5_sequences))
   
   # Apply Kozak strength rules
   kozak_strengths[context_nucleotide_plus4 == "G"] <- 2
   kozak_strengths[context_nucleotide_minus3 %in% c("A", "G")] <- 3
   kozak_strengths[(context_nucleotide_plus4 == "G") & (context_nucleotide_minus3 %in% c("A", "G"))] <- 4
   
   return(kozak_strengths)
}