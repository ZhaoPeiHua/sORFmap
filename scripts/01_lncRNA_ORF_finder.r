library(seqinr)
library(foreach)
library(doParallel)
library(stringr)

# Read and filter gene annotations
gene_annotations <- read.table("../data/input/v37.genetype.txt")
gene_annotations <- gene_annotations[gene_annotations$V7=="lncRNA",]

# Read and filter lncRNA sequences
lncrna_sequences <- read.fasta("../data/input/gencode.v37.transcripts.fa")
gene_ids <- str_extract(names(lncrna_sequences), "ENSG[0-9]+\\.[0-9]+")
lncrna_sequences <- lncrna_sequences[gene_ids %in% gene_annotations$V5]

total_sequences <- length(lncrna_sequences)
registerDoParallel(60)

all_orfs <- foreach(sequence_index = 1:total_sequences, .combine=rbind, .packages = "seqinr") %dopar% {
    orf_results <- c()
    orf_entry <- c()
    
    current_sequence <- lncrna_sequences[sequence_index]
    sequence_data <- unlist(current_sequence)
    dna_sequence <- paste(unname(sequence_data), collapse="")
    dna_sequence <- toupper(dna_sequence)
    
    sequence_name <- names(current_sequence)
    sequence_name <- strsplit(sequence_name,"\\|")
    sequence_name <- paste(sequence_name[[1]][c(1,2,6,7,8)], collapse="_")
    print(sequence_index)
    
    start_positions <- unlist(gregexpr('ATG', dna_sequence))
    stop_positions <- unlist(gregexpr('TAA|TGA|TAG', dna_sequence))
    
    if (any(start_positions != -1) && any(stop_positions != -1)) {   
        for (start_pos in start_positions) {
            for (stop_pos in stop_positions) {
                orf_length <- stop_pos - start_pos
                if ((orf_length %% 3 == 0) & (start_pos < stop_pos)) {
                    orf_id <- paste(c(sequence_name, start_pos, stop_pos+2), collapse="_")
                    coding_sequence <- substr(dna_sequence, start_pos, stop_pos+2)
                    utr5_sequence <- substr(dna_sequence, 1, start_pos-1)
                    utr3_sequence <- substr(dna_sequence, stop_pos+3, nchar(dna_sequence))
                    amino_acids <- translate(s2c(coding_sequence))
                    amino_acids <- paste(amino_acids, collapse="")
                    amino_acids <- gsub("\\*$", "", amino_acids)
                    
                    if (nchar(amino_acids) >= 10) {
                        orf_entry <- data.frame(orf_id, utr5_sequence, coding_sequence, utr3_sequence, amino_acids)
                        orf_results <- rbind(orf_results, orf_entry)
                    }
                    break
                }
            }
        }
    }
    return(orf_results)
}

write.table(all_orfs, "../data/intermediate/v37.LncRNA.translated.ORFs.txt", sep="\t", quote=F, row.names=F)
write.fasta(sequences = as.list(all_orfs$amino_acids), names = all_orfs$orf_id, file.out = "../data/intermediate/v37.LncRNA.translated.peptides.fasta")
