calculate_gc_percentage <- function(sequence) {
    # Convert sequence to lowercase to handle both upper and lower case
    sequence_chars <- tolower(strsplit(sequence, '')[[1]])
    
    # Calculate GC percentage
    gc_count <- sum(sequence_chars %in% c('g', 'c'))
    gc_percentage <- gc_count / length(sequence_chars)
    
    return(gc_percentage)
}