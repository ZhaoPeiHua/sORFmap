Internal.hexamerScore <- function(OneSeq, referFreq, k, step, alphabet){

        count.num <- seqinr::count(OneSeq, wordsize = k, by = step, alphabet = alphabet, freq = FALSE)

        seq.sum  <- 0
        for(i in names(count.num)[count.num != 0]){
                if(referFreq$ref.cds[[i]] != 0 & referFreq$ref.lnc[[i]] != 0) {
                        seq.sum <- seq.sum + (count.num[[i]] * log(referFreq$ref.cds[[i]] / referFreq$ref.lnc[[i]]))
                }
        }
        hex.score <- seq.sum / sum(count.num)

        hex.score
}

compute_hexamerScore <- function(Sequences, label = NULL, referFreq,
                                 k = 6, step = 1, alphabet = c("a", "c", "g", "t"),
                                 on.ORF = FALSE, auto.full = FALSE, parallel.cores = 2) {

        if(on.ORF & !all(alphabet %in% c("a", "t", "g", "c"))) {
                stop("Error: Calculation of hexamer score on ORF region can only be applied to DNA sequences!")
        }

        if (parallel.cores == 2) message("Users can try to set parallel.cores = -1 to use all cores!", "\n")
        message("Processing... ", Sys.time(), "\n")
        parallel.cores <- ifelse(parallel.cores == -1, parallel::detectCores(), parallel.cores)
        cl <- parallel::makeCluster(parallel.cores)
        parallel::clusterExport(cl, varlist = "find_ORF", envir = environment())

        if(on.ORF & !auto.full) {
                message("Calculating ORF region...")
                Sequences <- parallel::parLapply(cl, Sequences, function(x) {
                        orf.info <- find_ORF(x, max.only = TRUE)
                        if(orf.info[[1]] >= 12) orf <- seqinr::s2c(orf.info[[3]]) else orf <- NA
                        orf
                })
                Sequences <- Sequences[!is.na(Sequences)]
        }

        if(on.ORF & auto.full) {
                message("Calculating ORF region...")
                Sequences <- parallel::parLapply(cl, Sequences, function(x) {
                        orf.info <- find_ORF(x, max.only = TRUE)
                        if(orf.info[[1]] >= 12) orf <- seqinr::s2c(orf.info[[3]]) else orf <- x
                        orf
                })
        }

        message("Calculating hexamer score...")
        seqFreq <- parallel::parSapply(cl, Sequences, Internal.hexamerScore, k = k, step = step,
                                       alphabet = alphabet, referFreq = referFreq)
        parallel::stopCluster(cl)

        seqFreq.df <- data.frame(Hexamer.Score = seqFreq)

        if(!is.null(label)) seqFreq.df <- cbind(label = label, seqFreq.df)
        message("\n", "Completed. ", Sys.time(), "\n")

        seqFreq.df
}




find_ORF <- function(OneSeq, max.only = TRUE) {
        OneSeq <- unlist(seqinr::getSequence(OneSeq, as.string = TRUE))
        OneSeq <- gsub("\n", "", OneSeq)
        start_pos <- unlist(gregexpr("atg", OneSeq, ignore.case = TRUE))
        if(sum(start_pos) == -1) {
                if(max.only) {
                        ORF_info <- list(ORF.Max.Len = 0, ORF.Max.Cov = 0, ORF.Max.Seq = NA)
                } else {
                        ORF_info <- data.frame(ORF.Seq = NA, ORF.Start = 0, ORF.Stop = 0, ORF.Len = 0)
                }

        } else {
                stop_pos <- sapply(c("taa", "tag", "tga"), function(x){
                        pos <- unlist(gregexpr(x, OneSeq, ignore.case = TRUE))
                })
                stop_pos <- sort(unlist(stop_pos, use.names = FALSE))
                seq_length <- nchar(OneSeq)
                if(all(stop_pos == -1)) {
                        if(max.only) {
                                orf_start <- min(start_pos)
                                orf_stop  <- seq_length - ((seq_length - orf_start + 1) %% 3) - 2
                                max_len   <- orf_stop - orf_start + 3
                                max_cov   <- max_len / seq_length
                                orf_seq   <- substr(OneSeq, start = orf_start, stop = orf_stop + 2)
                                ORF_info  <- list(ORF.Max.Len = max_len, ORF.Max.Cov = max_cov, ORF.Max.Seq = orf_seq)
                        } else {
                                for(i in seq_along(start_pos)) {
                                        start.val  <- start_pos[i]
                                        stop.val   <- seq_length - ((seq_length - start.val + 1) %% 3) - 2
                                        length.val <- stop.val - start.val + 3
                                        cover.val  <- length.val / seq_length
                                        ORF.seq    <- substr(OneSeq, start = start.val, stop = stop.val + 2)
                                        output.tmp <- data.frame(ORF.Seq = ORF.seq, ORF.Start = start.val, ORF.Stop = stop.val,
                                                                 ORF.Len = length.val, ORF.Cov = cover.val, stringsAsFactors = FALSE)
                                        if(i == 1) {
                                                ORF_info <- output.tmp
                                        } else {
                                                ORF_info <- rbind(ORF_info, output.tmp)
                                        }

                                }
                        }

                } else {
                        orf_starts  <- c()
                        orf_stops   <- c()
                        orf_lengths <- c()
                        for(i in start_pos){
                                orf_flag <- 0
                                for(j in stop_pos){
                                        if(j < i) next
                                        diff_mod    <- (j - i) %% 3
                                        if(diff_mod != 0) next
                                        orf_starts  <- c(orf_starts, i)
                                        orf_stops   <- c(orf_stops, j)
                                        orf_lengths <- c(orf_lengths, j - i + 3)
                                        orf_flag    <- 1
                                        break
                                }
                                if(orf_flag == 0){
                                        orf_starts   <- c(orf_starts, i)
                                        orf_stop_tmp <- seq_length - ((seq_length - i + 1) %% 3) - 2
                                        orf_stops    <- c(orf_stops, orf_stop_tmp)
                                        orf_lengths  <- c(orf_lengths, orf_stop_tmp - i + 3)
                                }
                        }
                        if(max.only) {
                                max_len   <- max(orf_lengths)
                                max_cov   <- max_len / seq_length
                                max_index <- which(orf_lengths == max_len)
                                orf_start <- orf_starts[max_index]
                                orf_stop  <- orf_stops[max_index]
                                orf_seq   <- substr(OneSeq, start = orf_start, stop = orf_stop + 2)
                                ORF_info  <- list(ORF.Max.Len = max_len, ORF.Max.Cov = max_cov, ORF.Max.Seq = orf_seq)
                        } else {
                                for(i in seq_along(orf_starts)) {
                                        start.val  <- orf_starts[i]
                                        stop.val   <- orf_stops[i] + 2
                                        length.val <- orf_lengths[i]
                                        cover.val  <- length.val / seq_length
                                        ORF.seq    <- substr(OneSeq, start = start.val, stop = stop.val)
                                        output.tmp <- data.frame(ORF.Seq = ORF.seq, ORF.Start = start.val, ORF.Stop = stop.val,
                                                                 ORF.Len = length.val, ORF.Cov = cover.val, stringsAsFactors = FALSE)
                                        if(i == 1) {
                                                ORF_info <- output.tmp
                                        } else {
                                                ORF_info <- rbind(ORF_info, output.tmp)
                                        }

                                }
                        }
                }

        }
        ORF_info
}
