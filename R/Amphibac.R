#
#
#' ReadFasta
#' 
#' This function reads a FASTA file and formats it for use with the Amphibac search function.
#' 
#' @param fasta The path to the FASTA file.
#' @return A data frame containing species, isolating strain, and sequences.
#' @export
#this is to load in a fasta file and format it for use with the amphibac search function
ReadFasta <- function(fasta) {
  read_lines <- readLines(fasta)
  original_species <- c()
  isolating_strain <- c()
  sequences <- c()
  headers <- c()
  current_seq <- ""
  
  for (line in read_lines) {
    if (startsWith(line, '>')) {
      if (current_seq != '') {
        sequences <- c(sequences, current_seq)
        current_seq <- ''
      }
      headers <- c(headers, line)
      
      # Split header and handle cases without "-"
      cleaned_header <- gsub('>', '', line)
      split_header <- unlist(strsplit(cleaned_header, "-"))
      
      # Extract original species and isolating strain if available
      original_species <- c(original_species, split_header[1])
      isolating_strain <- c(isolating_strain, ifelse(length(split_header) > 1, split_header[2], NA))
      
    } else {
      current_seq <- paste0(current_seq, line)
    }
  }
  
  if (current_seq != "") {
    sequences <- c(sequences, current_seq)
  }
  
  fasta_df <- data.frame(
    species = original_species,
    isolating_strain = isolating_strain,
    Sequence = sequences,
    stringsAsFactors = FALSE
  )
  fasta_df <- fasta_df[, colSums(!is.na(fasta_df)) > 0]
  return(fasta_df)
}

df_trial <- Amphibac::ReadFasta('data/tagr_trial.fasta')

# new_df <- new_sequence_fasta('R_code/trial_fasta.fasta')
# 
# amphibac_trail <- new_sequence_fasta('R_code/Amphibian-skin_bacteria_16S_sequences.fasta')
# write.csv(amphibac_trail, 'amphibac_df.csv', row.names = F)


# match given dataframe to Amphibac database.
# file structure should include separate ID column and sequence to be called within the function
# perc is percent identity match desired within table output
#' AmphibacMatch
#'
#' This function matches user-provided sequences against a reference database.
#'
#' @param df A data frame containing 'ID' and 'Sequence' columns.
#' @param gap_start The gap opening penalty. Default is -15.
#' @param gap The gap extension penalty. Default is -5.
#' @param match The score for a match. Default is 10.
#' @param mismatch The penalty for a mismatch. Default is -15.
#' @param database The name of the database. Default is 'amphibac'.
#'
#' @return A data frame with alignment results.
#'
#' @import Biostrings
#'
#' @examples
#' user_df <- data.frame(ID = c("seq1", "seq2"), Sequence = c("ATCG", "AAGC"))
#' result <- AmphibacMatch(user_df)
#'
#' @export
#' 
AmphibacMatch <- function(df, gap_start = -15, gap = -5, match = 10, mismatch = -15, database = 'amphibac') {
  if (!requireNamespace("Biostrings", quietly = TRUE)) {
    stop("The Biostrings package is required. Please install it using install.packages('Biostrings').")
  }
  library(Biostrings)
  match_list <- list()
  amphibac_file <- system.file('data', 'amphibac_df.csv', package = 'Amphibac')
  
  if (!file.exists(amphibac_file)) {
    stop("Data file not found in the package.")
  }
  
  amphibac <- read.csv(amphibac_file, sep = ',')
  
  if (!all(c("ID", "Sequence") %in% colnames(df))) {
    stop("Input data frame must contain 'ID' and 'Sequence' columns.")
  }
  
  ids <- df$ID
  seqs <- as.character(df$Sequence)
  seqs <- DNAStringSet(seqs)
  
  for (i in seq_along(seqs)) {
    for (j in seq_along(amphibac$ref_seqs)) {
      seq1 <- seqs[i]
      seq2 <- DNAString(amphibac$ref_seqs[j])
      alignment <- pairwiseAlignment(seq1, seq2,
                                     gapOpening = gap_start,
                                     gapExtension = gap,
                                     substitutionMatrix = nucleotideSubstitutionMatrix(match = match, mismatch = mismatch, baseOnly = TRUE))
      alignment_score <- score(alignment)
      aligned_matches <- nmatch(alignment)
      alignment_length <- width(alignment)
      percent_identity <- (aligned_matches / alignment_length) * 100
      match_list[[length(match_list) + 1]] <- data.frame(ID1 = ids[i], 
                                                         ID2 = amphibac$ref_seqs[j],
                                                         Score = alignment_score,
                                                         PercentIdentity = percent_identity)
    }
  }
  
  match_df <- do.call(rbind, match_list)
  return(match_df)
}

# rm(list = ls())
# library(Amphibac)
trial <- Amphibac::AmphibacMatch(df_trial)

