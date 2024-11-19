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




# match given dataframe to Amphibac database.
# file structure should include separate ID column and sequence to be called within the function
# perc is percent identity match desired within table output
#' AmphibacMatch
#'
#' This function matches user-provided sequences against a reference database.
#'
#' @param df A data frame containing 'ID' and 'Sequence' columns.
#' @param perc numerical percentage cutoff for minimum percent identity match in results df. Default is 0.
#' @param score a numerical cutoff for alignment score in results df. Default is 0.
#' @param gap_start The gap opening penalty. Default is 0.
#' @param gap The gap extension penalty. Default is 3.
#' @param match The score for a match. Default is 1.
#' @param mismatch The penalty for a mismatch. Default is -3.
#' @param database The name of the database to be used to compare df sequences against. Default is 'amphibac'.
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
AmphibacMatch <- function(df, perc = 0, score = 0, gap_start = 0, gap = 3, 
                          match = 1, mismatch = -3, database = 'amphibac') {
  if (!requireNamespace("Biostrings", quietly = TRUE)) {
    stop("The Biostrings package is required. Please install it using install.packages('Biostrings').")
  }
  library(Biostrings)
  match_list <- list()
  load("data/Amphibac.rda")
  if (!all(c("ID", "Sequence") %in% colnames(df))) {
    stop("Input data frame must contain 'ID' and 'Sequence' columns.")
  }
  
  ids <- df$ID
  seqs <- as.character(df$Sequence)
  seqs <- DNAStringSet(seqs)
  ids2 <- Amphibac$species
  seqs2 <- DNAStringSet(Amphibac$ref_seq)
  
  for (i in seq_along(seqs)) {
    for (j in seq_along(seqs2)) {
      seq1 <- seqs[i]
      seq2 <- seqs2[j]
      alignment <- pwalign::pairwiseAlignment(seq1, seq2,
                                              gapOpening = gap_start,
                                              gapExtension = gap,
                                              substitutionMatrix = pwalign::nucleotideSubstitutionMatrix(match = match, mismatch = mismatch, baseOnly = FALSE))
      alignment_score <- score(alignment)
      aligned_matches <- nmatch(alignment)
      alignment_length <- nchar(alignment)
      percent_identity <- (aligned_matches / alignment_length) * 100
      match_list[[length(match_list) + 1]] <- data.frame(df_ids = ids[i],
                                                         seq = seqs[i],
                                                         Amphibac_ids = ids2[j],
                                                         ref_seq = seqs2[j],
                                                         Score = alignment_score,
                                                         Percent_Identity = percent_identity)
    }
  }
  
  match_df <- do.call(rbind, match_list)
  match_df <- match_df[order(match_df$Percent_Identity, decreasing = TRUE), , drop = FALSE]
  # match_df <- subset(match_df, Percent_Identity > perc & Score > score)
  
  as.data.frame(match_df)
}


#' Amphibac Dataset
#'
#' Dataframe containing anti-batrachochytrium dendrobatidis sequences and the organism they were isolated from
#'
#' @format a data frame with 1944 rows and 3 variables:
#' \describe{
#'   \item{species}{the species that the bacterial sequence was originally isolated from}
#'   \item{isolating_strain}{the name of the strain og the sequence}
#'   \item{ref_seq}{the sequence from the region of identification}
#'   ...
#' }
#' @source \url{http://somewhere.important.com/}
#' @docType data
#' @name Amphibac
#' @usage data(Amphibac)
#' @usage a data frame with 1944 rows and 3 variables
"Amphibac"



#' AmphiBac_2023.2_Facilitating Dataset
#'
#' Dataframe containing facilitory sequences and the organism they were isolated from
#'
#' @format a data frame with 597 rows and 2 variables:
#' \describe{
#'   \item{species}{the species that the bacterial sequence was originally isolated from}
#'   \item{isolating_strain}{the name of the strain og the sequence}
#'   \item{ref_seq}{the sequence from the region of identification}
#'   ...
#' }
#' @source \url{http://somewhere.important.com/}
#' @docType data
#' @name AmphiBac_2023.2_Facilitating
#' @usage data(AmphiBac_2023.2_Facilitating)
#' @usage a data frame with 597 rows and 2 variables
"AmphiBac_2023.2_Facilitating"




#' Amphibac Dataset
#'
#' Dataframe containing anti-batrachochytrium dendrobatidis sequences and the organism they were isolated from
#'
#' @format a data frame with 1944 rows and 3 variables:
#' \describe{
#'   \item{species}{the species that the bacterial sequence was originally isolated from}
#'   \item{isolating_strain}{the name of the strain og the sequence}
#'   \item{ref_seq}{the sequence from the region of identification}
#'   ...
#' }
#' @source \url{http://somewhere.important.com/}
#' @docType data
#' @name Amphibac
#' @usage data(Amphibac)
#' @usage a data frame with 1944 rows and 3 variables
"Amphibac"


#' AmphiBac_2023.2_Full Dataset
#'
#' Dataframe containing all sequences contained within the AmphiBac database and the organism they were isolated from
#'
#' @format a data frame with 7927 rows and 2 variables:
#' \describe{
#'   \item{species}{the species that the bacterial sequence was originally isolated from}
#'   \item{isolating_strain}{the name of the strain og the sequence}
#'   \item{ref_seq}{the sequence from the region of identification}
#'   ...
#' }
#' @source \url{http://somewhere.important.com/}
#' @docType data
#' @name AmphiBac_2023.2_Full
#' @usage data(AmphiBac_2023.2_Full)
#' @usage a data frame with 7927 rows and 2 variables
"AmphiBac_2023.2_Full"

######################

#' AmphiBac_2023.2_Inhibitory Dataset
#'
#' Dataframe containing all sequences contained within the AmphiBac database and the organism they were isolated from
#'
#' @format a data frame with 2518 rows and 2 variables:
#' \describe{
#'   \item{species}{the species that the bacterial sequence was originally isolated from}
#'   \item{isolating_strain}{the name of the strain og the sequence}
#'   \item{ref_seq}{the sequence from the region of identification}
#'   ...
#' }
#' @source \url{http://somewhere.important.com/}
#' @docType data
#' @name AmphiBac_2023.2_Inhibitory
#' @usage data(AmphiBac_2023.2_Inhibitory)
#' @usage a data frame with 2518 rows and 2 variables
"AmphiBac_2023.2_Inhibitory"

########################


#' AmphiBac_2023.2_StrictInhibitory Dataset
#'
#' Dataframe containing all sequences contained within the AmphiBac database and the organism they were isolated from
#'
#' @format a data frame with 2056 rows and 2 variables:
#' \describe{
#'   \item{species}{the species that the bacterial sequence was originally isolated from}
#'   \item{isolating_strain}{the name of the strain og the sequence}
#'   \item{ref_seq}{the sequence from the region of identification}
#'   ...
#' }
#' @source \url{http://somewhere.important.com/}
#' @docType data
#' @name AmphiBac_2023.2_StrictInhibitory
#' @usage data(AmphiBac_2023.2_StrictInhibitory)
#' @usage a data frame with 2056 rows and 2 variables
"AmphiBac_2023.2_StrictInhibitory"




