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
    ID = original_species,
    isolating_strain = isolating_strain,
    Sequence = sequences,
    stringsAsFactors = FALSE
  )
  fasta_df <- fasta_df[, colSums(!is.na(fasta_df)) > 0]
  return(fasta_df)
}

######################################

#
#
#' read.qza
#' 
#' This function reads a .qza file derived from QIIME2 and delivers a editable dataframe in the r environment
#' 
#' @param qza The path to the FASTA file.
#' @return A data frame from a QIIME2 derived relative frequency table
#' 
#' @import rbiom
#' 
#' @export
#this is to load in a fasta file and format it for use with the amphibac search function

read.qza <- function(file) {
  start <- proc.time()
  unzipped_file <- unzip(file)
  biom_path <- file.path(dirname(unzipped_file[1]), 'data')
  biom_ending <- '\\.biom'
  
  find_files_with_ending <- function(folder, file_ending) {
    # Validate the folder path
    if (!dir.exists(folder)) {
      stop("The specified folder does not exist.")
    }
    
    # Search for files with the given ending
    matching_files <- list.files(
      path = folder,
      pattern = paste0(file_ending, "$"), # Regex to match the file ending
      full.names = TRUE                  # Include full file paths
    )
    return(matching_files)
  }
  biom_table <- find_files_with_ending(biom_path, biom_ending)
  otu <- read.biom(biom_table)
  mat <- as.matrix(otu$count)
  mat2 <- as.data.frame(mat)
  final_mat <- mat2 %>%
    select_if(~ !is.numeric(.) || sum(.) != 0)
  #this package decompresses the folder and leaves the folder in the directory. need to figure out how to delete that automatically, or perhaps use a tmp
  unlink(dirname(unzipped_file[1]), recursive = TRUE)
  time <- (proc.time() - start)
  attr(final_mat, "time") <- time["elapsed"]
  return(final_mat)
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
AmphibacMatch <- function(df, perc = 0, k = 0, score = 0, gap_start = 0, gap = 3, 
                          match = 1, mismatch = -3, database = NULL, ver = '2023') {
  start <- proc.time()
  if (!requireNamespace("Biostrings", quietly = TRUE)) {
    stop("The Biostrings package is required. Please install it using install.packages('Biostrings').")
  }
  library(Biostrings)
  match_list <- list()
  
  load("data/AmphiBac_2023.2_Inhibitory.rda")
  
  if (!all(c("ID", "Sequence") %in% colnames(df))) {
    stop("Input data frame must contain 'ID' and 'Sequence' columns.")
  }
  
  #this is just to define what year should be used for the function
  ver_list <- c('2023',
                '2020')
  # Load appropriate database
  if (!is.null(ver)) {
    db_file <- paste0("data/AmphiBac_", ver, ".2_Inhibitory.rda")
    if (!file.exists(db_file)) {
      stop("The specified database version does not exist.")
    }
    load(db_file)  # Assuming 'Amphibac' object is loaded
  } else {
    stop("Please select a database version.")
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
      #percent identity has a few potential methods to calculate using the pid(x, type=) function in pwalign
      percent_identity <- pwalign::pid(alignment, type = 'PID1')
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
  match_df <- subset(match_df, Percent_Identity > perc & Score > score)
  
  if (k > 0) {
    match_df <- match_df %>%
      group_by(df_ids) %>%
      slice_head(n = k) %>%
      ungroup()
  }
  
  time <- (proc.time() - start)
  attr(match_df, "time") <- time["elapsed"]
  as.data.frame(match_df)
}

# match given dataframe to Amphibac database.
# file structure should include separate ID column and sequence to be called within the function
# perc is percent identity match desired within table output
#' AmphibacMatch_Facilitating
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
#' @param ver The year of the database which the user desires to use. Default is the latest version of the database available
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
AmphibacMatch_Facilitating <- function(df, perc = 0, k = 0, score = 0, gap_start = 0, gap = 3, 
                                       match = 1, mismatch = -3, ver = '2023') {
  start <- proc.time()
  # Load required package
  if (!requireNamespace("Biostrings", quietly = TRUE)) {
    stop("The 'Biostrings' package is required. Please install it using BiocManager::install('Biostrings').")
  }
  library(Biostrings)
  
  # Validate input data frame
  if (!all(c("ID", "Sequence") %in% colnames(df))) {
    stop("Input data frame must contain 'ID' and 'Sequence' columns.")
  }
  
  ### Validate version ###
  ver_list <- c('2023', '2020')
  
  #checking to make sure that the version is available if the user DOES define a particular dataset
  #this will default to the latest year, change default year when applicable
  if (!ver %in% ver_list) {
    stop(sprintf("The selected version '%s' is not available. Available versions are: %s.", ver, paste(ver_list, collapse = ", ")))
  }
  
  # Load the appropriate database
  database_path <- sprintf("data/AmphiBac_%s.2_Facilitating.rda", ver)
  if (!file.exists(database_path)) {
    stop(sprintf("Database file '%s' not found.", database_path))
  }
  load(database_path)
  
  # Prepare query and database sequences
  ids <- df$ID
  seqs <- DNAStringSet(as.character(df$Sequence))
  ids2 <- Amphibac$species
  seqs2 <- DNAStringSet(Amphibac$ref_seq)
  
  # Initialize results
  match_list <- list()
  
  # Perform pairwise alignment
  for (i in seq_along(seqs)) {
    for (j in seq_along(seqs2)) {
      alignment <- pairwiseAlignment(
        seqs[i], seqs2[j],
        gapOpening = gap_start,
        gapExtension = gap,
        substitutionMatrix = pwalign::nucleotideSubstitutionMatrix(match = match, mismatch = mismatch, baseOnly = FALSE)
      )
      
      # Compute alignment statistics
      alignment_score <- score(alignment)
      aligned_matches <- nmatch(alignment)
      alignment_length <- nchar(alignment)
      
      #calculate percent identity using pwalign pid
      percent_identity <- pwalign::pid(alignment, type = 'PID1')
      
      # Store results
      match_list[[length(match_list) + 1]] <- data.frame(
        Query_ID = ids[i],
        Query_Sequence = as.character(seqs[i]),
        Reference_ID = ids2[j],
        Reference_Sequence = as.character(seqs2[j]),
        Score = alignment_score,
        Percent_Identity = percent_identity
      )
    }
  }
  
  # Combine results into a data frame
  match_df <- do.call(rbind, match_list)
  match_df <- subset(match_df, Percent_Identity >= perc & Score >= score)
  match_df <- match_df[order(-match_df$Percent_Identity, -match_df$Score), ]
  
  #slice function to limit number of return results according to the number provided by the user
  if (k > 0) {
    match_df <- match_df %>%
      group_by(df_ids) %>%
      slice_head(n = k) %>%
      ungroup()
  }
  time <- (proc.time() - start)
  attr(match_df, "time") <- time["elapsed"]

  return(as.data.frame(match_df))
}

# match given dataframe to Amphibac database.
# file structure should include separate ID column and sequence to be called within the function
# perc is percent identity match desired within table output
#' AmphibacMatch_Strict
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
#' result <- AmphibacMatch_Strict(user_df)
#'
#' @export
#' 
AmphibacMatch_Strict <- function(df, perc = 0, k = 0, score = 0, rel_freq = NULL, gap_start = 0, gap = 3, 
                                 match = 1, mismatch = -3, database = NULL, ver = '2023') {
  start <- proc.time()
  
  # Check Biostrings package availability
  if (!requireNamespace("Biostrings", quietly = TRUE)) {
    stop("The Biostrings package is required. Please install it using install.packages('Biostrings').")
  }
  library(Biostrings)
  
  # Load the Amphibac database
  db_path <- system.file('data', 'AmphiBac_2023.2_StrictInhibitory', package = 'Amphibac')
  if (db_path == "") {
    stop("Database 'AmphiBac_2023.2_StrictInhibitory' not found in the Amphibac package.")
  }
  load(db_path)
  
  if (!all(c("ID", "Sequence") %in% colnames(df))) {
    stop("Input data frame must contain 'ID' and 'Sequence' columns.")
  }
  
  ver_list <- c('2023', '2020')
  if (!ver %in% ver_list) {
    stop(sprintf("Invalid version '%s'. Available versions: %s.", ver, paste(ver_list, collapse = ", ")))
  }
  
  ids <- df$ID
  seqs <- DNAStringSet(as.character(df$Sequence))
  ids2 <- Amphibac$species
  seqs2 <- DNAStringSet(Amphibac$ref_seq)
  
  match_list <- list()
  for (i in seq_along(seqs)) {
    for (j in seq_along(seqs2)) {
      alignment <- Biostrings::pairwiseAlignment(
        seqs[i], seqs2[j],
        gapOpening = gap_start,
        gapExtension = gap,
        substitutionMatrix = pwalign::nucleotideSubstitutionMatrix(match = match, mismatch = mismatch, baseOnly = FALSE)
      )
      alignment_score <- score(alignment)
      percent_identity <- pid(alignment, type = "PID1")
      
      match_list[[length(match_list) + 1]] <- data.frame(
        df_ids = ids[i],
        seq = as.character(seqs[i]),
        Amphibac_ids = ids2[j],
        ref_seq = as.character(seqs2[j]),
        Score = alignment_score,
        Percent_Identity = percent_identity
      )
      
      match_list <- match_list[sapply(match_list, function(x) x$Score >= k & x$Percent_Identity >= perc)]
    }
  }

  match_df <- do.call(rbind, match_list)
  match_df <- subset(match_df, Percent_Identity > perc & Score > score)
  match_df <- match_df[order(match_df$Percent_Identity, decreasing = TRUE), ]
  
  if (k > 0) {
    match_df <- match_df %>%
      group_by(df_ids) %>%
      slice_head(n = k) %>%
      ungroup()
  }
  
  amphibac.relab <- function(match_df, rel_freq_table) {
    rel_w_ids <- rel_freq_table %>% 
      rownames_to_column(var = 'OTU_ID') %>% 
      mutate(row_sum = rowSums(select(., -OTU_ID)))
    
    new_df <- rel_w_ids %>% 
      filter(OTU_ID %in% match_df$df_ids) %>%
      mutate(rel_sum = rowSums(select(., -c(OTU_ID, row_sum)))) %>%
      mutate(RelAbun = rel_sum / row_sum) %>%
      select(-row_sum) %>%
      column_to_rownames(var = 'OTU_ID')
    
    return(new_df)
  }
  
  if (!is.null(rel_freq)) {
    new_df <- amphibac.relab(match_df, rel_freq)
    time <- proc.time() - start
    attr(new_df, "time") <- time["elapsed"]
    return(new_df)
  }
  
  time <- proc.time() - start
  attr(match_df, "time") <- time["elapsed"]
  return(as.data.frame(match_df))
}

# match given dataframe to Amphibac database.
# file structure should include separate ID column and sequence to be called within the function
# perc is percent identity match desired within table output
#' AmphibacMatch_Strict_clustered
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
#' @import crayon
#' @import pbapply
#' @import tidyverse
#'
#' @examples
#' user_df <- data.frame(ID = c("seq1", "seq2"), Sequence = c("ATCG", "AAGC"))
#' result <- AmphibacMatch_Strict_clustered(user_df)
#'
#' @export
#' 
AmphibacMatch_Strict_clustered <- function(df, perc = 0, k = 0, score = 0, 
                                           gap_start = 5, gap = -2, match = 1, 
                                           mismatch = -3, ver = '2023', perc_threshold = 95) {
  
  start <- proc.time()
  
  # Check Biostrings package availability
  if (!requireNamespace("Biostrings", quietly = TRUE)) {
    stop("The Biostrings package is required. Please install it using install.packages('Biostrings').")
  }
  
  library(Biostrings)
  library(crayon)
  library(pbapply)
  library(tidyverse)
  
  # Ensure the input dataframe has the necessary columns
  if (!all(c("ID", "Sequence") %in% colnames(df))) {
    stop("Input data frame must contain 'ID' and 'Sequence' columns.")
  }
  
  # Convert sequences to DNAStringSet
  ids <- df$ID
  seqs <- DNAStringSet(as.character(df$Sequence))
  ids2 <- Amphibac$species
  seqs2 <- DNAStringSet(Amphibac$ref_seq)
  
  # Create a cluster with the desired number of cores
  cl <- makeCluster(detectCores() - 1)  # Use all but 1 core
  
  # Define the pairwise alignment function
  align_sequences <- function(i, j, seqs, seqs2, ids, ids2, gap_start, gap, match, mismatch, perc_threshold) {
    # Perform pairwise sequence alignment
    alignment <- pairwiseAlignment(
      seqs[i], seqs2[j],
      gapOpening = gap_start,
      gapExtension = gap,
      substitutionMatrix = nucleotideSubstitutionMatrix(match = match, mismatch = mismatch, baseOnly = FALSE)
    )
    
    alignment_score <- score(alignment)
    percent_identity <- pid(alignment, type = "PID1")
    
    # Only proceed if the percent identity is above the threshold
    if (percent_identity >= perc_threshold) {
      return(data.frame(
        Input_ID = ids[i],
        Input_Seq = as.character(seqs[i]),
        DB_ID = ids2[j],
        DB_Seq = as.character(seqs2[j]),
        Score = alignment_score,
        Percent_Identity = percent_identity
      ))
    } else {
      return(NULL)  # Skip this pair if percent identity is too low
    }
  }
  
  # Generate all pair indices (combinations of sequences to align)
  pair_indices <- expand.grid(seq_along(seqs), seq_along(seqs2))
  
  # Export variables and functions to the cluster workers
  clusterExport(cl, varlist = c("seqs", "seqs2", "ids", "ids2", "gap_start", "gap", 
                                "match", "mismatch", "align_sequences", "perc_threshold"),
                envir = environment())
  clusterEvalQ(cl, library(Biostrings))
  
  
  cat("Starting pairwise sequence alignment...\n")
  
  # Perform pairwise alignment in parallel
  results <- pblapply(1:nrow(pair_indices), function(idx) {
    i <- pair_indices[idx, 1]
    j <- pair_indices[idx, 2]
    if (idx %% 100 == 0) {  # Print every 100th iteration
      cat("Processed", idx, "pairs\n")
    }
    align_sequences(i, j, seqs, seqs2, ids, ids2, gap_start, gap, match, mismatch, perc_threshold)
  })
  
  # Stop the cluster after computation
  stopCluster(cl)
  
  match_list <- Filter(Negate(is.null), results)
  match_df <- do.call(rbind, match_list)
  match_df <- as.data.frame(match_df)
  match_df$Percent_Identity <- as.numeric(match_df$Percent_Identity)  # Ensure it's numeric
  match_df <- match_df[order(match_df$Percent_Identity, decreasing = TRUE), ]
  
  if (k > 0) {
    match_df <- match_df %>%
      group_by(Input_ID) %>%
      slice_head(n = k) %>%
      ungroup()
  } else {
    stop("There are no results which match to the Amphibac database based on input parameters.")
  }
  
  time <- proc.time() - start
  attr(match_df, "time") <- time["elapsed"]
  
  cat(green("Pairwise sequence alignment completed.\n"))
  
  return(as.data.frame(match_df))
}
############################


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




