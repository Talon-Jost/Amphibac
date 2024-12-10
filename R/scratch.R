

AmphibacMatch_Strict_clustered <- function(df, perc = 0, k = 0, score = 0, 
                                           gap_start = 5, gap = -2, match = 1, 
                                           mismatch = -3, ver = '2023', perc_threshold = 99) {
  
  start <- proc.time()
  
  # Check Biostrings package availability
  if (!requireNamespace("Biostrings", quietly = TRUE)) {
    stop("The Biostrings package is required. Please install it using BiocManager::install('Biostrings').")
  }
  
  library(parallel)
  library(Biostrings)
  library(pbapply)
  
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
  cl <- makeCluster(detectCores() - 1)
  on.exit(stopCluster(cl))  # Ensure the cluster stops when the function exits
  
  # Define the pairwise alignment function
  align_sequences <- function(i, seqs, seqs2, ids, ids2, gap_start, gap, match, mismatch, perc_threshold) {
    for (j in seq_along(seqs2)) {  # Loop over all database sequences
      # Perform pairwise sequence alignment
      alignment <- pairwiseAlignment(
        seqs[i], seqs2[j],
        gapOpening = gap_start,
        gapExtension = gap,
        substitutionMatrix = nucleotideSubstitutionMatrix(match = match, mismatch = mismatch, baseOnly = FALSE)
      )
      
      alignment_score <- score(alignment)
      percent_identity <- pid(alignment, type = "PID1")
      
      # Check if the percent identity meets the threshold
      if (percent_identity >= perc_threshold) {
        return(data.frame(
          Input_ID = ids[i],
          Input_Seq = as.character(seqs[i]),
          DB_ID = ids2[j],
          DB_Seq = as.character(seqs2[j]),
          Score = alignment_score,
          Percent_Identity = percent_identity
        ))
      }
    }
    return(NULL)  # No match found for this sequence
  }
  
  # Perform parallel alignment
  results <- pblapply(seq_along(seqs), function(i) {
    align_sequences(i, seqs, seqs2, ids, ids2, gap_start, gap, match, mismatch, perc_threshold)
  }, cl = cl)
  
  # Combine results into a single data frame
  final_results <- do.call(rbind, results)
  
  end <- proc.time() - start
  message("Alignment completed in ", round(end["elapsed"], 2), " seconds.")
  
  return(final_results)
}

AmphibacMatch_Strict_clustered <- function(df, perc = 0, score = 0, ver = '2023') {
  
  start <- proc.time()
  
  # Check Biostrings package availability
  required_packages <- c("Biostrings", "crayon", "pbapply", "tidyverse", "pryr")
  missing_packages <- required_packages[!sapply(required_packages, requireNamespace, quietly = TRUE)]
  if (length(missing_packages) > 0) {
    stop("The following packages are required but not installed: ", paste(missing_packages, collapse = ", "))
  }
  library(Biostrings)
  library(crayon)
  library(pbapply)
  library(tidyverse)
  library(pryr)
  
  # Ensure the input dataframe has the necessary columns
  if (!all(c("ID", "Sequence") %in% colnames(df))) {
    stop("Input data frame must contain 'ID' and 'Sequence' columns.")
  }
  
  # Set variables
  gap_start = 5 
  gap = -2 
  match = 1 
  mismatch = -3
  k = 1
  perc_threshold = 95
  
  #warning on df size if necessary
  if (nrow(df) >= 100) {
    message("The data frame has more than 100 rows. This can significantly increase computation time")
    message('Do you wish to continue?')
    continue <- readline(prompt = "Type 'yes' to continue or 'no' to stop: ")
    
    if (tolower(continue) != "yes") {
      stop("Process stopped by user.")
    }
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
    alignment <- pairwiseAlignment(
      seqs[i], seqs2[j],
      gapOpening = gap_start,
      gapExtension = gap,
      substitutionMatrix = nucleotideSubstitutionMatrix(match = match, mismatch = mismatch, baseOnly = FALSE)
    )
    alignment_score <- score(alignment)
    percent_identity <- pid(alignment, type = "PID1")
    
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
      return(NULL)
    }
  }
  
  # Generate all pair indices
  pair_indices <- expand.grid(seq_along(seqs), seq_along(seqs2))
  
  # Export variables to the cluster workers
  clusterExport(cl, varlist = c("seqs", "seqs2", "ids", "ids2", "gap_start", "gap", 
                                "match", "mismatch", "align_sequences"),
                envir = environment())
  clusterEvalQ(cl, library(Biostrings))
  
  cat(green("Starting pairwise sequence alignment...\n"))
  
  # Perform pairwise alignment in parallel with progress bar
  results <- pblapply(1:nrow(pair_indices), function(idx) {
    i <- pair_indices[idx, 1]
    j <- pair_indices[idx, 2]
    align_sequences(i, j, seqs, seqs2, ids, ids2, gap_start, gap, match, mismatch, perc_threshold)
  }, cl = cl)
  
  # Stop the cluster
  stopCluster(cl)
  
  # Process matches
  match_list <- Filter(Negate(is.null), results)
  if (length(match_list) == 0) {
    cat(red("No matches found.\n"))
    return(data.frame())
  }
  
  match_df <- do.call(rbind, match_list)
  match_df <- match_df[order(match_df$Percent_Identity, decreasing = TRUE), ]
  
  if (k > 0) {
    match_df <- match_df %>%
      group_by(Input_ID) %>%
      slice_head(n = k) %>%
      ungroup()
  }
  
  # Track system statistics
  time <- proc.time() - start
  memory_used <- as.numeric(pryr::mem_used()) / 1024^2  # Memory in MB
  
  cpu_load <- if (.Platform$OS.type == "unix") {
    as.numeric(system("top -bn1 | grep 'Cpu(s)' | awk '{print $2 + $4}'", intern = TRUE))
  } else {
    as.numeric(system("wmic cpu get loadpercentage", intern = TRUE)[2])
  }
  
  system_info <- list(
    R_version = R.version.string,
    Platform = R.version$platform,
    OS = Sys.info()["sysname"],
    CPU_cores = detectCores(),
    CPU_load = cpu_load,
    Memory_used_MB = memory_used,
    Elapsed_time = time["elapsed"]
  )
  
  # Generate log filename with random hex code
  log_dir <- file.path(system.file("runlog", package = "Amphibac"))
  if (!dir.exists(log_dir)) {
    dir.create(log_dir, recursive = TRUE)  # Create the directory if it doesn't exist
  }
  
  log_filename <- file.path(log_dir, paste0("system_info_", paste(sample(c(0:9, letters), 10, replace = TRUE), collapse = ""), ".log"))
  
  # Save system info to log file
  writeLines(c(
    paste("R Version:", system_info$R_version),
    paste("Platform:", system_info$Platform),
    paste("OS:", system_info$OS),
    paste("CPU Cores:", system_info$CPU_cores),
    paste("CPU Usage (%):", system_info$CPU_usage),
    paste("Memory Used (MB):", round(system_info$Memory_used_MB, 2)),
    paste("Elapsed Time (s):", round(system_info$Elapsed_time, 2))
  ), log_filename)
  
  # Attach attributes
  attr(match_df, "system_info") <- system_info
  attr(match_df, "match_summary") <- list(
    Total_matches = nrow(match_df),
    Avg_percent_identity = mean(match_df$Percent_Identity, na.rm = TRUE)
  )
  attr(match_df, "function_version") <- ver
  
  cat(green("Pairwise sequence alignment completed in", round(time["elapsed"], 2), "seconds.\n"))
  
  return(as.data.frame(match_df))
}


rm(list=ls())
library(Biostrings)
library(crayon)
library(pbapply)
library(tidyverse)
library(ps)
library(parallel)
library(Amphibac)

amp <- AmphiBac_2023.2_Full
amp.sub <- amp %>% 
  top_n(150, ID)
new.end <- AmphibacMatch_Strict_clustered(amp.sub)

#full_fasta <- Amphibac::ReadFasta('mantella_sequences.fasta')
#write.csv(full_fasta, 'Mantella_test_seqs.csv', row.names = F)

#fasta
full_fasta <- Amphibac::ReadFasta('filt_mantella.fasta')

#otu table
maybe <- read.qza('table_mantella.qza')

#metadata
mant.meta <- read.csv('metadata.tsv', sep = '\t')
mant.meta <- mant.meta %>% 
  rename(SampleID = 'X.SampleID')

#this should be the amphibac match df
maybe <- maybe %>% 
  top_n(30) %>% 
  rownames_to_column(var = "ID")

#rownames to filter to
sub_fasta.1 <- full_fasta %>% 
  column_to_rownames(var = 'ID')
sub.maybe <- maybe %>% 
  column_to_rownames(var = 'ID')

fas_id <- rownames(sub_fasta.1)
otu_id <- rownames(sub.maybe)


filtered_table <- sub.maybe[fas_id %in% otu_id,] %>% 
  filter(if_all(everything(), ~ !is.na(.)))
t.ft <- as.data.frame(t(filtered_table)) %>% 
  rownames_to_column(var = 'SampleID')
rar_depth <- 5200
prop <- t.ft %>% 
  mutate(PropAmpMatch = rowSums(across(where(is.numeric))) / rar.depth) %>% 
  select(SampleID, PropAmpMatch)

alt.meta <- tagr.meta %>% 
  left_join(t.ft, by = 'SampleID')



rar_depth = 5200
trial <- filter_func(full_fasta, meta.data = mant.meta, otu.table = maybe, rar.depth = rar_depth)


filter_func <- function(fasta, otu.table, meta.data, rar.depth) {
  
  #check to see if data is in proper order
  if (!"ID" %in% colnames(fasta)) {
    stop("Error in fasta: column 'ID' does not exist")
  }
  if (!"Sequence" %in% colnames(fasta)) {
    stop("Error in fasta: column 'Sequence' does not exist")
  }
  #check if metadata column is the correct name
  if (!"SampleID" %in% colnames(meta.data)) {
    stop("Error in metadata: column 'SampleID' does not exist")
  }
  
  
  #set column names to filter by. This otu table is what will be filtered
  otu_ids <- rownames(otu.table)
  fasta_ids <- fasta %>% 
    column_to_rownames(var = 'ID')
  #this will filter the otu table to only ids that were matched to the database because that's all we care about
  #we will use this to find the proportion of each sample's composition of the matches to the database
  filtered_otu.table <- otu.table[rownames(fasta_ids) %in% otu_ids,] %>%
    filter(if_all(everything(), ~ !is.na(.)))
  
  
  t.ft <- as.data.frame(t(filtered_otu.table)) %>%
    rownames_to_column(var = 'SampleID')
  
  # # #create a dataframe that we can match to the metadata for later for just the proportion
  if (!is.null(rar.depth)) {
    prop <- t.ft %>%
      mutate(PropAmpMatch = rowSums(across(where(is.numeric))) / rar.depth) %>%
      select(SampleID, PropAmpMatch)
  } else {
    prop <- NULL
  }
  
  ampcount.df <- t.ft %>%
    mutate(AmpAbundance = rowSums(across(-SampleID, ~ ifelse(. != 0, 1, 0)))) %>% 
    select(SampleID, AmpAbundance)

  #join to the metadata
  meta2.data <- meta.data %>%
    left_join(prop, by = 'SampleID') %>%
    left_join(ampcount.df, by = 'SampleID')

  #return
  as.data.frame(meta2.data)
  }

rar_depth = 10000
trial <- filter_func(full_fasta, meta.data = tagr.meta, otu.table = maybe, rar.depth = rar_depth)

#fasta
full_fasta1 <- Amphibac::ReadFasta('filt_mantella.fasta')

#otu table
maybe1 <- read.qza('tagr/rarefied_table.qza')
#changing to matches
sub.maybe1 <- maybe1 %>% 
  top_n(50)

#metadata
tagr.meta1 <- read.csv('tagr/na_no_metadata.tsv', sep = '\t')
tagr.meta1 <- tagr.meta1 %>% 
  rename(SampleID = 'X.SampleID')

otu_ids1 <- rownames(maybe1)
fas1_ids1 <- full_fasta1 %>% 
  column_to_rownames(var = 'ID')

filtered_otu.table1 <- sub.maybe1[rownames(fas1_ids1) %in% otu_ids1,] %>% 
  filter(if_all(everything(), ~ !is.na(.)))


filtered_otu.table1 <- sub.maybe1[fas1_ids1 %in% otu_ids1,] %>% 
  filter(if_all(everything(), ~ !is.na(.)))








