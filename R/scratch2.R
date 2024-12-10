rm(list=ls())
library(Biostrings)
library(crayon)
library(pbapply)
library(tidyverse)
library(ps)
library(parallel)
library(Amphibac)
library(rbiom)

full_amp <- function(fasta, ver = '', otu_table = NULL, met = '', rar_depth = 'num', perc = 99) {
  
  #bring in the OTU table
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
  
  
  
  
  otu_table <- read.qza(otu_table)
  
  
  ######1 filter otu table to active occurances, i.e. no zero sums
  otu.table <- as.data.frame(otu_table)
  #####check for non-zero
  otu.table.checked <- otu.table %>% 
    select_if(~ !is.numeric(.) || sum(.) != 0)
  
  ######2 match sequence list to filtered_otu in step 1
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
  
  
  
  fasta.df <- ReadFasta(fasta)
  
  
  
  
  
  
  #######3 run pwalign on matched list from step 2
  AmphibacMatch_Strict_clustered <- function(df, perc = 0, score = 0, ver = '2023') {
    
    #warning on df size if necessary
    if (nrow(df) >= 100) {
      message(red("The data frame has more than 100 rows. This can significantly increase computation time"))
      message('Do you wish to continue?')
      continue <- readline(prompt = "Type 'yes' to continue or 'no' to stop: ")
      
      if (tolower(continue) != "yes") {
        stop("Process stopped by user.")
      }
    }
    start <- proc.time()
    
    # Check required packages
    if (!requireNamespace("Biostrings", quietly = TRUE)) {
      stop("The Biostrings package is required. Please install it using BiocManager::install('Biostrings').")
    }
    if (!requireNamespace("crayon", quietly = TRUE)) {
      stop("The crayon package is required. Please install it using install.packages('crayon').")
    }
    if (!requireNamespace("pbapply", quietly = TRUE)) {
      stop("The pbapply package is required. Please install it using install.packages('pbapply').")
    }
    if (!requireNamespace("tidyverse", quietly = TRUE)) {
      stop("The tidyverse package is required. Please install it using install.packages('tidyverse').")
    }
    if (!requireNamespace("ps", quietly = TRUE)) {
      stop("The ps package is required. Please install it using install.packages('ps').")
    }
    if (!requireNamespace("parallel", quietly = TRUE)) {
      stop("The parallel package is required. Please install it using install.packages('parallel').")
    }
    if (!requireNamespace("rbiom", quietly = TRUE)) {
      stop("The rbiom package is required. Please install it using install.packages('rbiom').")
    }
    
    library(Biostrings)
    library(crayon)
    library(pbapply)
    library(tidyverse)
    library(parallel)
    library(ps)
    library(rbiom)
    
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
      align_sequences(i, j, seqs, seqs2, ids, ids2, gap_start, gap, match, mismatch, 95)
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
    
    if (k == 1) {
      match_df <- match_df %>%
        group_by(Input_ID) %>%
        slice_head(n = k) %>%
        ungroup()
    }
    
    # Calculate runtime, memory, and system information
    time <- proc.time() - start
    memory_used <- pryr::mem_used()  # Total memory used
    
    # Prepare system information
    system_info <- list(
      R_version = R.version.string,
      Platform = R.version$platform,
      OS = Sys.info()["sysname"],
      CPU_cores = detectCores(),
      Memory_used_MB = memory_used / (1024^2),  # Convert to MB
      Elapsed_time = time["elapsed"]
    )
    
    # Log information
    log_dir <- file.path(getwd(), "runlog")
    if (!dir.exists(log_dir)) {
      dir.create(log_dir)
    }
    
    # Create log file name with timestamp (YYYY-MM-DD_HH-MM-SS)
    timestamp <- format(Sys.time(), "%Y-%m-%d_%H-%M-%S")
    log_file <- file.path(log_dir, paste0("runlog_", timestamp, ".log"))
    
    # Open log file for writing
    log_con <- file(log_file, open = "w")
    
    # Write system information and match summary to log file
    cat("Log file created: ", log_file, "\n", file = log_con)
    cat("System Information:\n", file = log_con)
    cat(paste(names(system_info), system_info, sep = ": "), file = log_con, sep = "\n")
    cat("\nMatch Summary:\n", file = log_con)
    cat(paste(names(attr(match_df, "match_summary")), unlist(attr(match_df, "match_summary")), sep = ": "), file = log_con, sep = "\n")
    cat("\nFunction Version: ", ver, "\n", file = log_con)
    
    # Close log file
    close(log_con)
    
    cat(green("Pairwise sequence alignment completed in", round(time["elapsed"], 2), "seconds.\n"))
    
    return(as.data.frame(match_df))
  }
  
  alignment <- AmphibacMatch_Strict_clustered(fasta.df)
  
  new_fasta_df <- alignment %>% 
    select(c(Input_ID, Input_Seq)) %>% 
    rename(ID = 'Input_ID') %>% 
    rename(Sequence = 'Input_seq')
  
  
  
  
  #5 - optional, match to OTU for prop using rar_depth
  filter_function <- function(fasta, otu.table, meta.data, rar.depth) {
    
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
    t.ft <- t.ft %>%
      mutate(PropAmpMatch = rowSums(across(where(is.numeric))) / rarefaction.depth)
    
    # # #create a dataframe that we can match to the metadata for later for just the proportion
    if (!is.null(rar.depth)) {
      prop <- t.ft %>%
        mutate(PropAmpMatch = rowSums(across(where(is.numeric))) / rar.depth)
      prop <- dplyr::select(t.ft, SampleID, PropAmpMatch)
    } else {
      prop <- NULL
    }
    
    ampcount.df <- t.ft %>%
      mutate(AmpAbundance = rowSums(across(-SampleID, ~ ifelse(. != 0, 1, 0))))
    ampcount.df <- dplyr::select(ampcount.df, SampleID, AmpAbundance)
    
    #join to the metadata
    meta2.data <- meta.data %>%
      left_join(prop, by = 'SampleID') %>%
      left_join(ampcount.df, by = 'SampleID')
    
    #return
    as.data.frame(meta2.data) 
  }
  
  
  
  
  
  filter_function(fasta = alignment, otu.table = otu_table, rar.depth = rar_depth, meta.data = met)
}



full_fasta <- Amphibac::ReadFasta('filt_mantella.fasta')

#otu table
maybe <- read.qza('rarefied_table.qza')

#metadata
mant.meta <- read.csv('metadata.txt', sep = '\t')
mant.meta <- mant.meta %>% 
  rename(SampleID = 'X.SampleID')


first_try <- full_amp(fasta = "filt_mantella.fasta", ver = '2023', otu_table = 'rarefied_table.qza', met = mant.meta, rar_depth = 5000)


##graveyard
rar_table <- read.qza('rarefied_table.qza')
rarefaction.depth <- 5000
fasta.orig <- Amphibac::ReadFasta('filt_mantella.fasta')
metadata <- read.csv('metadata.txt', sep = '\t') %>% 
  rename(SampleID = 'X.SampleID')



trial_filt <- filter_function(fasta = fasta.orig, otu.table = rar_table, meta.data = metadata.2, rar.depth = rarefaction.depth)
otu_ids <- rownames(rar_table)
fasta_ids <- fasta.orig %>% 
  column_to_rownames(var = 'ID')
filtered_otu.table <- rar_table[rownames(fasta_ids) %in% otu_ids,] %>%
  filter(if_all(everything(), ~ !is.na(.)))
t.ft <- as.data.frame(t(filtered_otu.table)) %>%
  rownames_to_column(var = 'SampleID')
t.ft <- t.ft %>%
  mutate(PropAmpMatch = rowSums(across(where(is.numeric))) / rarefaction.depth)
prop <- dplyr::select(t.ft, SampleID, PropAmpMatch)
ampcount.df <- t.ft %>%
  mutate(AmpAbundance = rowSums(across(-SampleID, ~ ifelse(. != 0, 1, 0))))
ampcount.df <- dplyr::select(ampcount.df, SampleID, AmpAbundance)

metadata.2 <- subset(metadata, Genus == 'Mantella')
meta2.data <- metadata.2 %>%
  left_join(prop, by = 'SampleID') %>%
  left_join(ampcount.df, by = 'SampleID')
