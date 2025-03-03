rm(list = ls())

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


vsearch_amphibac <- function(cmd = "vsearch", sequences, amphibac, output_txt, output_csv) {
  args <- c(
    "--usearch_global", sequences,
    "--db", amphibac, 
    "--blast6out", output_txt,
    "--id", "1.0",
    "--strand", "both"
  )
  
  # Run the command
  result <- system2(cmd, args, stdout = TRUE, stderr = TRUE)
  
  # Check if the output file was created and read it
  if (file.exists(output_txt)) {
    blast_results <- read.delim(output_txt, header = FALSE)
    
    # Save to CSV
    write.csv(blast_results, file = output_csv, row.names = FALSE)
    
    message("Results saved to: ", output_csv)
    return(output_csv)
  } else {
    warning("VSEARCH did not produce an output file.")
    return(NULL)
  }
}

# Example usage
csv_path <- vsearch_amphibac("vsearch", "input_sequences.fasta", "amphibac_db.fasta", 
                             "new", "user_output.csv")
# Example usage
output <- vsearch_amphibac("vsearch", "dna-sequences.fasta", "edit_mantella_refdb.fasta", "new1.csv")
as.data.frame(output)


results <- read.csv('new1.csv', sep = '\t', header = F)
colnames(results) <- c('ID')

rartable <- read_qza('rarefied_table.qza')$data
rartable <- as.data.frame(t(rartable))

rartable$Reads <- rowSums(rartable, na.rm = T)

match.table <- rartable[rownames(rartable) %in% results$ID,]
match.table <- as.data.frame(match.table)
match.table2 <- as.data.frame(t(match.table))
match.table2$sum <- rowSums(match.table2, na.rm = T)
match.table2 <- match.table2 %>% 
  mutate(PropAmp = sum / 5000)

seq.fasta <- ReadFasta('dna-sequences.fasta')
fixed.fasta <- seq.fasta %>% 
  rename(ID = 'species')
mantella.fasta <- fixed.fasta[fixed.fasta$ID %in% rownames(match.table),]

write.fasta(sequences = as.list(mantella.fasta$Sequence),
            name = mantella.fasta$ID,
            file.out = 'vsearch_results.fasta')

write.csv(results, "AmphiBac.matches.csv")


##

# Example usage
mantella.subset.vsearch <- vsearch_amphibac("vsearch", "rar.mantella.equences.fasta", "AmphibBac_FullDatabase_2023.2.fasta", "mantella.amphibac.matches.csv")
mant.df <- read.csv('mantella.amphibac.matches.csv', header = F, sep = '\t')
library(qiime2R)
otu.table <- read_qza('mantella.rarefied.table.qza')$data
filt.otu <- otu.table[rownames(otu.table) %in% mant.df$V1,]
filt.otu <- as.data.frame(filt.otu)
filt.otu <- as.data.frame(t(filt.otu))  # Ensure proper transposition
filt.otu[] <- lapply(filt.otu, as.numeric)  # Convert all columns to numeric
filt.otu$PropAmphiBac <- rowSums(filt.otu) / 5000  # Compute rarefaction metric

