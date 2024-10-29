# library(devtools)
# library(roxygen2)


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
  
  return(fasta_df)
}


# new_df <- new_sequence_fasta('R_code/trial_fasta.fasta')
# 
# amphibac_trail <- new_sequence_fasta('R_code/Amphibian-skin_bacteria_16S_sequences.fasta')
# write.csv(amphibac_trail, 'amphibac_df.csv', row.names = F)


# match given dataframe to Amphibac database.
# file structure should include separate ID column and sequence to be called within the function
# perc is percent identity match desired within table output
# AmphibacMatch <- function(df, ids, seqs, perc) {
#   df <- c()
#   ids <- c()
#   perc <- c()
# }





