\name{ReadFasta}
\alias{ReadFasta}
\title{ReadFasta}
\description{Load in a fasta/fna file which contains separated sample ID and sequence by line. This will output as a callable dataframe within R}

\usage{
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
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fasta}{
%%     ~~Describe \code{fasta} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or standard data sets, see data().

## The function is currently defined as
function (fasta) 
{
    read_lines <- readLines(fasta)
    original_species <- c()
    isolating_strain <- c()
    sequences <- c()
    headers <- c()
    current_seq <- ""
    for (line in read_lines) {
        if (startsWith(line, ">")) {
            if (current_seq != "") {
                sequences <- c(sequences, current_seq)
                current_seq <- ""
            }
            headers <- c(headers, line)
            cleaned_header <- gsub(">", "", line)
            split_header <- unlist(strsplit(cleaned_header, "-"))
            original_species <- c(original_species, split_header[1])
            isolating_strain <- c(isolating_strain, ifelse(length(split_header) > 
                1, split_header[2], NA))
        }
        else {
            current_seq <- paste0(current_seq, line)
        }
    }
    if (current_seq != "") {
        sequences <- c(sequences, current_seq)
    }
    fasta_df <- data.frame(species = original_species, isolating_strain = isolating_strain, 
        Sequence = sequences, stringsAsFactors = FALSE)
    return(fasta_df)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
