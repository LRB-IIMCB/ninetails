% Generated by roxygen2: do not edit by hand
% Please edit documentation in
%   R/ninetails_training_dataset_production_functions.R
\name{prepare_trainingset}
\alias{prepare_trainingset}
\title{Filters out signals of given type of interest.}
\usage{
prepare_trainingset(
  nucleotide,
  nanopolish,
  sequencing_summary,
  workspace,
  num_cores = 1,
  basecall_group = "Basecall_1D_000",
  pass_only = TRUE
)
}
\arguments{
\item{nucleotide}{character. One of the following ["A"/"C"/"G","U"].
This parameter defines the type of filtering approach applied to the input
dataset in order to obtain the signals potentially carrying the
desired nucleotide context.}

\item{nanopolish}{character string. Full path of the .tsv file produced
by nanopolish polya function.}

\item{sequencing_summary}{character string. Full path of the .txt file
with sequencing summary.}

\item{workspace}{character string. Full path of the directory to search the
basecalled fast5 files in. The Fast5 files have to be multi-fast5 file.}

\item{num_cores}{numeric [1]. Number of physical cores to use in processing
the data. Do not exceed 1 less than the number of cores at your disposal.}

\item{basecall_group}{character string ["Basecall_1D_000"]. Name of the
level in the fast5 file hierarchy from which the data should be extracted.}

\item{pass_only}{logical [TRUE/FALSE]. If TRUE, only reads tagged by
nanopolish as "PASS" would be taken into consideration. Otherwise, reads
tagged as "PASS" & "SUFFCLIP" will be taken into account in analysis.
As a default, "TRUE" value is set.}
}
\value{
A list of filtered GAF matrices organized by the read ID_index
is returned. Always assign this returned list to a variable, otherwise
the long list will be printed to the console, which may crash your R session.
}
\description{
Filters out signals corresponding to a given category of reads
(containing A-nucleotides alone or particular types of non-adenosine
nucleotides: G, C or U) in order to prepare the set for neural network
training.
}
\examples{
\dontrun{

prepare_trainingset(nucleotide="A"
                    nanopolish = '/path/to/file',
                    sequencing_summary = '/path/to/file',
                    workspace = '/path/to/guppy/workspace',
                    num_cores = 10,
                    basecall_group = 'Basecall_1D_000',
                    pass_only=TRUE)
}

}
