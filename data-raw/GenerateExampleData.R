# Generate example trees

# Imports
library(dowser)
library(alakazam)

# load example AIRR tsv data
data(ExampleAirr)

ExampleAirr$sequence_alignment = gsub("-","", ExampleAirr$sequence_alignment)

# convert to fasta for use with IgBlast
dowser::dfToFasta(ExampleAirr[1:50,], "test.fasta", seq="sequence_alignment",
    columns=c("biopsy","timepoint"))

ExampleFasta = readFasta("test.fasta")

usethis::use_data(ExampleFasta, overwrite=TRUE)

