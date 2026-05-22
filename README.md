
AirrBLAST
-------------------------------------------------------------------------------

Airrblast is an *experimental*, lightweight R package for running AIRR sequence data (BCR/TCR) through
IgBlast. Intended to work within the Immcantation framework as an alternative to AssignGenes.py and 
MakeDb.py in Change-O.

Example usage:

# installation

library(devtools)

install_github("https://github.com/hoehnlab/airrblast")


# Load required packages
library(airrblast)

# igblast download folder
# download and uncompress appropriate IgBlast binaries from
# https://ftp.ncbi.nih.gov/blast/executables/igblast/release/LATEST/
# replace with the name of the directory downloaded
igblast_dir <- "ncbi-igblast-1.22.0"

# Download IMGT GENE-DB int directory "imgt"
# Do this only once
downloadIMGT("imgt", organism="human", verbose=1)

# Set up search databases for IgBlast using imgt database
# Do this only once
imgtToIgblast("imgt", "igblast", igblast_dir)


data(ExampleFasta)

# make example fasta file
writeFasta(ExampleFasta, "test.fasta")

# Align sequences using IgBlast
data <- assignGenes("test.fasta", 
	igblast=igblast_dir,
	refs="igblast",
	outfile="test.tsv",
	nproc=4)

# optional, read in output file using AIRR package
data <- airr::read_rearrangement("test.tsv")

# Add IMGT gaps
gapped_data <- addGaps(data, gapdb="imgt", organism="human", locus="Ig")
gapped_data <- dplyr::filter(gapped_data, !is.na(junction))
