
AirrBLAST
-------------------------------------------------------------------------------

Airrblast is an *experimental*, lightweight R package for running AIRR sequence data (BCR/TCR) through
IgBlast. Intended to work within the Immcantation framework as an alternative to AssignGenes.py and 
MakeDb.py in Change-O.

Example usage:

### installation
``` r
library(devtools)

install_github("https://github.com/hoehnlab/airrblast")
```


### IgBLAST setup
download and uncompress appropriate IgBlast binaries from
https://ftp.ncbi.nih.gov/blast/executables/igblast/release/LATEST/
replace with the name of the directory downloaded


### Set up IMGT directory and IgBlast files
```r
library(airrblast)

igblast_dir <- "ncbi-igblast-1.22.0"

# Download IMGT GENE-DB int directory "imgt"
# Do this only once
downloadIMGT("imgt", organism="human", verbose=1)

# Set up search databases for IgBlast using imgt database
# Do this only once
imgtToIgblast("imgt", "igblast", igblast_dir)
```

### Run IgBLAST on example data
```r
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
```

### Identifying clones, reconstructing germlines, building trees.

From here, you can perform any analyses supported by Immcantation within R. Here, we perform clonal clustering, reconstruct clonal germlines, build trees, and visualize them. These topics are handled in more detail on other pages. In particular:

[Clonal clustering with Scoper](https://scoper.readthedocs.io/en/stable/vignettes/Scoper-Vignette/)
[Building clonal germlines with Dowser](Germlines-Vignette.md)
[Building trees with Dowser](Building-Trees-Vignette.md)
[Plotting trees with Dowser](Plotting-Trees-Vignette.md)

For a full tour, check out the [Immcantation tutorials](https://immcantation.readthedocs.io/en/stable/getting_started/getting-started.html). 

```r
# note, these next steps use Scoper and Dowser functions
library(scoper)
library(dowser)

# Identify clones using a threshold of 0.1 for convenience
gapped_data <- hierarchicalClones(gapped_data, 0.1, summarize=FALSE)

# read in IMGT database for clonal germline reconstruction
# this is the same directory downloaded above. Isn't that nice?
imgt = readIMGT("imgt/human/vdj")

# reconstruct clonal germlines
gapped_data = createGermlines(gapped_data, imgt, trim_lengths=TRUE)

# formatClones for tree building
clones = formatClones(gapped_data)

# build trees
trees = getTrees(clones, build="pml")

# plot trees
plots <- plotTrees(trees)

# plot largest tree
plots[[1]]

# save all trees as a pdf
treesToPDF(plots, file="trees.pdf")

# save tree objects for later use
saveRDS(trees, "trees.rds")
```
