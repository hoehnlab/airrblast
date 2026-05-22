#' \code{downloadIMGT} Download IMGT GENE-DB databases
#' 
#' Loads all reference germlines from an Immcantation-formatted IMGT database.
#' 
# TODO: make auto-download or internal IMGT database
#' @param dir      directory to contain database
#' @param organism  vector of species to download: human, mouse, rat, rabbit, rhesus_monkey
#' @param timeout  max time allowed for each download (temporarily sets global timeout option)
#' @param verbose  if >0 print out more info
#' @return directory of downloaded IMGT BCR and TCR reference files
#' @details 
#' Recoding of this script: 
#' https://bitbucket.org/kleinstein/immcantation/src/master/scripts/fetch_imgtdb.sh
#' @export
downloadIMGT <- function(dir="imgt", organism=c("human", "mouse", "rat", "rabbit", "rhesus_monkey"),
                         timeout=300, verbose=0){
  
  old_timeout <- getOption("timeout")
  options(timeout = max(timeout, getOption("timeout")))
  
  species_list <- organism
  OUTDIR <- dir
  REPERTOIRE <- "imgt"
  DATE <- format(Sys.time(), "%Y.%m.%d")
  
  # Associative array (for BASH v3) where keys are species folder names and values are query strings
  # rabbit:Oryctolagus_cuniculus, rat:Rattus_norvegicus, rhesus_monkey:Macaca_mulatta
  SPECIES_QUERY=c("human"="Homo+sapiens",
                  "mouse"="Mus",
                  "rat"="Rattus+norvegicus",
                  "rabbit"="Oryctolagus+cuniculus",
                  "rhesus_monkey"="Macaca+mulatta")
  # Associative array (for BASH v3) with species name replacements
  SPECIES_REPLACE=c("human"="Homo sapiens/Homo_sapiens",
                    "mouse"="Mus musculus/Mus_musculus",
                    "rat"="Rattus norvegicus/Rattus_norvegicus",
                    "rabbit"="Oryctolagus cuniculus/Oryctolagus_cuniculus",
                    "rhesus_monkey"="Macaca mulatta/Macaca_mulatta")
  
  for(SPECIES in species_list){
    KEY <- SPECIES
    VALUE <- SPECIES_QUERY[SPECIES]
    REPLACE_VALUE <- SPECIES_REPLACE[SPECIES]
    
    print(paste("Downloading", KEY,"repertoires into",OUTDIR))
    
    # Create directories
    dir.create(file.path(OUTDIR, KEY, "vdj"), recursive=TRUE, showWarnings=FALSE)
    dir.create(file.path(OUTDIR, KEY, "vdj_aa"), recursive=TRUE, showWarnings=FALSE)
    dir.create(file.path(OUTDIR, KEY, "leader_vexon"), recursive=TRUE, showWarnings=FALSE)
    dir.create(file.path(OUTDIR, KEY, "leader"), recursive=TRUE, showWarnings=FALSE)
    dir.create(file.path(OUTDIR, KEY, "constant"), recursive=TRUE, showWarnings=FALSE)
    
    # Download functions
    download_and_process <- function(url, outfile, replace, verbose=0) {
      tmpfile <- paste0(outfile, ".tmp")
      if(verbose > 0){
        print(url)
      }
      utils::download.file(url, tmpfile, quiet=TRUE)
      lines <- readLines(tmpfile, warn=F)
      # fasta between last pre and /pre
      pre_start <- max(grep("<pre>", lines))
      pre_end <- max(grep("</pre>", lines))
      data <- lines[(pre_start + 1):(pre_end - 1)]
      
      # replace names with spaces with unerlines
      split <- strsplit(replace,split="/")[[1]]
      sub <- gsub(split[1], split[2], data)
      writeLines(sub, outfile)
      file.remove(tmpfile)
    }
    
    # Download VDJ regions
    cat("|---- Ig\n")
    for (CHAIN in c("IGHV", "IGHD", "IGHJ", "IGKV", "IGKJ", "IGLV", "IGLJ")) {
      URL <- paste0("https://www.imgt.org/genedb/GENElect?query=7.1+", CHAIN, "&species=", VALUE)
      FILE_NAME <- file.path(OUTDIR, KEY, "vdj", paste0(REPERTOIRE,"_", KEY, "_", CHAIN, ".fasta"))
      download_and_process(URL, FILE_NAME, REPLACE_VALUE, verbose=verbose)
    }
    
    # Download leader and V exon for Ig
    for (CHAIN in c("IGHV", "IGKV", "IGLV")) {
      URL <- paste0("https://www.imgt.org/genedb/GENElect?query=8.1+", CHAIN, "&species=", VALUE,"&IMGTlabel=L-PART1+V-EXON")
      FILE_NAME <- file.path(OUTDIR, KEY, "leader_vexon", paste0(REPERTOIRE,"_lv_", KEY, "_", CHAIN, ".fasta"))
      download_and_process(URL, FILE_NAME, REPLACE_VALUE, verbose=verbose)
    }
    
    # V amino acid for Ig
    for (CHAIN in c("IGHV", "IGKV", "IGLV")) {
      URL <- paste0("https://www.imgt.org/genedb/GENElect?query=7.3+",CHAIN,"&species=",VALUE)
      FILE_NAME <- file.path(OUTDIR, KEY, "vdj_aa", paste0(REPERTOIRE,"_aa_", KEY, "_", CHAIN, ".fasta"))
      download_and_process(URL, FILE_NAME, REPLACE_VALUE, verbose=verbose)
    }
    
    # VDJ TCR
    cat("|---- TCR\n")
    for (CHAIN in c("TRAV", "TRAJ", "TRBV", "TRBD", "TRBJ", "TRDV", "TRDD", "TRDJ", "TRGV", "TRGJ")) {
      URL <- paste0("https://www.imgt.org/genedb/GENElect?query=7.1+",CHAIN,"&species=",VALUE)
      FILE_NAME <- file.path(OUTDIR, KEY, "vdj", paste0(REPERTOIRE,"_", KEY, "_", CHAIN, ".fasta"))
      download_and_process(URL, FILE_NAME, REPLACE_VALUE, verbose=verbose)
    }
    
    # Download leader and V exon for Ig
    for (CHAIN in c("TRAV", "TRBV", "TRDV", "TRGV")) {
      URL <- paste0("https://www.imgt.org/genedb/GENElect?query=8.1+",CHAIN,"&species=",VALUE,"&IMGTlabel=L-PART1+V-EXON")
      FILE_NAME <- file.path(OUTDIR, KEY, "leader_vexon", paste0(REPERTOIRE,"_lv_", KEY, "_", CHAIN, ".fasta"))
      download_and_process(URL, FILE_NAME, REPLACE_VALUE, verbose=verbose)
    }
    
    # V amino acid for TCR
    for (CHAIN in c("TRAV", "TRBV", "TRDV", "TRGV")) {
      URL <- paste0("https://www.imgt.org/genedb/GENElect?query=7.3+",CHAIN,"&species=",VALUE)
      FILE_NAME <- file.path(OUTDIR, KEY, "vdj_aa", paste0(REPERTOIRE,"_aa_", KEY, "_", CHAIN, ".fasta"))
      download_and_process(URL, FILE_NAME, REPLACE_VALUE, verbose=verbose)
    }
    
    
    # Download leaders
    cat("|- Spliced leader regions\n")
    
    # Download leader and V exon for Ig
    cat("|---- Ig\n")
    for (CHAIN in c("IGH", "IGK", "IGL")) {
      URL <- paste0("https://www.imgt.org/genedb/GENElect?query=8.1+",CHAIN,"V&species=",VALUE,"&IMGTlabel=L-PART1+L-PART2")
      FILE_NAME <- file.path(OUTDIR, KEY, "leader", paste0(REPERTOIRE,"_", KEY, "_", CHAIN, "L.fasta"))
      download_and_process(URL, FILE_NAME, REPLACE_VALUE, verbose=verbose)
    }
    
    # Download leader and V exon for TCR
    cat("|---- TCR\n")
    for (CHAIN in c("TRA", "TRB", "TRG", "TRD")) {
      URL <- paste0("https://www.imgt.org/genedb/GENElect?query=8.1+",CHAIN,"V&species=",VALUE,"&IMGTlabel=L-PART1+L-PART2")
      FILE_NAME <- file.path(OUTDIR, KEY, "leader", paste0(REPERTOIRE,"_", KEY, "_", CHAIN, "L.fasta"))
      download_and_process(URL, FILE_NAME, REPLACE_VALUE, verbose=verbose)
    }
    
    # Download constant regions
    cat("|- Spliced constant regions\n")
    cat("|---- Ig\n")
    for (CHAIN in c("IGHC", "IGKC", "IGLC")) {
      # IMGT does not have artificially spliced IGKC / IGLC for multiple species
      if(CHAIN == "IGHC"){
        QUERY <- 14.1
      }else{
        QUERY <- 7.5
      }
      
      URL <- paste0("https://www.imgt.org/genedb/GENElect?query=",QUERY,"+",CHAIN,"&species=",VALUE)
      FILE_NAME <- file.path(OUTDIR, KEY, "constant", paste0(REPERTOIRE,"_", KEY, "_", CHAIN, ".fasta"))
      download_and_process(URL, FILE_NAME, REPLACE_VALUE, verbose=verbose)
    }
    
    cat("|---- TCR\n")
    for (CHAIN in c("TRAC", "TRBC", "TRGC", "TRDC")) {
      URL <- paste0("https://www.imgt.org/genedb/GENElect?query=14.1+",CHAIN,"&species=",VALUE)
      FILE_NAME <- file.path(OUTDIR, KEY, "constant", paste0(REPERTOIRE,"_", KEY, "_", CHAIN, ".fasta"))
      download_and_process(URL, FILE_NAME, REPLACE_VALUE, verbose=verbose)
    }
  }
  
  # Write download info
  INFO_FILE <- file.path(OUTDIR, "IMGT.yaml")
  info <- c(
    "source:  https://www.imgt.org/genedb",
    paste0("date:    ",DATE),
    "species:",
    paste0("    - ",paste0(species_list,":",SPECIES_QUERY[species_list]))
  )
  writeLines(info, con=INFO_FILE)
  
  options(timeout = old_timeout)
  
  return(dir)
}


#' \code{writeFasta} Write a fasta file of sequences given a 
#' named list of sequences
#' @param    seqs      named list of sequences (output from \code{readFasta})
#' @param    file      FASTA file for output
#'
#' @return   File of FASTA formatted sequences
#' @export
writeFasta <- function(seqs, file){
  if(!is.null(seqs)){
    out <- paste0(">", names(seqs), "\n", seqs)
    writeLines(out, con=file)
  }else{
    file.create(file)
  }
}

#' Read a fasta file into a list of sequences
#' Copied straight from Dowser to save a dependency
#' \code{readFasta} reads a fasta file
#' @param    file      FASTA file
#'
#' @return   List of sequences
#' @export
readFasta <- function(file){
  f <- readLines(file)
  if(length(f) == 1){
    return(NULL)
  }
  seqs <- list()
  id <- NA
  for(line in f){
    if(grepl("^>",line)){
      id <- gsub(">","",line)
      seqs[[id]] <- ""
    }else{
      if(is.na(id)){
        stop(paste("Error reading",file))
      }
      seqs[[id]] <- paste0(seqs[[id]],line)
    }
  }
  seqs
}

#' \code{imgtToIgblast} Format IMGT database to Igblast database
#' 
#' Cleans IMGT database and creates blast dbs for IgBlast
#' 
#' @param imgt  directory containing IMGT sequences ("dir" in downloadIMGT)
#' @param outdir   output directory for IgBlast databases, will contain /fasta and /database
#' @param igblast  IgBlast home directory location
#' @param organism   organism species to use (must be human, mouse, or rhesus_monkey)
#' @return outdir
#' @details 
#' Recoding of this script: 
#' https://bitbucket.org/kleinstein/immcantation/src/master/scripts/imgt2igblast.sh
#' @export
imgtToIgblast <- function(imgt, outdir, igblast, 
                          organism=c("human", "mouse", "rhesus_monkey")){
  
  OUTDIR <- outdir
  GERMDIR <- imgt
  OUTFASTA <- file.path(OUTDIR, "fasta")
  OUTDATA <- file.path(OUTDIR, "database")
  
  exec = file.path(igblast, "bin", "makeblastdb")
  exec <- path.expand(exec)
  if(file.access(exec, mode=1) == -1) {
    stop("The file ", exec, " cannot be executed.")
  }
  
  dir.create(file.path(OUTDIR), recursive=TRUE, showWarnings=FALSE)
  dir.create(file.path(OUTFASTA), recursive=TRUE, showWarnings=FALSE)
  dir.create(file.path(OUTDATA), recursive=TRUE, showWarnings=FALSE)
  
  species_list <- dir(GERMDIR)
  species_list <- species_list[species_list != "IMGT.yaml"]
  species_list <- species_list[species_list %in% organism]
  
  
  vdj <- dir(file.path(GERMDIR, species_list, c("vdj")))
  vdj <- vdj[grepl(".fasta$", vdj)]
  vdj_aa <- dir(file.path(GERMDIR, species_list, c("vdj_aa")))
  vdj_aa <- vdj_aa[grepl(".fasta$", vdj_aa)]
  constant <- dir(file.path(GERMDIR, species_list, c("constant")))
  constant <- constant[grepl(".fasta$", constant)]
  
  # Create fasta files of each species, chain and segment combination
  for(SPECIES in species_list){
    #dir(file.path(GERMDIR, SPECIES, "vdj"))
    for(CHAIN in c("IG","TR")){
      # VDJ nucleotides
      for(SEGMENT in c("V","D","J")){
        #F=$(echo imgt_${SPECIES}_${CHAIN}_${SEGMENT}.fasta | tr '[:upper:]' '[:lower:]')
        #cat ${GERMDIR}/${SPECIES}/vdj/imgt_${SPECIES}_${CHAIN}?${SEGMENT}.fasta > ${TMPDIR}/${F}
        # requires updated readFasta
        fastas <- list.files(file.path(GERMDIR, SPECIES, "vdj"))
        fastas <- fastas[grepl(paste0("imgt_", SPECIES, "_", CHAIN, ".", SEGMENT,".*.fasta$"), fastas)]
        f <- unlist(lapply(fastas, function(x)readFasta(file.path(GERMDIR, SPECIES, "vdj", x))))
        cf <- cleanSeqs(f)
        writeFasta(cf, file.path(OUTFASTA, tolower(paste0("imgt_", SPECIES, "_", CHAIN, "_", SEGMENT, ".fasta"))))
      }
      
      # C nucleotides
      #F=$(echo imgt_${SPECIES}_${CHAIN}_c.fasta | tr '[:upper:]' '[:lower:]')
      #cat ${GERMDIR}/${SPECIES}/constant/imgt_${SPECIES}_${CHAIN}?C.fasta > ${TMPDIR}/${F}
      fastas <- list.files(file.path(GERMDIR, SPECIES, "constant"))
      fastas <- fastas[grepl(paste0("imgt_", SPECIES, "_", CHAIN, ".C.*.fasta$"), fastas)]
      f <- unlist(lapply(fastas, function(x)readFasta(file.path(GERMDIR, SPECIES, "constant", x))))
      cf <- cleanSeqs(f)
      writeFasta(cf, file.path(OUTFASTA, tolower(paste0("imgt_", SPECIES, "_", CHAIN, "_C", ".fasta"))))
      
      # V amino acids
      #F=$(echo imgt_aa_${SPECIES}_${CHAIN}_v.fasta | tr '[:upper:]' '[:lower:]')
      #cat ${GERMDIR}/${SPECIES}/vdj_aa/imgt_aa_${SPECIES}_${CHAIN}?V.fasta > ${TMPDIR}/${F}
      fastas <- list.files(file.path(GERMDIR, SPECIES, "vdj_aa"))
      fastas <- fastas[grepl(paste0("imgt_aa_", SPECIES, "_", CHAIN, ".V.*.fasta$"), fastas)]
      f <- unlist(lapply(fastas, function(x)readFasta(file.path(GERMDIR, SPECIES, "vdj_aa", x))))
      cf <- cleanSeqs(f)
      writeFasta(cf, file.path(OUTFASTA, tolower(paste0("imgt_aa_", SPECIES, "_", CHAIN, "_V", ".fasta"))))
    }
  }
  
  fastas <- list.files(OUTFASTA, full.names=FALSE)
  NT_FILES <- fastas[!grepl("imgt_aa_", fastas)]
  AA_FILES <- fastas[grepl("imgt_aa_", fastas)]
  
  for(file in NT_FILES){
    stem <- strsplit(file, split="\\.")[[1]][1]
    command <- paste("-parse_seqids -dbtype nucl -in",
                     file.path(OUTFASTA, file), "-out", file.path(OUTDATA, stem))
    params <- list(exec, command, stdout="", stderr="")
    status <- tryCatch(do.call(base::system2, params), error=function(e){
      return(e)
    }, warning=function(w){
      return(w)
    })
    if(status != 0){
      cat(paste(exec, command, "\n"))
      warning(status, paste(exec, command, "\n"))
    }
  }
  
  for(file in AA_FILES){
    stem <- strsplit(file, split="\\.")[[1]][1]
    command <- paste("-parse_seqids -dbtype prot -in",
                     file.path(OUTFASTA, file), "-out", file.path(OUTDATA, stem))
    params <- list(exec, command, stdout="", stderr="")
    status <- tryCatch(do.call(base::system2, params), error=function(e){
      return(e)
    }, warning=function(w){
      return(w)
    })
    if(status != 0){
      cat(paste(exec, command, "\n"))
      warning(status, paste(exec, command, "\n"))
    }
  }
  return(OUTDIR)
}

# remove IMGT gaps, uppercase sequences, and use second delimiter for name
cleanSeqs <- function(seqs, rm_gaps=TRUE){
  if(length(seqs) == 0){
    return(NULL)
  }
  if(rm_gaps){
    seqs <- toupper(gsub("\\.","",seqs))
  }else{
    seqs <- toupper(seqs)
  }
  names(seqs) <- sapply(strsplit(names(seqs), split="\\|"), function(x)x[2])
  # if repeated IDs, use the first one
  # might be consistent with cleanIMGT.py?
  cseqs <- list()
  for(x in 1:length(seqs)){
    if(is.null(cseqs[[names(seqs)[x]]])){
      cseqs[[names(seqs)[x]]] <- seqs[[x]]
    }
  }
  cseqs
}

#' \code{assignGenes} Runs IgBlast on a fasta file
#'  
#' @param file      Fasta file of Ig or TR sequences
#' @param igblast   Location of IgBlast program directory, containing /bin/igblastn
#' @param refs     Reference directory of sequences (see downloadIMGT)
#' @param igdata  Internal data directory for IgBlast (defaults to \code{igblast}/internal_data)
#' @param organism    Organism (human, mouse, rhesus_monkey)
#' @param domain_system Currently only IMGT supported
#' @param outfile       Name of AIRR rearrangement file (must end in TSV)
#' @param nproc     Nummber of threads to use
#' @param db_prefix   File prefix for reference fastas
#' @param locus     Ig or TR
#' @param set_igdata  Set IGDATA environment variable?
#' @param return    Return data.frame of output?
#' @param verbose   Print extra info?
#' @return AIRR-rearrangement formatted data frame
#' @details 
#' Runs IgBlast, similar to AssignGenes.py in Changeo.
#' Must have IgBlast downloaded, precompiled binaries recommended:
#' https://ncbi.github.io/igblast/
#' 
#' Note for M1/M2 Mac OS, may be necessary to install IgBlast .dmg and 
#' manually set igdata. Otherwise most issues stem from the internal_data
#' directory location.
#' 
#' @export
assignGenes <- function(
    file,
    igblast,
    refs,
    igdata = NULL,
    organism = "human",
    domain_system = "imgt",
    outfile = NULL,
    nproc = 1,
    db_prefix="imgt",
    locus = "Ig",
    set_igdata=TRUE,
    return=TRUE,
    verbose=TRUE){
  
  if(!organism %in% c("human", "mouse", "rhesus_monkey")){
    stop(paste("organism must be either human, mouse, rhesus_monkey"))
  }
  
  if(!locus %in% c("Ig", "TR")){
    stop(paste("locus must be either Ig or TR"))
  }
  
  exec <- file.path(igblast, "bin", "igblastn")
  exec <- path.expand(exec)
  if(file.access(exec, mode=1) == -1) {
    stop("The file ", exec, " cannot be executed.")
  }
  
  if(is.null(igdata)){
    print(paste0("igdata not specified, using ",igblast," for IGDATA"))
    igdata <- igblast
  }
  if(!is.null(igdata) && set_igdata){    
    cat(paste0("Setting IGDATA to ", igdata,"\n"))
    id <- file.path(igdata, "internal_data")
    opt <- file.path(igdata, "optional_file")
    if(!dir.exists(id)){
      stop(paste("Error setting IGDATA, directory:", id, "does not exist. See ?assignGenes for help."))
    }
    if(!dir.exists(opt)){
      stop(paste("Error setting IGDATA, directory:", opt, "does not exist. See ?assignGenes for help."))
    }
    Sys.setenv(IGDATA = igdata)
  }
  
  db_dir <- file.path(refs, "database")
  fasta_dir <- file.path(refs, "fasta")
  if(!dir.exists(db_dir)){
    stop(paste(db_dir, "does not exist."))
  }
  if(!dir.exists(fasta_dir)){
    stop(paste(db_dir, "does not exist."))
  }
  
  if(!grepl("\\.fasta$", file) && !grepl("\\.fa$", file)){
    stop(paste(file, "is not a fasta file (at least, it doesn't end in .fasta or .fa)"))
  }
  
  if(is.null(outfile)){
    filestem <- gsub("\\.fasta$|\\.fa$","", file)
  }else{
    if(!grepl("\\.tsv$", outfile)){
      stop(outfile, " must be a .tsv file")
    }
    filestem <- gsub("\\.tsv", "", outfile)
  }
  
  blastout <- paste0(filestem, ".tsv")
  #airrout <- paste0(filestem, ".tsv")
  
  db_v <- file.path(db_dir, paste(db_prefix, organism, tolower(locus), "v", sep="_"))
  db_d <- file.path(db_dir, paste(db_prefix, organism, tolower(locus), "d", sep="_"))
  db_j <- file.path(db_dir, paste(db_prefix, organism, tolower(locus), "j", sep="_"))
  aux <- file.path(igdata, "optional_file", paste0(organism,"_gl.aux"))
  
  
  stem <- strsplit(file, split="\\.")[[1]][1]
  command <- paste(
    paste0("-germline_db_V ", db_v),
    paste0("-germline_db_D ", db_d),
    paste0("-germline_db_J ", db_j),
    paste0("-auxiliary_data ", aux),
    paste0("-domain_system ", domain_system," -ig_seqtype ",locus," -organism ", organism),
    #paste0("-outfmt '7 std qseq sseq btop'"),
    paste0("-outfmt 19"),
    paste0("-query ", file),
    paste0("-num_threads ", nproc),
    paste0("-out ", blastout)
  )
  if(verbose){
    print(paste(exec, command))
  }
  params <- list(exec, command, stdout="", stderr="")
  status <- tryCatch(do.call(base::system2, params), error=function(e){
    return(e)
  }, warning=function(w){
    return(w)
  })
  if(status != 0){
    cat(paste(exec, command, "\n"))
    stop("error running IgBlast")
  }
  
  if(return){
    return(airr::read_rearrangement(blastout))
  }
}

#' \code{addGaps} Add IMGT gaps to IgBlast output
#'  
#' @param db      AIRR-rearrangment formatted dataframe outputted from \link{assignGenes}
#' @param gapdb   Root directory of gapped fasta files (\code{dir} in \link{downloadIMGT})
#' @param organism  Organism (human, mouse, rhesus_monkey)
#' @param locus   Ig or TR
#' @param gapped_d  Include IMGT gaps in D regions? Only applicable to 
#' @return AIRR-rearrangement formatted data frame in which sequence_alignment and germline_alignment
#' have been updated with IMGT gaps
#' @details 
#' Similar functionality to MakeDb.py in Change-O. 
#' 
#' @export
addGaps <- function(db, gapdb, organism="human", locus="Ig", gapped_d=FALSE){
  
  if(!organism %in% c("human", "mouse", "rhesus_monkey")){
    stop(paste("organism must be either human, mouse, rhesus_monkey"))
  }
  
  if(!locus %in% c("Ig", "TR")){
    stop(paste("locus must be either Ig or TR"))
  }
  
  gap_files <- list.files(file.path(gapdb, organism, "vdj"))
  if(length(gap_files) == 0){
    stop(paste("No files found in",file.path(gapdb, organism, "vdj")))
  }
  vgap_files <- gap_files[grepl(paste0("_",toupper(locus),".V\\.fasta$"), gap_files)]
  gaps <- unlist(lapply(vgap_files, function(x)readFasta(file.path(gapdb, organism, "vdj", x))))
  gaps <- cleanSeqs(gaps, rm_gaps=FALSE)
  
  # IGHD3-10*02 has a gap :-(
  # need to account for that
  dgap_files <- gap_files[grepl(paste0("_",toupper(locus),".D\\.fasta$"), gap_files)]
  dgaps <- unlist(lapply(dgap_files, function(x)readFasta(file.path(gapdb, organism, "vdj", x))))
  dgaps <- cleanSeqs(dgaps, rm_gaps=FALSE)
  
  
  results <- dplyr::tibble()
  for(i in 1:nrow(db)){
    row <- db[i,]
    
    insertion_sites <- 0
    
    if(is.na(row$np1_length) || is.na(row$np1)){
      row$np1_length <- 0
      row$np1 <- ""
    }
    if(is.na(row$np2_length) || is.na(row$np2)){
      row$np2_length <- 0
      row$np2 <- ""
    }
    if(is.na(row$d_sequence_alignment) || is.na(row$d_germline_alignment)){
      row$d_sequence_alignment <- ""
      row$d_germline_alignment <- ""
    }
    if(is.na(row$j_sequence_alignment) || is.na(row$j_germline_alignment)){
      row$j_sequence_alignment <- ""
      row$j_germline_alignment <- ""
    }
    if(is.na(row$v_sequence_alignment)){
      warning(row$sequence_id, " v alignment not found")
      next
    }
    
    if(paste0(row$v_sequence_alignment, row$np1, row$d_sequence_alignment, row$np2, row$j_sequence_alignment) != row$sequence_alignment){
      warning(row$sequence_id, " alignments don't add up")
      next
    }
    
    vcall <- strsplit(row$v_call, split=",")[[1]][1]
    dcall <- strsplit(row$d_call, split=",")[[1]][1]
    jcall <- strsplit(row$j_call, split=",")[[1]][1]
    
    # add IMGT gaps
    # remove germline-relative insertations
    qv <- strsplit(as.character(row$v_sequence_alignment),"")[[1]]
    gv <- strsplit(as.character(row$v_germline_alignment),"")[[1]]
    
    # remove insertions relative to germline
    v_del_pos <- gv == "-"
    gv <- gv[!v_del_pos]
    qv <- qv[!v_del_pos]
    
    qj <- strsplit(as.character(row$j_sequence_alignment),"")[[1]]
    gj <- strsplit(as.character(row$j_germline_alignment),"")[[1]]
    
    # remove insertions relative to germline
    j_del_pos <- gj == "-"
    gj <- gj[!j_del_pos]
    qj <- qj[!j_del_pos]
    
    insertion_sites <- sum(v_del_pos) + sum(j_del_pos)
    
    qd <- strsplit(as.character(row$d_sequence_alignment),"")[[1]]
    gd <- strsplit(as.character(row$d_germline_alignment),"")[[1]]
    
    # remove insertions relative to germline
    d_del_pos <- gd == "-"
    gd <- gd[!d_del_pos]
    qd <- qd[!d_del_pos]
    
    # v gene with IMGT gaps
    gappedv <- strsplit(gaps[[vcall]], split="")[[1]]
    
    # positions without gaps in IMGT
    no_gap_pos <- 1:length(gappedv)
    no_gap_pos <- no_gap_pos[gappedv != "."]
    
    # positions without gaps present in alignment
    v_germ_gap <- rep(".", length(gappedv))
    v_quer_gap <- rep(".", length(gappedv))
    
    # add non gap characters to non-gap sites
    gstart <- row$v_germline_start
    no_gap_pos_g <- no_gap_pos[1:row$v_germline_end] #gapped germline is always full length
    no_gap_pos <- no_gap_pos[gstart:(gstart + length(gv) - 1)]
    
    v_germ_gap[no_gap_pos_g] <- gappedv[no_gap_pos_g]
    v_quer_gap[no_gap_pos] <- qv
    
    if(length(v_germ_gap[no_gap_pos]) != length(gv)){
      stop("Error processing gapped V segment: ", row$sequence_id)
    }
    if(length(v_quer_gap[no_gap_pos]) != length(qv)){
      stop("Error processing gapped V segment: ", row$sequence_id)
    }
    
    # remove trailing gaps
    max_char <- max(which(v_germ_gap != "."))
    v_germ_gap <- v_germ_gap[1:max_char]
    v_quer_gap <- v_quer_gap[1:max_char]
    
    # correct germline coordinate positions
    vgaps_added <- sum(v_quer_gap == ".")
    vpos_del <- sum(v_del_pos)
    jpos_del <- sum(j_del_pos)
    
    row$v_germline_end <- max(no_gap_pos_g)
    # always start at 1 when adding IMGT gaps
    row$v_germline_start <- 1
    
    if(!is.na(dcall) && gapped_d){
      # d gene with IMGT gaps
      gappedd <- strsplit(dgaps[[dcall]], split="")[[1]]
      
      # positions without gaps in IMGT
      no_gap_pos <- 1:length(gappedd)
      no_gap_pos <- no_gap_pos[gappedd != "."]
      
      # positions without gaps present in alignment
      d_germ_gap <- rep(".", length(gappedd))
      d_quer_gap <- rep(".", length(gappedd))
      
      # add non gap characters to non-gap sites
      gdstart <- row$d_germline_start
      no_gap_pos_gd <- no_gap_pos[gdstart:row$d_germline_end] #non-gap positions in germline 
      no_gap_pos <- no_gap_pos[gdstart:row$d_germline_end] #non-gap positions in query
      
      d_germ_gap[no_gap_pos_gd] <- gappedd[no_gap_pos_gd]
      d_quer_gap[no_gap_pos] <- qd
      
      if(length(d_germ_gap[no_gap_pos]) != length(gd)){
        stop("Error processing gapped D segment: ", row$sequence_id)
      }
      if(length(d_quer_gap[no_gap_pos]) != length(qd)){
        stop("Error processing gapped D segment: ", row$sequence_id)
      }
      
      # remove trailing gaps
      max_char <- max(which(d_germ_gap != "."))
      d_germ_gap <- d_germ_gap[1:max_char]
      d_quer_gap <- d_quer_gap[1:max_char]
      
      # remove trailing gaps
      min_char <- min(which(d_germ_gap != "."))
      if(min_char > 1){
        d_germ_gap <- d_germ_gap[(min_char):length(d_germ_gap)]
        d_quer_gap <- d_quer_gap[(min_char):length(d_quer_gap)]
      }
      
      dgaps_added <- sum(d_quer_gap == ".")
      dpos_del <- sum(d_del_pos)
      row$d_germline_end <- max(no_gap_pos_gd) + dgaps_added
      
    }else{
      d_germ_gap <- strsplit(row$d_germline_alignment, split="")[[1]]
      d_quer_gap <- strsplit(row$d_sequence_alignment, split="")[[1]]
      dgaps_added <- 0
      dpos_del <- 0
    }
    
    vdj_quer_gapped <- paste0(
      paste(v_quer_gap,collapse=""),
      row$np1, 
      paste(d_quer_gap,collapse=""),
      row$np2,
      paste(qj,collapse=""))
    
    vdj_germ_gapped <- paste0(
      paste(v_germ_gap,collapse=""),
      paste0(rep("N",row$np1_length),collapse=""),
      paste(d_germ_gap,collapse=""),
      paste0(rep("N",row$np2_length),collapse=""),
      paste(gj,collapse=""))
    
    row$sequence_alignment <- vdj_quer_gapped
    row$germline_alignment <- vdj_germ_gapped
    
    # check region values
    # https://www.imgt.org/IMGTScientificChart/Numbering/IMGTIGVLsuperfamily.html
    if(!is.na(row$fwr1)){
      row$fwr1_start <- 1
      row$fwr1_end <- 26*3
      if(gsub("\\.|-","",substr(row$sequence_alignment, row$fwr1_start, row$fwr1_end)) != row$fwr1){
        if(insertion_sites == 0){
          warning(row$sequence_id, " fwr1 doesn't match")
          #next
        }else{
          warning(row$sequence_id, " fwr1 doesn't match but insertions present")
        }
      }#else{
      row$fwr1 <- substr(row$sequence_alignment, row$fwr1_start, row$fwr1_end)
      #}
    }
    if(!is.na(row$cdr1)){
      row$cdr1_start <- 26*3+1
      row$cdr1_end <- 38*3
      if(gsub("\\.|-","",substr(row$sequence_alignment, row$cdr1_start, row$cdr1_end)) != row$cdr1){
        if(insertion_sites == 0){
          warning(row$sequence_id, " cdr1 doesn't match")
          #next
        }else{
          warning(row$sequence_id, " cdr1 doesn't match but insertions present")
        }
      }#else{
      row$cdr1 <- substr(row$sequence_alignment, row$cdr1_start, row$cdr1_end)
      #}
    }
    if(!is.na(row$fwr2)){
      row$fwr2_start <- 38*3+1
      row$fwr2_end <- 55*3
      if(gsub("\\.|-","",substr(row$sequence_alignment, row$fwr2_start, row$fwr2_end)) != row$fwr2){
        if(insertion_sites == 0){
          warning(row$sequence_id, " fwr2 doesn't match")
          #next
        }else{
          warning(row$sequence_id, " fwr2 doesn't match but insertions present")
        }
      }#else{
      row$fwr2 <- substr(row$sequence_alignment, row$fwr2_start, row$fwr2_end)
      #}
    }
    if(!is.na(row$cdr2)){
      row$cdr2_start <- 55*3+1
      row$cdr2_end <- 65*3
      if(gsub("\\.|\\-","",substr(row$sequence_alignment, row$cdr2_start, row$cdr2_end)) != row$cdr2){
        if(insertion_sites == 0){
          warning(row$sequence_id, " cdr2 doesn't match")
          #next
        }else{
          warning(row$sequence_id, " cdr2 doesn't match but insertions present")
        }
      }#else{
      row$cdr2 <- substr(row$sequence_alignment, row$cdr2_start, row$cdr2_end)
      #}
    }
    if(!is.na(row$fwr3)){
      row$fwr3_start <- 65*3+1
      row$fwr3_end <- 104*3
      if(gsub("\\.|-","",substr(row$sequence_alignment, row$fwr3_start, row$fwr3_end)) != row$fwr3){
        if(insertion_sites == 0){
          warning(row$sequence_id, " fwr3 doesn't match")
          #next
        }else{
          warning(row$sequence_id, " fwr3 doesn't match but insertions present")
        }
      }#else{
      row$fwr3 <- substr(row$sequence_alignment, row$fwr3_start, row$fwr3_end)
      #}
    }
    if(!is.na(row$junction)){
      row$cdr3_start <- 104*3+1
      row$cdr3_end <- 104*3+1 + row$junction_length -6 - 1 + dgaps_added - dpos_del
      row$fwr4_start <- 104*3+1 + row$junction_length -6 - 1 + 1 + dgaps_added - dpos_del
      row$fwr4_end <- nchar(row$sequence_alignment)
      if(gsub("\\.","",substr(row$sequence_alignment, row$cdr3_start, row$cdr3_end)) != row$cdr3){
        if(insertion_sites == 0){
          warning(row$sequence_id, " cdr3 doesn't match")
          #next
        }else{
          warning(row$sequence_id, " cdr3 doesn't match but insertions present")
        }
      }#else{
      row$cdr3 <- gsub("\\.","",substr(row$sequence_alignment, row$cdr3_start, row$cdr3_end))
      #}
      #if(substr(row$sequence_alignment, row$fwr4_start, row$fwr4_end) != row$fwr4){
      # stop(row$sequence_id, " fwr4 doesn't match")
      #}
      row$fwr4 <- gsub("\\.","",substr(row$sequence_alignment, row$fwr4_start, row$fwr4_end))
      if(gsub("\\.","",substr(row$sequence_alignment, 310, 
                              310+row$junction_length-1 + dgaps_added - dpos_del)) != row$junction){
        if(insertion_sites == 0){
          warning(row$sequence_id, " junction doesn't match")
          next
        }else{
          warning(row$sequence_id, " junction doesn't - insertions present")
          next
        }
      }else{
        row$junction <- gsub("\\.","",substr(row$sequence_alignment, 310, 
                                             310+row$junction_length-1 + dgaps_added - dpos_del))
        row$junction_aa <- alakazam::translateDNA(row$junction)
      }
    }
    rm_cols <- c(
      "v_sequence_alignment", "d_sequence_alignment", "j_sequence_alignment",
      "v_germline_alignment", "d_germline_alignment", "j_germline_alignment",
      "v_alignment_start", "d_alignment_start", "j_alignment_start",
      "v_alignment_end", "d_alignment_end", "j_alignment_end",
      "v_sequence_alignment_aa", "d_sequence_alignment_aa", "j_sequence_alignment_aa",
      "v_germline_alignment_aa", "d_germline_alignment_aa", "j_germline_alignment_aa"
    )
    row <- row[,!names(row) %in% rm_cols]
    results <- dplyr::bind_rows(results, row)
  }
  return(results)
}


#' \code{annotate10X} Add 10X annotations to BCR data table
#'  
#' @param db       AIRR-rearrangment formatted dataframe outputted from \link{assignGenes}
#' @param annot    data frame of 10X annotations (read.csv(filtered_contig_annotations.csv))
#' @param data_id  sequence column in db
#' @param annot_id sequence id in annot
#' @param postfix  add postfix to new columns
#' @return AIRR-rearrangement formatted data frame with new columns corresponding to 10X annotations
#' @details 
#' Mapping of 10X annotation names to returned object columns:
#' 'barcode'= 'cell_id',
#' 'c_gene'=  'c_call',
#' 'chain'=   'locus',
#' 'reads'=   'conscount',
#' 'umis'=    'umicount',
#' 'v_gene'=  'v_call_10x',
#' 'd_gene'=  'd_call_10x',
#' 'j_gene'=  'j_call_10x',
#' 'cdr3_nt'= 'junction_10x',
#' 'cdr3'=    'junction_10x_aa'
#' 
#' @export
annotate10X <- function(db, annot, data_id="sequence_id", 
  annot_id="contig_id", postfix="_10x"){

  attr_map <- c(
  'barcode'= 'cell_id',
  'c_gene'=  'c_call',
  'chain'=   'locus',
  'reads'=   'conscount',
  'umis'=    'umicount',
  'v_gene'=  'v_call_10x',
  'd_gene'=  'd_call_10x',
  'j_gene'=  'j_call_10x',
  'cdr3_nt'= 'junction_10x',
  'cdr3'=    'junction_10x_aa')

  if(!data_id %in% names(db)){
    stop(data_id, " not found in data")
  }
  if(!annot_id %in% names(annot)){
    stop(annot_id, " not found in annot")
  }
  m = match(db[[data_id]], annot[[annot_id]])
  if(sum(is.na(m)) > 0){
    if(sum(is.na(m)) == length(m)){
      stop("No matching values found between data columns")
    }else{
      warning(sum(is.na(m)), 
        " annotations couldn't be mapped to data object")
    }
  }
  for(attr in names(attr_map)){
    db[[attr_map[attr]]] <- annot[[attr]][m]
  }

  return(db)
}
