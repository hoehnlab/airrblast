##' @importFrom utils packageDescription
.onAttach <- function(libname, pkgname) {
  air_ref <- paste('Jensen CG, Sumner JA, Kleinstein SH, Hoehn KB (2024).',
  ' Inferring B cell phylogenies from paired H and L chain BCR sequences with Dowser.',
  'Journal of Immunology. https://doi.org/10.4049/jimmunol.2300851.')
  ref <- paste0("If you are using airrblast in published research please cite ", air_ref)
  if(!is.null(ref)) packageStartupMessage(ref)
}
