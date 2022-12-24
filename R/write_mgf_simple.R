#' A simple MGF writer 
#' 
#' @param mz Vector of mz values. Must be same length as intensity vector. 
#' @param intensity Vector of intensity values. Must be same length as mz vector. 
#' @param outpath Path to write file.
#' @param mass The mass type of the peptide. Default provided.
#' @param title Name for scan. Default provided. 
#' @param pepmass The mass of the peptide. Default provided. 
#' @param charge The charge of the peptide. Default provided. 
#' @param rt The retention time of the peptide. Default provided. 
#' @param scans The scan number of the peptide. Default provided. 
#'     
#' @export
write_mgf_simple <- function(mz, 
                             intensity, 
                             outpath,
                             mass = "Summed", 
                             title = "MS/MS gfma zoospore H2B Ac", 
                             pepmass = "29095.7",
                             charge = "10+", 
                             rt = "10",
                             scans = "1") {
  
  ##################
  ## CHECK INPUTS ##
  ##################
  
  # Check mz and intensity vector lengths
  if (length(mz) != length(intensity)) {
    stop("mz and intensity must be the same length.")
  }
  
  ##################
  ## WRITE OUTPUT ##
  ##################
  
  # Write lines
  cat(paste0("MASS=", mass, "\nBEGIN IONS", "\nTITLE=", title, "\nPEPMASS=",
                    pepmass, "\nCHARGE=", charge, "\nRTINSECONDS=", rt, "\nSCANS=",
                    scans), file = outpath, append = T)
  
  # Iterate through mz and intensity vectors
  for (line in 1:length(mz)) {
    cat(paste0("\n", mz[line], "\t", intensity[line]), file = outpath, append = T)
  }
  
  # Add last line
  cat("\nEND IONS", file = outpath, append = T)
  
}