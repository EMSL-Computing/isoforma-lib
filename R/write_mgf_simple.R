#' A simple MGF writer 
#' 
#' @description Some tools require MGF files. This function was written to accommodate
#'    those external tools by providing a simple MGF writer.
#' 
#' @param MZ (numeric) Vector of mz values. Must be same length as intensity vector. 
#' @param Intensity (numeric) Vector of intensity values. Must be same length as mz vector. 
#' @param Outpath (character) Path to write file.
#' @param Mass (character) The mass type of the peptide. Default provided.
#' @param Title (character) Name for scan. Default provided. 
#' @param Pepmass (character) The mass of the peptide. Default provided. 
#' @param Charge (character) The charge of the peptide. Default provided. 
#' @param RT (character) The retention time of the peptide. Default provided. 
#' @param Scans (character) The scan number of the peptide. Default provided. 
#' 
#' @returns Writes an mgf file to the specified location in outpath
#' 
#' @examples
#' \dontrun{
#' write_mgf_simple(MZ = c(1,2,3), Intensity = c(1,2,3), Outpath = "~/Downloads/Test.MGF")
#' }
#'     
#' @export
write_mgf_simple <- function(MZ, 
                             Intensity, 
                             Outpath,
                             Mass = "Summed", 
                             Title = "MS/MS gfma zoospore H2B Ac", 
                             Pepmass = "29095.7",
                             Charge = "10+", 
                             RT = "10",
                             Scans = "1") {
  
  ##################
  ## CHECK INPUTS ##
  ##################
  
  # Check mz and intensity vector lengths
  if (length(MZ) != length(Intensity)) {
    stop("MZ and Intensity must be the same length.")
  }
  
  ##################
  ## WRITE OUTPUT ##
  ##################
  
  # Write lines
  cat(paste0("MASS=", Mass, "\nBEGIN IONS", "\nTITLE=", Title, "\nPEPMASS=",
                    Pepmass, "\nCHARGE=", Charge, "\nRTINSECONDS=", RT, "\nSCANS=",
                    Scans), file = Outpath, append = T)
  
  # Iterate through mz and intensity vectors
  for (line in 1:length(MZ)) {
    cat(paste0("\n", MZ[line], "\t", Intensity[line]), file = Outpath, append = T)
  }
  
  # Add last line
  cat("\nEND IONS", file = Outpath, append = T)
  
}