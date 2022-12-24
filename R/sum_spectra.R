#' Function to sum spectrum based on scan number
#' 
#' @param ScanMetadata Object of the scan_metadata class from get_scan_metadata 
#'     in pspecterlib. Required.
#' @param ScanNumbers A vector of scan numbers contained in the ScanMetadata 
#'     object. Can be calculated by pull_scan_numbers. Required. 
#' @param MZRound Decimal place the MZ values should be rounded to for binning. 
#'     Default is 3.
#' @param Percentage The percentage of peaks to take. Default is 25%. 
#'  
#' @examples
#' \dontrun{
#' 
#' # Generate the ScanMetadata object 
#' ScanMetadata <- pspecterlib::get_scan_metadata(
#'     MSPath = "/Users/degn400/Desktop/IsoForma_Test/Sorghum-Histone0622162L11.mzML"
#' )
#' 
#' # Run function
#' sum_spectra(ScanMetdata = ScanMetadata,
#'             ScanNumbers = c(1236, 1237, 1240, 1316, 1497, 1520, 1534))
#' }
#' @export
sum_spectra <- function(ScanMetadata, 
                        ScanNumbers,
                        MZRound = 3,
                        Percentage = 25) {
  
  ##################
  ## CHECK INPUTS ##
  ##################
  
  # Check that ScanMetadata is of the appropriate class
  if ("scan_metadata" %in% class(ScanMetadata) == FALSE) {
    stop("ScanMetadata must be a scan_metadata object from get_scan_metadata in the pspecterlib package.")
  }
  
  # Check that ScanNumberVector is numeric
  if (!is.numeric(ScanNumbers)) {
    stop("ScanNumberVector must be numeric.")
  } 
  
  ## Check that RelativeAbundancePercentage is numeric
  if (!is.numeric(MZRound)) {
    stop("RelativeAbundancePercentage must be numeric.")
  }
  RelativeAbundancePercentage <- abs(floor(MZRound))
  
  ##########################
  ## PULL AND SUM SPECTRA ##
  ##########################
  
  # Pull all the spectra - skipping those not in the range
  Spectra <- do.call(rbind, lapply(ScanNumbers, function(scan_number) {
    tryCatch({pspecterlib::get_peak_data(ScanMetadata, scan_number) }, 
             error = function(e) {message(e); return(NULL)})
  })) 
  
  # Stop if no spectra
  if (nrow(Spectra) == 0) {
    stop("No peaks detected. Try different scan numbers within the appropriate range.")
  }
  
  # Bin by the thousands place, rounded, and sum intensities
  Spectra$`M/Z` <- round(Spectra$`M/Z`, MZRound)
  Spectra <- Spectra %>%
    dplyr::group_by(`M/Z`) %>%
    dplyr::summarize(sum(Intensity)) 
  
  # Rename columns
  colnames(Spectra) <- c("M/Z", "Intensity")
  
  # Determine percentage
  Percent <- abs(Percentage) / 100
  if (Percent > 1) {Percent <- 1}
  if (Percent == 0) {Percent <- 1}
  
  # Take top percentage of peaks
  Spectra <- Spectra[order(-Spectra$Intensity)[1:ceiling(nrow(Spectra) * Percent)],] %>%
    dplyr::arrange(`M/Z`)
  
  return(Spectra)
  
}