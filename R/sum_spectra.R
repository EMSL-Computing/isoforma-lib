#' Function to sum spectrum based on scan number selections, or a list of pspecterlib::peak_data objects
#' 
#' @details Sum spectra based on peak selection from the pspecterlib ScanMetadata object, 
#'    or provide a list of pspecterlib::peak_data objects to sum. 
#' 
#' @param ScanMetadata (scan_metadata object) Object of the scan_metadata class from get_scan_metadata 
#'     in pspecterlib. Required if PeakDataList is NULL.
#' @param ScanNumbers (numeric) Scan numbers contained in the ScanMetadata 
#'     object. Can be calculated by pull_scan_numbers. Required if PeakDataList is NULL.
#' @param PeakDataList (peak_data object) A list of peak_data scan objects. Required if ScanMetadata and
#'     ScanNumbers are NULL.  
#' @param PPMRound (numeric) Number to round to in PPM. Default is 5.
#'  
#' @examples
#' \dontrun{
#' 
#' #################################################
#' ## EXAMPLE 1: PULL FROM A SCAN METADATA OBJECT ##
#' #################################################
#'             
#'             
#' #######################################
#' ## EXAMPLE 2: SUM PEAK DATA TOGETHER ##               
#' #######################################
#' 
#' 
#' 
#' 
#' }
#' @export
sum_spectra <- function(ScanMetadata = NULL, 
                        ScanNumbers = NULL,
                        PeakDataList = NULL,
                        PPMRound = 5) {
  
  ##################
  ## CHECK INPUTS ##
  ##################
  
  # Make sure all 3 parameters are not NULL
  if (is.null(ScanMetadata) | is.null(ScanNumbers)) {
    
    if (is.null(PeakDataList)) {
      stop("Either ScanMetadata and ScanNumbers should be provided, or a PeakDataList.")
    }
    
  }
  
  # Check inputs
  if (is.null(PeakDataList)) {
    
    # Both scan metadata and scan numbers should be provided 
    if (is.null(ScanMetadata) | is.null(ScanNumbers)) {
      stop("Both ScanMetadata and ScanNumbers should not be NULL.")
    }
    
    # Check that ScanMetadata is of the appropriate class
    if ("scan_metadata" %in% class(ScanMetadata) == FALSE) {
      stop("ScanMetadata must be a scan_metadata object from get_scan_metadata in the pspecterlib package.")
    }
    
    # Check that ScanNumberVector is numeric
    if (!is.numeric(ScanNumbers)) {
      stop("ScanNumbers must be numeric.")
    }
    
  } else {
    
    # All PeakDataList objects should be peak_data
    if (all(unlist(lapply(PeakDataList, function(x) {inherits(x, "peak_data")}))) == FALSE) {
      stop("All objects in the PeakDataList must be Peak")
    }
    
  }
  
  # PPM Round should not be a value that will slow down the functions too much
  if (!is.numeric(PPMRound) | PPMRound < 0.1) {
    stop("PPMRound must be a numeric of length 1 greater than 0.1")
  }
  
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
  
  # Function to sum data by the ppm value
  ppm_round <- function(vals) {
    
    lapply(vals, function(x) {
      
      if (x >= 1e7) {stop(paste0("Unexpectedly large number passed to rounding function: ", x))}
      
      # Get the number of decimal points
      deci_num <- (x %>% as.integer() %>% nchar()) - 1
      
      # Now use plyr's round any
      adj <- plyr::round_any(x, PPMRound / 10^6 * 10^deci_num)
      
      return(adj)
      
    }) %>% unlist()
    
  }
  
  # Iterate through MS1_Scans, pull peak data, and round MZ values
  Spectra_Rounded <- do.call(rbind, lapply(Spectra, function(peakData) {
    class(peakData) <- c("data.table", "data.frame")
    peakData <- peakData %>%
      dplyr::mutate(
        `M/Z` = ppm_round(`M/Z`)
      ) %>%
      dplyr::select(`M/Z`, Abundance)
    return(peakData)
  }))
  
  class(Spectra_Rounded) <- c("peak_data", "data.table", "data.frame")
  
  # Set an isoforma peak summing attribute 
  attr(Spectra_Rounded, "isoforma")$PeakSummed <- TRUE
  
  return(Spectra_Rounded)
  
}