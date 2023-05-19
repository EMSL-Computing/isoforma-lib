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
#' @param MinimumAbundance (numeric) The minimum abundance value to keep. Acceptable values
#'    range from 0 to 100. Default is 0.01.
#' @param PPMRound (numeric) Number to round to in PPM. Default is 5.
#' 
#' @returns (peak_data) A peak_data object with summed MS2 spectra
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
#' # Load 3 peaks to sum 
#' PeakDataList <- list(
#'   readRDS(system.file("extdata", "PeakData_1to1to1_1.RDS", package = "isoforma")),
#'   readRDS(system.file("extdata", "PeakData_1to1to1_2.RDS", package = "isoforma")),
#'   readRDS(system.file("extdata", "PeakData_1to1to1_3.RDS", package = "isoforma"))
#' )
#' 
#' sum_ms2_spectra(PeakDataList = PeakDataList)
#' 
#' }
#' @export
sum_ms2_spectra <- function(ScanMetadata = NULL, 
                            ScanNumbers = NULL,
                            PeakDataList = NULL,
                            PPMRound = 5,
                            MinimumAbundance = 0.01) {
  
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
  
  # Check that the minimum abundance filter is reasonable
  if (!is.numeric(MinimumAbundance) || (MinimumAbundance < 0 | MinimumAbundance > 100)) {
    stop("MinimumAbundance should be a single numeric between 0 and 100.")
  }
  
  ##########################
  ## PULL AND SUM SPECTRA ##
  ##########################
  
  if (is.null(PeakDataList)) {
    
    # Pull all the spectra - skipping those not in the range
    Spectra <- do.call(dplyr::bind_rows, lapply(ScanNumbers, function(scan_number) {
      tryCatch({pspecterlib::get_peak_data(ScanMetadata, scan_number) }, 
               error = function(e) {message(e); return(NULL)})
    })) 
    
    # Stop if no spectra
    if (nrow(Spectra) == 0) {
      stop("No peaks detected. Try different scan numbers within the appropriate range.")
    }
    
  } else {
    
    # Fix peak data list classes
    PeakDataList <- lapply(PeakDataList, function(peak_data) {
      class(peak_data) <- c("data.table", "data.frame")
      return(peak_data)
    })
    
    Spectra <- do.call(dplyr::bind_rows, PeakDataList)
    
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
  Spectra_Rounded <- Spectra %>%
    dplyr::mutate(`M/Z` = ppm_round(`M/Z`)) %>%
    dplyr::group_by(`M/Z`) %>%
    dplyr::summarise(
      Intensity = sum(Intensity)
    ) %>%
    dplyr::mutate(ScaleAbundance = Intensity / max(Intensity) * 100) %>%
    dplyr::filter(ScaleAbundance > MinimumAbundance)
  
  # Make Final Spectra
  FinalSpectra <- pspecterlib::make_peak_data(
    MZ = Spectra_Rounded$`M/Z`,
    Intensity = Spectra_Rounded$Intensity
  )
  
  # Set an isoforma peak summing attribute 
  attr(FinalSpectra, "isoforma")$PeakSummed <- TRUE
  
  return(FinalSpectra)
  
}