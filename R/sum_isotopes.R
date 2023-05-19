#' Sum isotopes within a matched_peaks object
#' 
#' @details Intensities are divided by their respective charge states and summed together.
#' 
#' @param IsoformaFragments A matched peaks isoforma object from fragments_per_ptm. Required. 
#' 
#' @examples
#' \dontrun{
#' 
#' 
#' }
#' @export
sum_isotopes <- function(IsoformaFragments) {
  
  ##################
  ## CHECK INPUTS ##
  ##################
  
  # Check matched_peaks_isoforma object
  if (class(IsoformaFragments) != "matched_peaks_isoforma") {
    stop("object must be a matched_peaks_isoforma object from fragments_per_ptm.")
  }
  
  # Check that charges is a numeric
  if (length(Charges) == 0 | !is.numeric(Charges)) {
    stop("Charges must be numeric.")
  }
  Charges <- abs(floor(Charges))
  
  ##################
  ## SUM ISOTOPES ##
  ##################
  
  # Sum up MatchedPeaks isotopes, add count of isotopes
  SummedIsotopes <- do.call(rbind, lapply(1:length(IsoformaFragments), function(el) {

    # Pull out fragments data 
    Fragments <- IsoformaFragments[[el]]
    
    # If there's NA's, let the user know and remove
    if (is.na(Fragments$Ion) %>% any()) {
      warning("NAs ions detected in matched_peaks. This should not happen. They were removed.")
      Fragments <- Fragments[is.na(Fragments$Ion) == FALSE,]
    }
    
    browser()
    
    SummedIntensities <- Fragments %>%
      dplyr::group_by(Ion) %>%
      tidyr::nest() %>%
      dplyr::mutate(
        `Summed Intensity` = purrr::map(data, function(x) {sum(x$`Intensity Experimental`)}) %>% unlist(),
        `Peak Count` = NA,
        Proteoform = names(IsoformaFragments)[el]
      ) %>%
      dplyr::select(-data)
    
    
    # Return Summed Isotopes
    return(SummedIntensities)
    
  })) %>% data.table::data.table()
  
  # Make the object
  SummedIsotopes <- list("SummedIsotopes" = SummedIsotopes)
  
  # Add class and attribute
  class(SummedIsotopes) <- "summed_isotopes_isoforma"
  
  # Return fragments
  return(SummedIsotopes)
  
}