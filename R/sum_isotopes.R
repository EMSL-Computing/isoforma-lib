#' Sum isotopes within a matched_peaks object
#' 
#' @details Intensities are divided by their respective charge states and summed together.
#' 
#' @param IsoformaFragments (matched_peaks_isoforma) A matched peaks isoforma object from fragments_per_ptm. Required. 
#' 
#' @returns (summed_isotopes_isoforma) A data.table with 3 columns for the Ion,
#'     Summed Intensity, and the Proteoform. 
#' 
#' @examples
#' \dontrun{
#' 
#' # Load fragment data
#' MultipleFragments <- readRDS(system.file("extdata", "Fragments_1to1to1.RDS", package = "isoforma"))
#' 
#' # Sum isotopes
#' sum_isotopes(MultipleFragments)
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
    
    # Change class to apply dplyr functions
    class(Fragments) <- c("data.frame", "data.table")
    Fragments <- data.table::data.table(Fragments)

    # Select IonType, Intensity, and Z
    SummedIntensities <- Fragments %>%
      dplyr::select(c(Ion, `Intensity Experimental`, Z)) %>%
      dplyr::mutate(
        `Intensity Experimental` = `Intensity Experimental` / Z
      ) %>%
      dplyr::group_by(Ion) %>%
      dplyr::summarise(
        `Summed Intensity` = sum(`Intensity Experimental`),
        Proteoform = names(IsoformaFragments)[el]
      )
    
    # Return Summed Isotopes
    return(SummedIntensities)
    
  })) %>% data.table::data.table()
  
  # Add class and attribute
  class(SummedIsotopes) <- c(class(SummedIsotopes), "summed_isotopes_isoforma")
  
  # Return fragments
  return(SummedIsotopes)
  
}