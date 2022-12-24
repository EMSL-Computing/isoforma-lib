#' Sum isotopes within a matched_peaks object
#' 
#' @param IsoformaFragments A matched peaks isoforma object from fragments_per_ptm. Required. 
#' @param Charges The charge values to sum. Default is 1-100. 
#' 
#' @examples
#' \dontrun{
#' 
#' 
#' }
#' @export
sum_isotopes <- function(IsoformaFragments, Charges = 1:100) {
  
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
    
    # Break up into ions,  and sum intensities 
    #SummedIntensities <- Fragments %>%
    #  subset(Z %in% Charges) %>% 
    #  dplyr::select(c(Ion, Z, `Intensity Experimental`, `PPM Error`, Isotope)) %>%
    #  dplyr::mutate(Label = paste0(Ion, "_", Z, "_", Isotope)) %>%
    #  dplyr::group_by(Label) %>%
    #  dplyr::slice(which.min(abs(`PPM Error`))) %>% 
    #  dplyr::select(-c(`PPM Error`, Isotope)) %>%
    #  dplyr::ungroup() %>%
    #  dplyr::mutate(Label = paste0(Ion, "_", Z)) %>%
    #  dplyr::group_by(Label) %>%
    #  tidyr::nest() %>%
    #  dplyr::mutate(
    #    `Summed Intensity` = purrr::map(data, function(x) {sum(x$`Intensity Experimental`)}) %>% unlist(),
    #    `Peak Count` = purrr::map(data, function(x) {length(x$`Intensity Experimental`)}) %>% unlist()
    #  ) %>%
    #  dplyr::mutate(Ion = lapply(Label, function(x) {strsplit(x, "_") %>% unlist() %>% head(1)}) %>% unlist()) %>%
    #  dplyr::ungroup() %>%
    #  dplyr::select(-c(data, Label)) %>%
    #  dplyr::group_by(Ion) %>%
    #  tidyr::nest() %>%
    #  dplyr::mutate(
    #    `Summed Intensity` = purrr::map(data, function(x) {x[which.max(x$`Summed Intensity`), "Summed Intensity"] %>% unlist()}) %>% unlist(),
    #    `Peak Count` = purrr::map(data, function(x) {x[which.max(x$`Summed Intensity`), "Peak Count"] %>% unlist()}) %>% unlist()
    #  ) %>%
    #  dplyr::select(-data) %>%
    #  dplyr::mutate(
    #    Proteoform = names(IsoformaFragments)[el]
    #  )
    
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