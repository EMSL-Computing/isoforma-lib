#' Return an abundance matrix based on summed isotopes 
#' 
#' @description Returns a data.frame of class "abundance_matrix_isoforma" where each
#'    row is a fragment, and each column is a proteoform summed abundance.
#' 
#' @param IonGroup (character) Must be an ion group to subset data down to (a-c or x-z). Required. 
#' @param SummedIsotopes (summed isotopes) The dataframe of summed_isotopes from sum_isotopes. Required.
#' 
#' @return (abundance matrix) A data.frame where each row is a fragment and 
#'    each column is a summed abundance. 
#'    
#' @examples 
#' \dontrun{
#' # Load summed isotopes 
#' SumIso <- readRDS(system.file("extdata", "SummedIsotopes_1to1to1.RDS", package = "isoforma"))
#'
#' # Generate the matrix
#' abundance_matrix(SummedIsotopes = SumIso, IonGroup = "c")
#' }
#' @export
abundance_matrix <- function(SummedIsotopes,
                             IonGroup) {
  
  ##################
  ## CHECK INPUTS ##
  ##################
      
  # Check that SummedIsotopes is a summed_isotopes_isoforma object 
  if (inherits(SummedIsotopes, "summed_isotopes_isoforma") == FALSE) {
    stop("SummedIsotopes must be a summed_isotopes_isoforma object from sum_isotopes.")
  }
  
  # Check the ion group designation
  if (!is.character(IonGroup) || length(IonGroup) > 1 || IonGroup %in% c("a", "b", "c", "x", "y", "z") == FALSE) {
    stop("IonGroup must be a, b, c, x, y, or z listed as a single character.")
  }
  
  #################
  ## MAKE MATRIX ##
  #################
  
  # Update SummedIsotopes
  class(SummedIsotopes) <- c("data.table", "data.frame")
  
  # Remove unmodified sequence, non-relevant ions, get ion position numbers, pivot wider
  AbundanceMatrix <- SummedIsotopes %>% 
    dplyr::filter(Proteoform != "UnmodifiedSequence" & grepl(IonGroup, Ion)) %>%
    dplyr::mutate(
      Ion = gsub("[^[:digit:]., ]", "", Ion) %>% as.numeric()
    ) %>% 
    dplyr::arrange(Proteoform, Ion) %>%
    tidyr::pivot_wider(id_cols = Ion, names_from = Proteoform, values_from = `Summed Intensity`) %>%
    dplyr::arrange(Ion) %>%
    dplyr::mutate(Ion = paste0(IonGroup, Ion)) %>% 
    data.table::data.table()
  
  ####################
  ## RETURN RESULTS ##
  ####################
  
  # Add class
  class(AbundanceMatrix) <- c(class(AbundanceMatrix), "abundance_matrix_isoforma")

  # Return results
  return(AbundanceMatrix)
    
}