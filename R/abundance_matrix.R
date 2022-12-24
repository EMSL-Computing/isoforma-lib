#' Return an abundance matrix based on summed isotopes 
#' 
#' @param IonGroup Must be an ion group to subset data down to (a-c or x-z). Required. 
#' @param SummedIsotopes The dataframe of summed_isotopes from sum_isotopes. Required.
#' 
#' @export
abundance_matrix <- function(IonGroup, SummedIsotopes, .Threshold = NULL) {
  
  ##################
  ## CHECK INPUTS ##
  ##################
  
  # Check the ion group designation
  if (!is.character(IonGroup) || length(IonGroup) > 1 || IonGroup %in% c("a", "b", "c", "x", "y", "z") == FALSE) {
    stop("IonGroup must be a, b, c, x, y, or z listed as a single character.")
  }
      
  # Check that SummedIsotopes is a summed_isotopes_isoforma object 
  if (class(SummedIsotopes) != "summed_isotopes_isoforma") {
    stop("SummedIsotopes must be a summed_isotopes_isoforma object from sum_isotopes.")
  }
  
  #################
  ## MAKE MATRIX ##
  #################
  
  # Get proteoform names and remove the unmodified 
  Names <- SummedIsotopes$SummedIsotopes$Proteoform %>% unique() %>% .[.!="Unmodified"]

  # Get changing positions
  Positions <- do.call(rbind, lapply(Names, function(Name) {
    Name %>% gsub(pattern = "[^[:digit:]., ]", replacement = "") %>% strsplit(" ") %>% 
      unlist() %>% .[.!=""] %>% as.numeric()
  })) 
  
  # Filter SummedIsotopes down to Names and the appropriate Ion Group, then convert to wide format
  AbundanceMatrix <- SummedIsotopes$SummedIsotopes %>%
    dplyr::filter(Proteoform %in% Names) %>%
    dplyr::filter(grepl(IonGroup, Ion))  %>%
    dplyr::mutate(
      Pos = gsub(pattern = "[^[:digit:]., ]", replacement = "", Ion) %>% as.numeric()
    ) %>%
    dplyr::filter(Pos %in% min(Positions):max(Positions)) %>%
    dplyr::select(-c(Ion, `Peak Count`)) %>%
    stats::reshape(idvar = "Pos", timevar = "Proteoform", direction = "wide") %>%
    dplyr::arrange(Pos) 
  
  if (nrow(AbundanceMatrix) == 0) {stop("No abundance data detected.")}
  
  # Identify missing positions
  MissingPositions <- which(1:20 %in% AbundanceMatrix$Pos == FALSE)
  
  # Add missing positions
  if (length(MissingPositions) > 0) {
    AbundanceMatrix <- dplyr::bind_rows(AbundanceMatrix, 
                                        data.table::data.table(Pos = MissingPositions)) %>% 
      dplyr::arrange(as.numeric(Pos))
  }
  
  # Get ion names and rename columns 
  AbundanceMatrix <- AbundanceMatrix %>%
    dplyr::mutate(Pos = paste0(IonGroup, Pos)) %>%
    dplyr::rename_all(function(x) {c("Ions", 
      lapply(Names, function(x) {
        Names[which(grepl(x, colnames(AbundanceMatrix))) - 1]
      }) %>% unlist()
    )})
  
  # Remove peaks not in threshold
  if (!is.null(.Threshold)) {
    AbundanceMatrix[AbundanceMatrix <= .Threshold] <- NA
  }
  
  ####################
  ## RETURN RESULTS ##
  ####################
  
  # Put the abundance matrix in a list
  AbundanceMatrix <- list(AbundanceMatrix = AbundanceMatrix)
  
  # Add class
  class(AbundanceMatrix) <- "abundance_matrix_isoforma"

  # Return results
  return(AbundanceMatrix)
    
}