#' Calculate the fragments per PTM from the modifications_isoforma object
#' 
#' @param Sequences A vector of valid proforma sequences. Required. 
#' @param SummedSpectra A data.table with a "M/Z" and "Intensity" column, preferably
#'     calculated with the function sum_spectra. Required.
#' @param PrecursorCharge A single numeric to determine the Precursor Charge state
#'     (maximum charge) to consider fragments in the calculation. Required. 
#' @param ActivationMethod A string to determine ions. "HCD", "CID", and "ETD" 
#'     select b and y; a, b, and y; and c and z, respectively. Anything else
#'     returns all 6 ions. Required. 
#' @param CorrelationScore A minimum Correlation Score to filter peaks with at 
#'     least 2 isotopes. Default is 0.7. Optional.
#' @param PPMThreshold A minimum PPM Threshold for matching calculated and exprimental
#'     peaks. Default is 10. Optional
#' @param Messages A TRUE/FALSE to indicate whether status messages should be printed.
#'     Default is FALSE.
#' @param ... Additional parameters to pass to pspecter's get_matched_peaks function.
#' 
#' @importFrom foreach %dopar% foreach
#'     
#' @export
fragments_per_ptm <- function(Sequences,
                              SummedSpectra,
                              PrecursorCharge,
                              ActivationMethod,
                              CorrelationScore = 0.7,
                              PPMThreshold = 10, 
                              Messages = FALSE,
                              ...) {
  
  ##################
  ## CHECK INPUTS ##
  ##################
  
  # Check that sequence is a character string first
  if (!is.character(Sequences)) {
    stop("Sequences must be proforma strings.")
  }
  
  # Check that summed spectra has an "M/Z" and "Intensity" in its colnames
  if (!("M/Z" %in% colnames(SummedSpectra) & "Intensity" %in% colnames(SummedSpectra))) {
    stop("SummedSpectra must have both an M/Z and Intensity column. See ?sum_spectra.")
  }
  
  # Check that summed spectra has at least 10 peaks. 
  if (nrow(SummedSpectra) <= 10) {
    stop("SummedSpectra must have at least 10 peaks. There is not enough data to reliably calculate fragments.")
  }
  
  # Check that Precursor Charge is a single number
  if (length(PrecursorCharge) > 1 | !is.numeric(PrecursorCharge)) {
    stop("PrecursorCharge must be a single number.")
  }
  PrecursorCharge <- round(abs(PrecursorCharge))
  
  # Ensure activation method is a single string
  if (length(ActivationMethod) > 1 | !is.character(ActivationMethod)) {
    stop("ActivationMethod should be a single string.")
  }
  
  # Check if Messages is a logical 
  if (!is.logical(Messages)) {
    stop("Messages must be a TRUE/FALSE.")
  }
  
  ############################
  ## CREATE SHORTHAND NAMES ##
  ############################
  
  # Get names of modifications
  ModNames <- lapply(Sequences, function(seq) {
    
    # Determine its name 
    nameDF <- pspecterlib::convert_proforma(seq)
    
    # If it's just the sequence, return original sequence
    if (!is.data.frame(nameDF)) {return("UnmodifiedSequence")}
    
    # Iterate through modifications data.table
    Name <- lapply(1:nrow(nameDF), function (row) {
      Pos <- nameDF$`N Position`[row] %>% unlist()
      paste0(nameDF$Name[row], "@", substr(seq, Pos, Pos), nameDF$`N Position`[row])
    }) %>% paste(collapse = " & ")
    
    # Return the name
    return(Name)
    
  }) %>% unlist()
  
  #########################
  ## CALCULATE FRAGMENTS ##
  #########################
  
  # Message user if applicable
  if (Messages) {message("Starting fragment calculations in parallel")}
  
  # Pull the ions 
  theIons <- .determine_ions(ActivationMethod)
  
  # Generate a peaks dataframe
  thePeaks <- pspecterlib::make_peak_data(MZ = SummedSpectra$`M/Z`,
                                          Intensity = SummedSpectra$Intensity)
  
  # Implement parallel computing for speed 
  doParallel::registerDoParallel(parallel::detectCores()) 
  
  # Calculate all fragments 
  Fragments <- foreach(seq = Sequences) %dopar% {
    
    return(
      pspecterlib::get_matched_peaks(
        PPMThreshold = PPMThreshold, 
        IonGroups = theIons,
        MinimumAbundance = 0.01,
        CorrelationScore = CorrelationScore,
        AlternativeSequence = seq,
        AlternativeSpectrum = thePeaks, 
        AlternativeCharge = PrecursorCharge
      )
    )
    
  }
  
  # Name Fragments
  names(Fragments) <- c(ModNames)
  
  # Add class
  class(Fragments) <- "matched_peaks_isoforma"
  
  # Return object
  return(Fragments)
      
}