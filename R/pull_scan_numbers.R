Glossary <- data.table::fread(
  system.file("extdata", "Unimod_v20220602.csv", package = "pspecterlib")
)

## Define a function for matching peaks in isoforma
.count_matches <- function(PeakData, SeqIsotopes, MS2Window) {
  
  # Change peakdata class
  class(PeakData) <- c("data.table", "data.frame")
  
  # Filter down to the MZ range
  FilterPeaks <- PeakData %>%
    dplyr::filter(`M/Z` >= min(SeqIsotopes$`M/Z`) & `M/Z` <= max(SeqIsotopes$`M/Z`))
  
  if (nrow(FilterPeaks) == 0) {
    return(c("MS1 Matches" = 0, "MS1 Non-Matches" = 0,
             "MS1 Matches in MS2 Window" = 0))
  }
  
  FilterPeaks <- FilterPeaks %>%
    dplyr::rename(`Experimental M/Z` = `M/Z`, `Experimental Intensity` = `Intensity`) %>%
    dplyr::mutate(`Experimental Abundance` = `Experimental Intensity` / max(`Experimental Intensity`) * 100)
  
  # If nothing
  if (nrow(FilterPeaks) == 0) {return(c("Matched" = 0, "Not Matched" = nrow(SeqIsotopes)))}
  
  # Add PPM error 
  SeqIsotopes$`M/Z Tolerance` <- SeqIsotopes$`M/Z` * (5 / 1e6)
  
  # Match closest within 5 ppm
  LeftIndex <- findInterval(SeqIsotopes$`M/Z`, FilterPeaks$`Experimental M/Z`, all.inside = TRUE)
  RightIndex <- LeftIndex + 1
  SeqIsotopes$`Left Difference` <- abs(FilterPeaks$`Experimental M/Z`[LeftIndex] - SeqIsotopes$`M/Z`)
  SeqIsotopes$`Right Difference` <- abs(FilterPeaks$`Experimental M/Z`[RightIndex] - SeqIsotopes$`M/Z`)
  SeqIsotopes$`Closest Index` <- LeftIndex
  
  # Set closest index as right side one, if difference is smaller:
  RightIndexBest <- which(SeqIsotopes$`Right Difference` < SeqIsotopes$`Left Difference`)
  SeqIsotopes$`Closest Index`[RightIndexBest] <- SeqIsotopes$`Closest Index`[RightIndexBest] + 1
  SeqIsotopes$`M/Z Difference` <- abs(FilterPeaks$`Experimental M/Z`[SeqIsotopes$`Closest Index`] - SeqIsotopes$`M/Z`)
  
  # Add PPM error
  SeqIsotopes$`PPM Error` <- SeqIsotopes$`M/Z Difference` / SeqIsotopes$`M/Z` * 1e6
  
  # Count up the matches
  Matches <- SeqIsotopes %>% dplyr::filter(`PPM Error` <= 15)
  Count <- Matches %>% nrow()
  
  # Count up matches in MS2 window
  CountMS2 <- sum(Matches$`M/Z` >= min(MS2Window) & Matches$`M/Z` <= max(MS2Window))
  
  return(c("MS1 Matches" = Count, "MS1 Non-Matches" = nrow(SeqIsotopes) - Count,
           "MS1 Matches in MS2 Window" = CountMS2))
  
}

#' Return a vector of scan numbers based on a Precursor Mass, mass tolerance, and 
#'    isolation window filtering
#' 
#' @param Sequence (character) A valid amino acid sequence. Any non-traditional symbols are
#'     removed and "M." is ignored. Required.
#' @param Modification (character) An IsoForma modifications annotation written as 
#'     "PTM,Residue(Positions)[Number of Modifications]". Required. 
#' @param ScanMetadata (pspecterlib scan metadata object) An object of the scan_metadata class from pspecterlib 
#'     get_scan_metadata. Required.
#' @param MassWindow (numeric) A m/z window for acceptable precursor masses. Default is 5.
#' @param MinMS1Matches (numeric) Minimum number of MS1 peaks from the isotopic distribution
#'     of the full peptide sequence that fall within the MS2 window. Default is 3. 
#' @param ActivationMethod (character) The method used by the MS instrument to dissociate fragments.
#'     Default is ETD.
#' @param AsDataframe (logical) A boolean to indicate whether you would like the values used
#'     to determine the ScanNumbers included or not. Default is FALSE.
#'     
#' @returns A data.frame or list of MS2 scan numbers to sum 
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
#' pull_scan_numbers(Sequence = "M.SGRGKGGKGLGKGGAKRHRKVLRDNIQGITKPAIRRLARRGGVKRISGLIYEETRGVLKIFLENVIRDAVTYTEHARRKTVTAMDVVYALKRQGRTLYGFGG",
#'                   Modification = "Acetyl,X(1*,5,8,12,16,20)[2]",
#'                   ScanMetdata = ScanMetadata,
#'                   PrecursorMass = 11355.38,
#'                   RTStart = 30,
#'                   RTEnd = 60)
#' 
#' }
#' @export
pull_scan_numbers <- function(Sequence,
                              Modification,
                              ScanMetadata,
                              RTStart,
                              RTEnd, 
                              MassWindow = 10,
                              MinMS1Matches = 3,
                              ActivationMethod = "ETD", 
                              AsDataframe = FALSE,
                              MinAbundance = 5) {
  
  ##################
  ## CHECK INPUTS ##
  ##################
  
  # Check that sequence is a character string first
  if (!is.character(Sequence)) {
    stop("Sequence must be a string.")
  }
  
  # Remove M. and any "." from the sequence 
  Sequence <- Sequence %>% 
    gsub(pattern = "M.", replacement = "", fixed = TRUE) %>%
    gsub(pattern = ".", replacement = "", fixed = TRUE)
  
  # Use pspecterlib's sequence check
  if (!pspecterlib::is_sequence(Sequence)) {
    stop("Sequence is not an acceptable peptide/protein sequence. See ?pspecterlib::is_sequence for more details.")
  }
  
  # Check that Modification is a string
  if (!is.character(Modification)) {
    stop("Modification must be a string.")
  }
  
  # Check that the modification is in the glossary
  PTM <- Modification %>% strsplit(",") %>% unlist() %>% head(1)
  if (!(PTM %in% Glossary$Modification)) {
    stop(paste("The modification", PTM, "is not in the backend glossary."))
  }
  
  # Check that ScanMetadata is of the appropriate class
  if ("scan_metadata" %in% class(ScanMetadata) == FALSE) {
    stop("ScanMetadata must be a scan_metadata object from get_scan_metadata in the pspecterlib package.")
  }
  
  # A check to ensure that each single mass value is correct
  .NotSingleNumeric <- function(x) {length(x) != 1 | !is.numeric(x)}
  
  
  # Ensure the MassWindow is a single numeric value 
  if (.NotSingleNumeric(MassWindow)) {
    stop("MassWindow must be a single numeric value.")
  }
  
  # Ensure that both RT start and end are also single numeric values
  if (.NotSingleNumeric(RTStart)) {stop("RTStart must be a single numeric value.")}
  if (.NotSingleNumeric(RTEnd)) {stop("RTEnd must be a single numeric value.")}
  
  ##############################
  ## CALCULATE PRECURSOR MASS ##
  ##############################
  
  # Make a function to get molecular formula 
  get_mol_form <- function(PTM, freq) {
    
    # Extract molecular formula
    MolForm <- Glossary[Glossary$Modification == PTM, 4:ncol(Glossary)] %>%
      dplyr::select(colnames(.)[!is.na(.)]) %>%
      paste0(colnames(.), ., collapse = "") %>%
      pspecterlib::as.molform()
    
    # Multiply molecular formula
    MolForm <- pspecterlib::multiply_molforms(MolForm, freq)
    
    return(MolForm)

  }
  
  # Get modifications and count 
  TotalForm <- NULL
  for (x in unlist(strsplit(Modification, ";"))) {
    Mod <- strsplit(x, ",") %>% unlist() %>% head(1)
    Count <- strsplit(x, "[", fixed = T) %>% unlist() %>% 
        tail(1) %>% gsub(pattern = "]", replacement = "", fixed = T) %>% as.numeric()
    if (is.null(TotalForm)) {TotalForm <- get_mol_form(Mod, Count)} else {
      TotalForm <- pspecterlib::add_molforms(TotalForm, get_mol_form(Mod, Count), CapNegatives = FALSE)
    }
  }
  
  # Now get the total molecular formula 
  MolForm <- pspecterlib::add_molforms(pspecterlib::get_aa_molform(Sequence), TotalForm)
  
  # Get full list of isotoping information for the molecular formula 
  IsotopingData <- pspecterlib::calculate_iso_profile(MolForm, min_abundance = MinAbundance)
  
  # Get the exact mass
  exactmass <- pspecterlib::get_mw(MolForm)
  
  #############################################################################
  ## DETERMINE SCAN NUMBERS WITHIN M/Z AND RT WINDOWS WITH ACTIVATION METHOD ##
  #############################################################################
  
  # Rewrite class of scan metadata 
  class(ScanMetadata) <- c("data.table", "data.frame")
  
  # Pull MS2 scan numbers with a mass and retention time within the designated window
  MS2ScanNumbers <- ScanMetadata %>%
    dplyr::filter(`MS Level` == 2) %>%
    dplyr::mutate(
      "Precursor Charge" = ifelse(`Precursor Charge` == 0, 1, `Precursor Charge`),
      "Precursor Mass" = `Precursor M/Z` * `Precursor Charge` - (1.00727647 * `Precursor Charge`)
    ) %>%
    dplyr::select(c(`Scan Number`, `Precursor Mass`, `Retention Time`, `Activation Method`, 
             `Precursor M/Z`, `Precursor Charge`, `Precursor Scan`)) %>%
    dplyr::filter(`Precursor Mass` >= exactmass - MassWindow & 
           `Precursor Mass` <= exactmass + MassWindow &
           `Retention Time` >= RTStart,
           `Retention Time` <= RTEnd,
           `Activation Method` == "ETD")
  
  # If no MS2ScanNumbers, return error 
  if (nrow(MS2ScanNumbers) == 0) {
    message("No MS2 Scan Numbers detected. Please extend retention time and mass windows.")
    return(NULL)
  }
  
  ##################################
  ## COUNT MATCHES WITHIN WINDOWS ##
  ##################################
  
  # Get header information
  Header <- mzR::header(attr(ScanMetadata, "pspecter")$mzRpwiz)
  
  # Count up matches 
  CountMatches <- do.call(rbind, lapply(1:nrow(MS2ScanNumbers), function(row) {
      
    # Get the MS2 data
    MS2ScanNumber <- MS2ScanNumbers[row, "Scan Number"] %>% unlist()
    PreCharge <- ScanMetadata[ScanMetadata$`Scan Number` == MS2ScanNumber, "Precursor Charge"] %>% unlist()
    SeqIsotopes <- IsotopingData
    colnames(SeqIsotopes)[1] <- c("M/Z")
    
    SeqIsotopes$`M/Z` <- SeqIsotopes$`M/Z` / PreCharge  + 1.007276
    
    # Pull Header data 
    IsolationMZ <- Header %>% dplyr::filter(acquisitionNum == MS2ScanNumber) %>% dplyr::select(isolationWindowTargetMZ) %>% unlist()
    IsolationWindow <- c(IsolationMZ - Header[Header$acquisitionNum == MS2ScanNumber, "isolationWindowLowerOffset"], 
                         IsolationMZ + Header[Header$acquisitionNum == MS2ScanNumber, "isolationWindowUpperOffset"])
    
    # Get MS1 data
    MS1ScanNumber <- MS2ScanNumbers[row, "Precursor Scan"] %>% unlist()
    class(ScanMetadata) <- c(class(ScanMetadata), "scan_metadata")
    MS1Peaks <- pspecterlib::get_peak_data(ScanMetadata, MS1ScanNumber)
    MS1_isomatch <- .count_matches(MS1Peaks, SeqIsotopes, IsolationWindow)
    
    return(
      c(MS2ScanNumber, 
        PreCharge,
        "MS1 Matches" = MS1_isomatch[["MS1 Matches"]],
        "MS1 Non-Matches" = MS1_isomatch[["MS1 Non-Matches"]],
        "MS1 Matches in MS2 Window" = MS1_isomatch[["MS1 Matches in MS2 Window"]]
      )
    )
    
  })) %>% data.table::data.table()
  
  # Add to MS2 scan numbers
  MS2ScanNumbers[,c("MS1 Matches", "MS1 Non-Matches", "MS1 Matches in MS2 Window")] <- CountMatches[,c(3:5)]
  
  # Filter down to proportions
  MS2ScanNumbers <- MS2ScanNumbers %>%
    dplyr::mutate(
      MS1Total = `MS1 Matches` + `MS1 Non-Matches`, 
      MS2Total = max(`MS1 Matches in MS2 Window`),
      Use = `MS1 Matches in MS2 Window` > MinMS1Matches &
        `MS1 Matches` >= floor((MS2ScanNumbers$`MS1 Matches`[1] + MS2ScanNumbers$`MS1 Non-Matches`[1]) * 0.4)
    ) %>%
    dplyr::select(-c(MS1Total, MS2Total))
  
  ####################
  ## RETURN RESULTS ##
  ####################
  
  if (AsDataframe) {
    
    return(list(MS2ScanNumbers, exactmass))
    
  } else {return( MS2ScanNumbers[MS2ScanNumbers$Use, "Scan Number"])}
  
}