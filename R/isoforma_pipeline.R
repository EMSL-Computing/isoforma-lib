#' Calculates fragments, sums isotopes and charge states, generates the abundance matrix,
#'     and calculates the proportions of each proteoform
#'
#' @description A wrapper function to run fragments_per_ptm, sum_isotopes, abundance_matrix,
#'    and calculate_proportions.
#'
#' @param Sequences (character) A vector of valid ProForma sequences. Required.
#' @param SummedSpectra (peak_data) A data.table with a "M/Z" and "Intensity" column,
#'     preferably calculated with the function sum_spectra. Required.
#' @param MaxFragmentCharge (numeric) A single numeric to determine the Precursor
#'     Charge state (maximum charge) to consider fragments in the calculation. Required.
#' @param ActivationMethod (character) A string to determine ions. "HCD", "CID", and
#'     "ETD" select b and y; a, b, and y; and c and z, respectively. Anything else
#'     returns all 6 ions. Required.
#' @param IonGroup (character) Must be an ion group to subset data down to (a-c or x-z). Required.
#' @param CorrelationScore (numeric) A minimum Correlation Score to filter peaks with
#'     at least 2 isotopes. Default is 0.7. Optional.
#' @param PPMThreshold (numeric) A minimum PPM Threshold for matching calculated
#'     and exprimental peaks. Default is 10. Optional.
#' @param IsotopeAlgorithm (character) "isopat" uses the isopat package to calculate
#'     isotopes, while "Rdisop" uses the Rdisop package. Though more accurate,
#'     Rdisop has been known to crash on Windows computers when called iteratively
#'     more than 1000 times. Default is Rdisop, though isopat is an alternative.
#' @param MinimumAbundance (numeric) Minimum abundance for calculating isotopes. Default is 0.01.
#' @param Message (logic) A TRUE/FALSE to indicate whether status messages should be printed. Default is FALSE.
#' @param Top	(integer) The top N proportions to return, ranked from higest to lowest. Default is 8.
#'
#' @returns (list) A list with 5 objects: matched peaks isoforma, summed isotopes isoforma,
#'     abundance matrix isoforma, the proportions dataframe, and a ggplot of the proportions
#'
#' @examples
#' \dontrun{
#' 
#' ###################################################
#' ## BRUNNER VALINE DATASET - KNOWN CONCENTRATIONS ##
#' ###################################################
#' 
#' # Load 3 peaks to sum
#' PeakDataList <- list(
#'  readRDS(system.file("extdata", "PeakData_1to1to1_1.RDS", package = "isoforma")),
#'  readRDS(system.file("extdata", "PeakData_1to1to1_2.RDS", package = "isoforma")),
#'  readRDS(system.file("extdata", "PeakData_1to1to1_3.RDS", package = "isoforma"))
#' )
#'
#' # Sum peaks
#' PeakSum <- sum_ms2_spectra(PeakDataList = PeakDataList)
#'
#' # Get proteoform strings
#' MultipleMods <- pspecterlib::multiple_modifications(
#'     Sequence = "LQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG",
#'     Modification = "6.018427,V(17,26,70)[1]",
#'     ReturnUnmodified = TRUE
#' )
#'
#' # Run main pipeline function
#' IsoForma <- isoforma_pipeline(
#'     Sequences = MultipleMods,
#'     SummedSpectra = PeakSum,
#'     MaxFragmentCharge = 11,
#'     ActivationMethod = "ETD",
#'     IonGroup = "c",
#'     CorrelationScore = 0,
#'     IsotopeAlgorithm = "isopat", # Rdisop is preferred, is faster, and is more accurate, but it tends to crash on Windows
#'     Messages = TRUE
#' )
#'
#' ################################################
#' ## PASAVENTO DATASET - UNKNOWN CONCENTRATIONS ##
#' ################################################
#' 
#' # Load raw mzML data 
#' xml_data <- pspecterlib::get_scan_metadata(MSPath = system.file("extdata", "Example.mzML", package = "isoforma"))
#' 
#' # Pull scan numbers
#' Sequence <- "SGRGKGGKGLGKGGAKRHRKVLRDNIQGITKPAIRRLARRGGVKRISGLIYEETRTVLKTFLENVIRDSVTYTEHARRKTVTAMDVVYALKRQGRTLYGFGG"
#' Modifications <- "Acetyl,X(1^,5,8,12,16)[2];Methyl,K(79)[1];Oxidation,M(84)[1]"
#' Modified_Sequences <- pspecterlib::multiple_modifications(Sequence, Modification, ReturnUnmodified = TRUE)
#' Scan_Numbers <- pull_scan_numbers(Sequence = Modified_Sequences[5], ScanMetadata = xml_data, RTStart = 100, RTEnd = 110)
#' 
#' # Sum Peaks
#' 
#' 
#' 
#'
#'
#' }
#'
#'
#' @export
isoforma_pipeline <- function(Sequences,
                              SummedSpectra,
                              MaxFragmentCharge,
                              ActivationMethod,
                              IonGroup,
                              CorrelationScore = 0.7,
                              PPMThreshold = 10,
                              IsotopeAlgorithm = "Rdisop",
                              MinimumAbundance = 0.01,
                              Messages = FALSE,
                              Top = 8) {

  ##################
  ## RUN ISOFORMA ##
  ##################

  # Steps 1 and 2, which are peak selection and peak summing, can be accomplished
  # outside of this function using pull_scan_numbers() and sum_ms2_spectra().

  # 3. Calculate fragments per PTM
  FragmentsPerPTM <- fragments_per_ptm(Sequences = Sequences,
                                       SummedSpectra = SummedSpectra,
                                       MaxFragmentCharge = MaxFragmentCharge,
                                       ActivationMethod = ActivationMethod,
                                       CorrelationScore = CorrelationScore,
                                       PPMThreshold = PPMThreshold,
                                       IsotopeAlgorithm = IsotopeAlgorithm,
                                       MinimumAbundance = MinimumAbundance,
                                       Messages = Messages)

  # 4. Sum isotopes
  SummedIsotopes <- sum_isotopes(FragmentsPerPTM)

  # 5. Calculate abundance matrix
  AbundanceMatrix <- abundance_matrix(SummedIsotopes = SummedIsotopes,
                                      IonGroup = IonGroup)

  # 6. Calculate Proportions
  Proportions <- calculate_proportions(
    AbundanceMatrix = AbundanceMatrix,
    Top = Top,
    IncludePlot = TRUE
  )

  # Return all objects
  return(
    list(
      FragmentsPerPTM,
      SummedIsotopes,
      AbundanceMatrix,
      Proportions[[1]],
      Proportions[[2]]
    )
  )

}
