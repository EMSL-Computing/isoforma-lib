#' Runs the 
#' 
#' @description Sum spectra based on peak selection from the pspecterlib ScanMetadata object, 
#'    or provide a list of pspecterlib::peak_data objects to sum. 
#' 
#' @param UnmodifiedSequence (character) The unmodified peptide sequence. Required.
#' @param Modifications (character) List of modifications to test. See ?pspecterlib::multiple_modifications for example. Required.
#' @param RTStart (numeric) The retention time at which to start summing applicable MS2 spectra. Required if 
#'    PeakDataList is NULL. See examples. 
#' @param RTEnd (numeric) The retention time at which to stop summing applicable MS2 spectra. Required if PeakDataList 
#'    is NULL. See examples. 
#' @param MSPath (character) The path to the mzML or raw file. Required if PeakDataList is NULL. See examples.
#' @param PeakDataList (list of peak data objects) A list of peak_data scan objects. 
#'    Required if no MSPath is supplied. See examples for more details.
#' @param OutputPath (character) The path to output the isoforma results folder. 
#' @param IsotopeAlgorithm (character) "isopat" uses the isopat package to calculate isotopes,
#'     while "Rdisop" uses the Rdisop package. Though more accurate, Rdisop has been
#'     known to crash on Windows computers when called iteratively more than 1000 times.
#'     Default is Rdisop, though isopat is an alternative.
#' @param IonGroup (character) Must be an ion group to subset data down to (a-c or x-z). Default is "c".
#' @param MassWindow (numeric) m/z values to sum isotopes over. Written Daltons and is applied as -/+ 10 Da. Default is 10.
#' @param PPMRound (numeric)
#' 
#' @returns Many objects
#' 
#' @examples
#' \dontrun{
#' ####################################
#' ## EXAMPLE 1: SKIP PEAK SELECTION ##
#' ####################################
#' 
#' # Load 3 peaks to sum 
#' PeakDataList <- list(
#'  readRDS(system.file("extdata", "PeakData_1to1to1_1.RDS", package = "isoforma")),
#'  readRDS(system.file("extdata", "PeakData_1to1to1_2.RDS", package = "isoforma")),
#'  readRDS(system.file("extdata", "PeakData_1to1to1_3.RDS", package = "isoforma"))
#' )
#' 
#' # Run main isoforma function
#' isoforma_pipeline(
#'    UnmodifiedSequence = "LQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG",
#'    Modification =  "6.018427,V(17,26,70)[1]",
#'    PeakDataList = PeakDataList,
#'    OutputPath = "~/Downloads/Isoforma_Brunner_Example", 
#'    IsotopeAlgorithm = "Rdisop", 
#'    PrecursorCharge = 11, 
#'    ActivationMethod = "ETD",
#'    IonGroup = "c",
#' 
#' 
#' 
#' 
#' )
#' 
#' 
#' 
#' #########################################
#' ## EXAMPLE 2: AUTOMATED PEAK SELECTION ##
#' #########################################
#' 
#' }
#' 
#'
#' @export
isoforma_pipeline <- function(unmodifiedsequence,
                              modifications,
                              rt_start,
                              rt_end,
                              mzML_path,
                              output_path,
                              IsotopeAlgorithm,
                              IonGroup = "c",
                              mass_window = 10,
                              mz_round = 3,
                              correlation_score = 0.7) {

  #################
  ## LOAD INPUTS ##
  #################

  # Check that unmodified_sequence is a valid input
  if (!pspecterlib::is_sequence(unmodified_sequence)) {
    stop(paste(unmodified_sequence, "is not an acceptable unmodified amino acid sequence."))
  }

  # Modifications should be a string
  if (!is.character(modifications)) {
    stop("modifications must be a string")
  }

  # Load mzML file
  if (grepl(".mzML", mzML_path, fixed = TRUE) == FALSE) {
    stop("mzML_path must be a mzML file.")
  }
  ScanMetadata <- pspecterlib::get_scan_metadata(mzML_path)

  ##################
  ## RUN ISOFORMA ##
  ##################

  # 0. Create output folder
  #out <- file.path(output_path, paste(sub(".mzml", "", basename(mzML_path), ignore.case = TRUE), "isoforma_out", sep = "_"))
  out <- output_path

  if (!dir.exists(out)) {dir.create(out)}

  # 1. Output scan number selection
  ScanNums <- pull_scan_numbers(Sequence = unmodified_sequence,
                                Modification = modifications,
                                ScanMetadata = ScanMetadata,
                                RTStart = rt_start,
                                RTEnd = rt_end,
                                AsDataframe = TRUE,
                                MassWindow = mass_window)

  ScanNumbers <- ScanNums[[1]]
  ExactMass <- ScanNums[[2]]
  if (is.null(ScanNumbers) | any(TRUE %in% ScanNumbers$Use) == FALSE) {
    utils::write.table("No MS2 Scan Numbers detected.", file.path(out, "Message.txt"))
    return(NULL)
  }
  utils::write.csv(ScanNumbers, file = file.path(out, "ScanNumbers.csv"), row.names = F, quote = F)

  # 2. Sum Spectra
  SumSpectra <- sum_spectra(ScanMetadata = ScanMetadata,
                            ScanNumbers = ScanNumbers$`Scan Number`[ScanNumbers$Use],
                            MZRound = mz_round,
                            Percent = 100)
  utils::write.csv(SumSpectra, file = file.path(out, "SumSpectra.csv"), row.names = F, quote = F)
  write_mgf_simple(round(SumSpectra$`M/Z`, 4), round(SumSpectra$Intensity, 4), file.path(out, "SumSpectra.mgf"),
                   title = strsplit(out, "/", fixed = T) %>% unlist() %>% tail(1),
                   pepmass = ExactMass, charge = paste0(max(ScanNumbers$`Precursor Charge`), "+"),
                   rt = mean(ScanNumbers$`Retention Time`), scans = paste0(ScanNumbers$`Scan Number`, collapse = " & "))

  # 3. Calculate fragments per PTM
  FragmentsPerPTM <- fragments_per_ptm(Sequences = c(unmodified_sequence, create_modifications_object(unmodified_sequence, modifications)),
                                       SummedSpectra = SumSpectra,
                                       PrecursorCharge = max(ScanNumbers$`Precursor Charge`),
                                       ActivationMethod = unlist(ScanNumbers[ScanNumbers$Use, "Activation Method"])[1],
                                       CorrelationScore = correlation_score,
                                       IsotopeAlgorithm = IsotopeAlgorithm,
                                       PPMThreshold = 5,
                                       Messages = TRUE)

  plot <- ptm_heatmap(FragmentsPerPTM)
  ggplot2::ggsave(file.path(out, "annotated_ptms_heatmap.png"), plot)
  plot2 <- annotated_spectrum_ptms_plot(SumSpectra, FragmentsPerPTM)
  htmlwidgets::saveWidget(plot2, file.path(out, "Spectra.html"))

  # 4. Sum isotopes
  SummedIsotopes <- sum_isotopes(FragmentsPerPTM)

  # 5. Calculate abundance matrix
  AbundanceMatrix <- abundance_matrix(IonGroup, SummedIsotopes)
  utils::write.csv(AbundanceMatrix$AbundanceMatrix, file.path(out, "AbundanceMatrix.csv"), quote = F, row.names = F)

  # 6. Calculate Proportions
  Proportions <- calculate_proportions(AbundanceMatrix, modifications, IncludePlot = TRUE)
  utils::write.csv(Proportions[[1]], file.path(out, "Proportions.csv"), quote = F, row.names = F)
  ggplot2::ggsave(file.path(out, "Proportions.png"), Proportions[[2]])

  # Clear out scan metadata
  rm(ScanMetadata)

}
