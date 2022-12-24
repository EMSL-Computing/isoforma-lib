#' Run the main isoforma function
#'
#' @param parameters_path The path to the parameters file. Required.
#' @param mzML_path The path to the mzML file. Required.
#' @param output_path The path to output the isoforma results folder. Default is same as mzML_path.
#' 
#' @export
isoforma_pipeline <- function(parameters_path,
                              mzML_path,
                              output_path) {
    
    #################
    ## LOAD INPUTS ##
    #################
    
    # Load parameters file
    if (!file.exists(parameters_path) | grepl(parameters_path, ".xlsx", fixed = T)) {
      stop("parameters_path must exist and be an xlsx file.")
    }
    parameters <- read.xlsx(parameters_path, 1)
    if (all(c("Parameter", "Setting", "Description") %in% colnames(parameters)) == FALSE) {
      stop("parameters_path file must have the following columns: Parameter, Setting, Description")
    }
    if (all(c("Mass Window", "MZRound", "Correlation Score") %in% parameters$Parameter) == FALSE) {
      stop("parameters_path file must have the following Parameters: Mass Window, MZRound, Correlation Score")
    }
    
    # Load mzML file
    if (grepl(".mzML", mzML_path, fixed = TRUE) == FALSE) {
      stop("mzML_path must be a mzML file.")
    }
    ScanMetadata <- get_scan_metadata(mzML_path)
    
    ##################
    ## RUN ISOFORMA ##
    ##################
    
    # 0. Create output folder
    #out <- file.path(output_path, paste(sub(".mzml", "", basename(mzML_path), ignore.case = TRUE), "isoforma_out", sep = "_"))
    out <- output_path
    
    if (!dir.exists(out)) {dir.create(out)}
    
    # 1. Output scan number selection
    ScanNums <- pull_scan_numbers(Sequence = targets$UnmodifiedSequence,
                                  Modification = targets$Modifications,
                                  ScanMetadata = ScanMetadata,
                                  RTStart = targets$RTStart,
                                  RTEnd = targets$RTEnd,
                                  AsDataframe = TRUE,
                                  MassWindow = parameters[parameters$Parameter == "Mass Window", "Setting"] %>% as.numeric())
    
    ScanNumbers <- ScanNums[[1]]
    ExactMass <- ScanNums[[2]]
    if (is.null(ScanNumbers) | any(TRUE %in% ScanNumbers$Use) == FALSE) {
      write.table("No MS2 Scan Numbers detected.", file.path(out, "Message.txt"))
      return(NULL)
    }
    write.csv(ScanNumbers, file = file.path(out, "ScanNumbers.csv"), row.names = F, quote = F)
    
    # 2. Sum Spectra
    SumSpectra <- sum_spectra(ScanMetadata = ScanMetadata,
                              ScanNumbers = ScanNumbers$`Scan Number`[ScanNumbers$Use],
                              MZRound = parameters[parameters$Parameter == "MZRound", "Setting"] %>% as.numeric(),
                              Percent = 100)
    write.csv(SumSpectra, file = file.path(out, "SumSpectra.csv"), row.names = F, quote = F)
    write_mgf_simple(round(SumSpectra$`M/Z`, 4), round(SumSpectra$Intensity, 4), file.path(out, "SumSpectra.mgf"),
                     title = strsplit(out, "/", fixed = T) %>% unlist() %>% tail(1),
                     pepmass = ExactMass, charge = paste0(max(ScanNumbers$`Precursor Charge`), "+"), 
                     rt = mean(ScanNumbers$`Retention Time`), scans = paste0(ScanNumbers$`Scan Number`, collapse = " & "))
    
    # 3. Calculate fragments per PTM
    FragmentsPerPTM <- fragments_per_ptm(Sequences = c(targets$UnmodifiedSequence, create_modifications_object(targets$UnmodifiedSequence, targets$Modifications)),
                                         SummedSpectra = SumSpectra,
                                         PrecursorCharge = max(ScanNumbers$`Precursor Charge`),
                                         ActivationMethod = unlist(ScanNumbers[ScanNumbers$Use, "Activation Method"])[1],
                                         CorrelationScore = parameters[parameters$Parameter == "Correlation Score", "Setting"] %>% as.numeric(),
                                         PPMThreshold = 5,
                                         Messages = TRUE)
    
    plot <- ptm_heatmap(FragmentsPerPTM)
    ggplot2::ggsave(file.path(out, "annotated_ptms_heatmap.png"), plot)
    plot2 <- annotated_spectrum_ptms_plot(SumSpectra, FragmentsPerPTM)
    htmlwidgets::saveWidget(plot2, file.path(out, "Spectra.html"))
    
    
    # 4. Sum isotopes
    SummedIsotopes <- sum_isotopes(FragmentsPerPTM)
    
    # 5. Calculate abundance matrix
    AbundanceMatrix <- abundance_matrix("c", SummedIsotopes)
    write.csv(AbundanceMatrix$AbundanceMatrix, file.path(out, "AbundanceMatrix.csv"), quote = F, row.names = F)
    
    # 6. Calculate Proportions
    Proportions <- calculate_proportions(AbundanceMatrix, targets$Modifications, IncludePlot = TRUE)
    write.csv(Proportions[[1]], file.path(out, "Proportions.csv"), quote = F, row.names = F)
    ggplot2::ggsave(file.path(out, "Proportions.png"), Proportions[[2]])
    
    # Clear out scan metadata 
    rm(ScanMetadata)
    
  }