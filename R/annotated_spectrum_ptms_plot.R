#' Overlays multiple PTM fragment identifications on one spectra
#'
#' @param SummedSpectra (peak data) A data.table with a "M/Z" and "Intensity" column, preferably
#'     calculated with the function sum_spectra. Required.
#' @param IsoformaFragments (matched peaks isoforma) A matched peaks isoforma object from fragments_per_ptm. Required. 
#' 
#' @returns (plotly graphic) An overlay of each identified peak over a spectrum
#' 
#' @examples 
#' \dontrun{
#' # Load the summed spectra
#' SumPeaks <- readRDS(system.file("extdata", "SummedPeaks_1to1to1.RDS", package = "isoforma"))
#' 
#' # Load fragment data
#' Fragments <- readRDS(system.file("extdata", "Fragments_1to1to1.RDS", package = "isoforma"))
#' 
#' annotated_spectrum_ptms_plot(SummedSpectra = SumPeaks, IsoformaFragments = Fragments)
#' }
#' @export
annotated_spectrum_ptms_plot <- function(SummedSpectra, 
                                         IsoformaFragments) {
  
  ##################
  ## CHECK INPUTS ##
  ##################
  
  # Check that summed spectra has an "M/Z" and "Intensity" in its colnames
  if (!("M/Z" %in% colnames(SummedSpectra) & "Intensity" %in% colnames(SummedSpectra))) {
    stop("SummedSpectra must have both an M/Z and Intensity column. See ?sum_spectra.")
  }
  
  # Check that PTMFragments is a "matched_peaks_isoforma" object
  if (class(IsoformaFragments) != "matched_peaks_isoforma") {
    stop("PTMFragments must be a matched_peaks_isoforma object from fragments_per_ptm.")
  }
  
  ################
  ## MAKE PLOTS ##
  ################
  
  class(SummedSpectra) <- "data.frame"
  
  # Pull out Peaks
  Peaks <- SummedSpectra %>% dplyr::select(-Abundance)
  
  # Generate zero center
  Peaks0 <- data.table::data.table("M/Z" = c(Peaks$`M/Z` - 1e-12, Peaks$`M/Z` + 1e-12),
                                   "Intensity" = 0)
  
  # Bind and order
  Peaks <- rbind(Peaks0, Peaks)
  Peaks <- Peaks[order(Peaks$`M/Z`),]
  
  # Generate blank plotly
  p <- plotly::plot_ly()
  
  # Generate base plotly with PeakData
  p <- plotly::add_trace(p, x = Peaks$`M/Z`, y = Peaks$Intensity,
                         type = "scatter", mode = "lines+markers", line = list(color = "black"),
                         marker = list(color = "black", opacity = 0), name = "Spectrum", hoverinfo = "text",
                         hovertext = paste("M/Z:", round(Peaks$`M/Z`, 3), "<br>Intensity:",
                                           round(Peaks$Intensity)))
  
  # Define a vector of 24 Colors
  ColorVec <- c("#800000", "#9A6324", "#808000", "#469990", "#000075", "#E6194B",
                "#F58231", "#3CB44B", "#42D4F4", "#4363D8", "#911EB4", "#A9A9A9",
                "#800000", "#9A6324", "#808000", "#469990", "#000075", "#E6194B",
                "#F58231", "#3CB44B", "#42D4F4", "#4363D8", "#911EB4", "#A9A9A9")
  
  # Rename first column of Peaks
  colnames(Peaks)[1] <- "M/Z Experimental"
  
  # Iterate through each fragment
  for (pos in 1:length(IsoformaFragments)) {
    
    # Pull fragments
    FragDF <- IsoformaFragments[[pos]]
    
    # Subset fragment dataframe
    FragDF <- FragDF[FragDF$Type == "c",]
    
    if (nrow(FragDF) > 0) {
      
      class(FragDF) <- c("data.frame", "data.table")
     
      # Merge and Zero Center
      FragDF <- merge(FragDF, Peaks, by = "M/Z Experimental")
      Frag0 <- data.table::data.table("M/Z Experimental" = c(FragDF$`M/Z Experimental` - 1e-12, 
                                                             FragDF$`M/Z Experimental` + 1e-12), "Intensity" = 0)
      
      # Bind and order
      FragDF <- dplyr::bind_rows(Frag0, FragDF)
      FragDF <- FragDF[order(FragDF$`M/Z Experimental`),]
      
      # Make plot
      p <- plotly::add_trace(p, x = FragDF$`M/Z Experimental`, y = FragDF$Intensity,
                             type = "scatter", mode = "lines+markers", line = list(color = ColorVec[pos]),
                             marker = list(color = ColorVec[pos], opacity = 0), name = names(IsoformaFragments)[pos], 
                             hoverinfo = "text", hovertext = paste("M/Z:", round(FragDF$`M/Z Experimental`, 3), 
                                                                   "<br>Intensity:", round(FragDF$Intensity), "<br>Ion:", paste0(FragDF$Ion, 
                                                                                                                                 "<sup>", FragDF$Z, "</sup> ", FragDF$Isotope, sep = ""))) 
       
    }
    
  }
  
  # Fix plot layout
  p <- p %>% plotly::layout(yaxis = list(title = "Intensity"),
                            legend = list(orientation = "h", borderwidth = 3))
  
  return(p %>% plotly::toWebGL())
  
}