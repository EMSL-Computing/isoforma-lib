#' Generate a heatmap of monoisotopic fragments for each PTM 
#' 
#' @param IsoformaFragments A matched peaks isoforma object from fragments_per_ptm. Required. 
#'
#' @examples
#' \dontrun{
#' 
#' 
#' 
#' }
#' @export    
ptm_heatmap <- function(IsoformaFragments) {
  
  Heatmap <- do.call(rbind, lapply(1:length(IsoformaFragments), function(el) {
    
    # Get object name and matched_peaks objects
    Name <- names(IsoformaFragments)[el]
    PTM_Frag <- IsoformaFragments[[el]]
    
    # Make a datatable with blanks to fill out
    Heatmap <- data.table::data.table(
      "Group" = expand.grid(PTM_Frag$`General Type` %>% unique(), 
                            1:nchar(attr(PTM_Frag, "pspecter")$Sequence)) %>% data.table::as.data.table() %>% 
        dplyr::mutate("Group" = paste0(Var1, "_", Var2)) %>% dplyr::select(Group) %>% unlist()
    )
    
    class(PTM_Frag) <- c("data.table", "data.frame")
    
    # Pull out the minimum ppm error per ion type and position, merge with blanks,
    # split into ion type groups
    Heatmap <- PTM_Frag %>%
      dplyr::mutate(
        Group = paste0(`General Type`, "_", Position)
      ) %>%
      dplyr::select(Group, `PPM Error`) %>%
      dplyr::group_by(Group) %>%
      tidyr::nest() %>% 
      dplyr::mutate(
        # Mowei approved this 
        "PPM Error" = purrr::map(data, function(x) {max(x$`PPM Error`)}) %>% unlist()
      ) %>%
      dplyr::select(-data) %>%
      merge(Heatmap, "Group", all = TRUE) %>% 
      dplyr::mutate(
        "Sample" = paste0(Name, "_", lapply(Group, function(x) {strsplit(x, "_") %>% unlist() %>% head(1)}) %>% unlist()),
        "Position" = lapply(Group, function(x) {strsplit(x, "_") %>% unlist() %>% tail(1)}) %>% unlist() %>% as.numeric()
      ) %>%
      dplyr::select(-Group)
    
    return(Heatmap)
    
  }))
  
  # Determine order 
  theOrder <- data.table::data.table(
    "Sample" = unique(Heatmap$Sample),
    "Position" = unique(Heatmap$Sample) %>% gsub(pattern = "\\D", replacement = "") %>% as.numeric(), 
    "Ion" = unique(Heatmap$Sample) %>% lapply(function(x) {strsplit(x, "_") %>% unlist() %>% tail(1)}) %>% unlist()
  ) %>% 
    dplyr::arrange(Position) %>%
    dplyr::arrange(Ion) %>% 
    dplyr::select(Sample) %>%
    unlist()
  
  # Order data
  Heatmap$Sample <- factor(Heatmap$Sample, levels = theOrder)
  
  # Make plot
  plot <- ggplot2::ggplot(Heatmap, ggplot2::aes(x = Position, y = Sample, fill = `PPM Error`)) +
    ggplot2::geom_tile() + ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1)) +
    ggplot2::scale_fill_gradientn(colours = c("blue", "grey", "red"),
                                  values = c(0,0.5,1), na.value = "black") +
    ggplot2::xlab("Position") + ggplot2::ylab("PTM")  
  
  return(plot)
  
}