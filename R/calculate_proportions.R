#' Generate the proportions of fragments from an abundance matrix
#' 
#' @description Returns a list where the first object is the calculated proportions,
#'     and the second is barplot with the top n calculated proportions and their 
#'     respective error bars. 
#' 
#' @param AbundanceMatrix (abundance matrix isoforma) An abundance matrix from abundance_matrix(). Required. 
#' @param IncludePlot (logical) TRUE/FALSE to indicate whether a plot should also be returned. Default is TRUE.
#' @param Top (integer) The top N proportions to return, ranked from higest to lowest. Default is 8.
#' 
#' @return (list) The first item is a data.table and the second is a plot
#' 
#' @examples
#' \dontrun{
#' # Load summed isotopes 
#' SumIso <- readRDS(system.file("extdata", "SummedIsotopes_1to1to1.RDS", package = "isoforma"))
#'
#' # Generate the matrix
#' AbundanceMat <- abundance_matrix(SummedIsotopes = SumIso, IonGroup = "c")
#' 
#' # Calculate Proportions
#' calculate_proportions(AbundanceMat)
#' }
#' 
#' @export
calculate_proportions <- function(AbundanceMatrix, 
                                  IncludePlot = TRUE, 
                                  Top = 8) {
  
  ##################
  ## CHECK INPUTS ##
  ##################
  
  # Check the ion group designation
  if (inherits(AbundanceMatrix, "abundance_matrix_isoforma") == FALSE) {
    stop("AbundanceMatrix must be calculated from the abundance_matrix function")
  }
  
  # Check that include plot is a logical
  if (!is.logical(IncludePlot)) {
    stop("IncludePlot must be TRUE or FALSE.")
  }
  
  # Top should be a number
  if (!is.numeric(Top)) {
    stop("Top should be a single numeric")
  }
  Top <- round(abs(Top))
  
  ###########################
  ## CALCULATE PROPORTIONS ##
  ###########################

  class(AbundanceMatrix) <- c("data.frame", "data.table")
  AbundanceMatrix <- data.table::data.table(AbundanceMatrix)
  
  # 1. Order columns by increasing position
  SummedPosition <- data.table::data.table(
    "Name" = colnames(AbundanceMatrix)[colnames(AbundanceMatrix) != "Ion"]
  ) %>%
    dplyr::mutate(
      Positions = lapply(Name, function(theName) {
        strsplit(theName, "@| & ") %>% 
          unlist() %>% 
          .[c(FALSE, TRUE)] %>% 
          gsub(pattern = "[[:alpha:]]", replacement = "") %>% 
          as.numeric()
      }), 
      MinPosition = lapply(Positions, function(x) {min(x)}) %>% unlist(),
      MaxPosition = lapply(Positions, function(x) {max(x)}) %>% unlist(),
      PositionSum = lapply(Positions, function(x) {sum(x)}) %>% unlist()
    ) %>%
    dplyr::arrange(MinPosition, MaxPosition)
  theOrder <- c("Ion", unlist(SummedPosition$Name))
  AbundanceMatrix <- AbundanceMatrix %>% dplyr::select(theOrder)
  
  # 2. Determine comparison ranges 
  ComparisonRanges <- data.table::data.table(
    Name = SummedPosition$Name, 
    Start = SummedPosition$MinPosition,
    Stop = c(SummedPosition$MaxPosition[2:nrow(SummedPosition)]-1, 
             max(AbundanceMatrix$Ion %>% gsub(pattern = "[[:alpha:]]", replacement = "") %>% as.numeric()))
  )
  
  # 3. Divide each intensity by the last intensity in the row
  AbunMat <- AbundanceMatrix[,2:ncol(AbundanceMatrix)] 
  End <- ncol(AbunMat)
  AbunMat <- as.matrix(AbunMat) / (AbunMat + AbunMat[,End])
  AbunMat <- data.table::data.table(AbunMat)
  AbunMat$Ion <- AbundanceMatrix$Ion %>% gsub(pattern = "[[:alpha:]]", replacement = "")
  AbunMat <- AbunMat %>% dplyr::relocate(Ion)
  
  # 4. Calculate Medians
  Medians <- do.call(dplyr::bind_rows, lapply(ComparisonRanges$Name[1:nrow(ComparisonRanges)-1], function(Name) {
    
    # Define range
    Range <- unlist(ComparisonRanges[ComparisonRanges$Name == Name, "Start"]):unlist(ComparisonRanges[ComparisonRanges$Name == Name, "Stop"])
    Values <- AbunMat[AbunMat$Ion %in% Range, Name] %>% .[!is.na(.) & . != 0.5] 
    
    return(
      c(Name = Name, 
        Pre = Values %>% median() %>% round(8),
        Error = Values %>% sd() %>% round(8)
      )
    )
    
  }))
  Medians <- dplyr::bind_rows(Medians, c(Name = ComparisonRanges$Name[nrow(ComparisonRanges)], Pre = NA, Error = NA))
  
  # 5. Ensure each value is higher than the next
  Pre <- Medians
  Pre$Pre <- as.numeric(Pre$Pre)
  PreNums <- which(!is.na(Pre$Pre))
  
  if (length(PreNums) > 1) {
    PreTest <- Pre[PreNums,] %>%
      dplyr::mutate(
        Larger = Pre > dplyr::lag(Pre)
      ) 
    PreTest$Larger[1] <- TRUE
    PreTest <- PreTest[PreTest$Larger,]
  } else if (length(PreNums) == 1) {
    PreTest <- Pre[PreNums,]
  } else {
    message("No proportions found")
    return(NULL)
  }
  
  # 6. Finally, calculate the real proportion
  if (nrow(PreTest) > 1) {
    PreTest <- PreTest %>% 
      dplyr::mutate(
        Proportion = Pre - dplyr::lag(Pre)
      )
    PreTest$Proportion[1] <- PreTest$Pre[1]
    PreTest <- PreTest %>% dplyr::mutate(
      Proportion = ifelse(Proportion <= 0, NA, Proportion)
    )
  } else if (nrow(PreTest) == 1) {
    PreTest$Proportion <- PreTest$Pre
  }
  
  Proportions <- dplyr::bind_rows(
    data.table::data.table(
      "Name" = unlist(Pre$Name)[unlist(Pre$Name) %in% PreTest$Name == FALSE],
      "Proportion" = NA
    ),
    PreTest[,c("Name", "Proportion")]
  )

  
  if (sum(Proportions$Proportion, na.rm = T) > 1) {
    Proportions <- Proportions %>%
      dplyr::mutate(
        Proportion = ifelse(is.na(Proportion), 0, Proportion),
        CumSum = cumsum(Proportion),
        Proportion = ifelse(CumSum > 1 | Proportion < 0.0001, NA, Proportion)
      ) %>% 
      dplyr::select(-CumSum)
  }
  
  # 8. If there is only one NA left in the dataset, we can calculate it 
  if (sum(is.na(Proportions$Proportion)) == 1) {
    Proportions$Proportion[is.na(Proportions$Proportion)] <- 1 - sum(Proportions$Proportion, na.rm = T)
  }
  
  # 9. Order and add a standard error
  Proportions$Name <- factor(Proportions$Name, levels = ComparisonRanges$Name)
  Proportions <- Proportions %>% dplyr::arrange(Name)
  Proportions$Error <- Medians$Error
  Proportions <- Proportions %>% dplyr::rename(Modification = Name)

  # Make percentage plot if necessary
  if (IncludePlot) {
    
    # Modify data frame for better plotting
    PercentPlot <- Proportions
    PercentPlot$x <- "x"
    PercentPlot$Modification <- factor(PercentPlot$Modification, levels = PercentPlot$Modification)
    
    if (nrow(PercentPlot) > Top) {
      PercentPlot <- PercentPlot %>%
        dplyr::group_by(Proportion) %>%
        dplyr::arrange(-Proportion) 
      PercentPlot <- PercentPlot[1:Top,]
      PercentPlot <- rbind(PercentPlot, list("Modification" = "Other", "Percentages" = 1 - sum(PercentPlot$Proportion), 
         "x" = "x"))  
    }
    
    PercentPlot$Error <- as.numeric(PercentPlot$Error)
    
    Plot <- ggplot2::ggplot(PercentPlot, ggplot2::aes(x = x, y = Proportion, fill = Modification)) +
      ggplot2::geom_bar(stat = "identity", position = "dodge") + ggplot2::theme_bw() +
      ggplot2::geom_errorbar(ggplot2::aes(ymin = Proportion - Error, ymax = Proportion + Error), width = 0.2, position = ggplot2::position_dodge(.9)) +
      ggplot2::theme(axis.title.x = ggplot2::element_blank(), axis.text.x = ggplot2::element_blank()) 
    
    return(list(Proportions, Plot))
    
  } else {return(Proportions)}
  
}