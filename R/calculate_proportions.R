#' Generate the proportions of fragments from an abundance matrix
#' 
#' @param AbundanceMatrix An abundance matrix from abundance_matrix(). Required. 
#' @param IncludePlot TRUE/FALSE to indicate whether a plot should also be returned. Default is TRUE.
#' 
#' @export
calculate_proportions <- function(AbundanceMatrix, 
                                  IncludePlot = TRUE, 
                                  top = 8) {
  
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

  # Make percentage plot if necessary
  if (IncludePlot) {
    
    # Modify data frame for better plotting
    PercentPlot <- Proportions
    PercentPlot$x <- "x"
    PercentPlot$Name <- factor(PercentPlot$Name, levels = PercentPlot$Name)
    
    if (nrow(PercentPlot) > top) {
      PercentPlot <- PercentPlot %>%
        dplyr::group_by(Proportion) %>%
        dplyr::arrange(-Proportion) 
      PercentPlot <- PercentPlot[1:top,]
      PercentPlot <- rbind(PercentPlot, list("Name" = "Other", "Percentages" = 1 - sum(PercentPlot$Proportion), 
         "x" = "x"))  
    }
    
    PercentPlot$Error <- as.numeric(PercentPlot$Error)
    
    Plot <- ggplot2::ggplot(PercentPlot, ggplot2::aes(x = x, y = Proportion, fill = Name)) +
      ggplot2::geom_bar(stat = "identity", position = "dodge") + ggplot2::theme_bw() +
      ggplot2::geom_errorbar(ggplot2::aes(ymin = Proportion - Error, ymax = Proportion + Error), width = 0.2, position = ggplot2::position_dodge(.9)) +
      ggplot2::theme(axis.title.x = ggplot2::element_blank(), axis.text.x = ggplot2::element_blank()) 
    
    return(list(Proportions, Plot))
    
  } else {return(Proportions)}
  
}