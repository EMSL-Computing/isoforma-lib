#' Generate the proportions of fragments from an abundance matrix
#' 
#' @param AbundanceMatrix An abundance matrix from abundance_matrix(). Required. 
#' @param IncludePlot TRUE/FALSE to indicate whether a plot should also be returned. Default is TRUE.
#' 
#' @export
calculate_proportions <- function(AbundanceMatrix, IncludePlot = TRUE, top = 8) {
  
  ##################
  ## CHECK INPUTS ##
  ##################
  
  # Check the ion group designation
  if (class(AbundanceMatrix) != "abundance_matrix_isoforma") {
    stop("AbundanceMatrix must be calculated from the abundance_matrix function")
  }
  
  # Check that include plot is a logical
  if (!is.logical(IncludePlot)) {
    stop("IncludePlot must be TRUE or FALSE.")
  }
  
  ###########################
  ## CALCULATE PROPORTIONS ##
  ###########################
  
  # Extract the matrix
  Abun <- AbundanceMatrix$AbundanceMatrix
  
  # 1. Order columns by increasing position
  SummedPosition <- data.table::data.table(
    "Name" = colnames(Abun)[colnames(Abun) != "Ions"]
  ) %>%
    dplyr::mutate(
      Positions = lapply(Name, function(theName) {theName %>% strsplit(" & ") %>% unlist() %>% 
        lapply(function(x) {gsub("[^[:digit:]]", "", x)}) %>% 
        unlist() %>%
        as.numeric()}), 
      MinPosition = lapply(Positions, function(x) {min(x)}) %>% unlist(),
      MaxPosition = lapply(Positions, function(x) {max(x)}) %>% unlist(),
      PositionSum = lapply(Positions, function(x) {sum(x)}) %>% unlist()
    ) %>%
    dplyr::arrange(MinPosition, MaxPosition)
  theOrder <- c("Ions", unlist(SummedPosition$Name))
  Abun <- Abun[,theOrder]
  
  # 2. Divide each intensity by the last intensity in the row
  AbunMat <- Abun[,2:ncol(Abun)] 
  End <- ncol(AbunMat)
  AbunMat <- as.matrix(AbunMat) / (AbunMat + AbunMat[,End])
  AbunMat <- data.table::data.table(AbunMat)
  
  # 3. Filter out options that contain only NA, 0.5, or NA and 0.5 
  toFilter <- lapply(apply(AbunMat, 1, unique), function(x) {
    x[is.na(x)] <- 0
    zeroTest <- 0 %in% x & length(x) == 1
    zeroFiveTest <- 0.5 %in% x & length(x) == 1
    bothTest <- 0 %in% x & 0.5 %in% x & length(x) == 2
    return(any(zeroTest, zeroFiveTest) | bothTest)
  }) %>% unlist()
  toKeep <- which(toFilter == FALSE)
  AbunMat <- AbunMat[toKeep,]
  
  # 4. Keep only leftmost not 0.5
  leftmost_not0.5 <- lapply(1:nrow(AbunMat), function(x) {
    which(AbunMat[x,] == 0.5) %>% head(1) - 1
  }) %>% unlist()
  
  for (row in 1:nrow(AbunMat)) {
    leftmost <- leftmost_not0.5[row]
    if (is.na(leftmost)) {return(NULL)} 
    else if (leftmost == 1) {
      AbunMat[row,] <- AbunMat[row,]
    } else if (leftmost == 2) {
      AbunMat[row, 1] <- NA
    } else {
      AbunMat[row, 1:(leftmost-1)] <- NA
    }
  }
  
  # 5. Convert 0.5 to NA and get row medians
  AbunMat[AbunMat == 0.5] <- NA
  median2 <- function(x) {
    val <- median(x, na.rm = T)
    if (is.infinite(val)) {return(NA)} else (return(val))
  }
  Medians <- apply(AbunMat, 2, median2)
  
  # 6. Ensure each value is higher than the next
  Pre <- data.table::data.table(
    Modification = names(Medians),
    Pre = Medians
  ) 
  
  PreNums <-  which(!is.na(Pre$Pre))
  
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
  
  # 7. Finally, calculate the real proportion
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
      "Modification" = unlist(Pre$Modification)[unlist(Pre$Modification) %in% PreTest$Modification == FALSE],
      "Proportion" = NA
    ),
    PreTest[,c("Modification", "Proportion")]
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

  # Make percentage plot if necessary
  if (IncludePlot) {
    
    # Modify data frame for better plotting
    PercentPlot <- Proportions
    PercentPlot$x <- "x"
    PercentPlot$Modification <- factor(PercentPlot$Modification, levels = PercentPlot$Modification)
    
    if (nrow(PercentPlot) > top) {
      PercentPlot <- PercentPlot %>%
        dplyr::group_by(Proportion) %>%
        dplyr::arrange(-Proportion) 
      PercentPlot <- PercentPlot[1:top,]
      PercentPlot <- rbind(PercentPlot, list("Modification" = "Other", "Percentages" = 1 - sum(PercentPlot$Proportion), 
         "x" = "x"))  
    }
    
    Plot <- ggplot2::ggplot(PercentPlot, ggplot2::aes(x = x, y = Proportion, fill = Modification)) +
      ggplot2::geom_bar(stat = "identity", position = "stack") + ggplot2::theme_bw() +
      ggplot2::theme(axis.title.x = ggplot2::element_blank(), axis.text.x = ggplot2::element_blank()) +
      ggplot2::ylim(c(0,1))
    
    return(list(Proportions, Plot))
    
  } else {return(Proportions)}
  
}