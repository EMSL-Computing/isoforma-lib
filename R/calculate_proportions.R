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
  
  # 1. Iterate through each column and calculate the proportions, removing all 0.5 and na
  Prop_Early <- lapply(2:(ncol(AbundanceMatrix) - 1), function(col) {
    Props <- AbundanceMatrix[,col] / (AbundanceMatrix[,col] + AbundanceMatrix[,(col+1)])
    return(Props[!is.na(Props) & Props != 0.5])
  })
  
  # 2. Calculate means and confidence interval using a negative binomial distribution
  Means <- do.call(dplyr::bind_rows, lapply(Prop_Early, function(x) {
    
    Values <- x %>% unlist()
    
    nloglikbeta <- function(mu, sig) {
      alpha = mu^2*(1-mu)/sig^2-mu
      beta = alpha*(1/mu-1)
      return(suppressWarnings(-1 * sum(dbeta(x, alpha, beta, log=TRUE))))
    }
    mle_fit <- suppressWarnings(stats4::mle(nloglikbeta, start = list(mu = mean(Values), sig = sd(Values))))
    theConfint <- stats4::confint(mle_fit)
    
    return(
      c(Pre = mle_fit@coef[1][[1]] %>% round(8),
        LowerCI = theConfint[1] %>% round(8),
        UpperCI = theConfint[3] %>% round(8)
      )
    )
    
  })) %>%
    dplyr::mutate(Name = colnames(AbundanceMatrix)[2:(ncol(AbundanceMatrix) - 1)]) %>%
    dplyr::relocate(Name)
  
  
  Means <- dplyr::bind_rows(Means, data.frame(
    Name = colnames(AbundanceMatrix)[ncol(AbundanceMatrix)], 
    Pre = 1, 
    LowerCI = 1 - as.numeric(Means$UpperCI[nrow(Means)]) + as.numeric(Means$Pre[nrow(Means)]),
    UpperCI = 1 - as.numeric(Means$LowerCI[nrow(Means)]) + as.numeric(Means$Pre[nrow(Means)])
  ))

  # 3. Ensure each value is higher than the next
  Pre <- Means
  Pre$Pre <- as.numeric(Pre$Pre)
  Pre$LowerCI <- as.numeric(Pre$LowerCI)
  Pre$UpperCI <- as.numeric(Pre$UpperCI)
  PreNums <- which(!is.na(Pre$Pre))
  
  if (length(PreNums) > 1) {
    PreTest <- Pre[PreNums,] %>%
      dplyr::mutate(
        Larger = Pre > dplyr::lag(Pre) 
      ) 
    PreTest$Larger[1] <- TRUE
    PreTest$Larger[nrow(PreTest)] <- TRUE
    PreTest <- PreTest[PreTest$Larger,]
  } else if (length(PreNums) == 1) {
    PreTest <- Pre[PreNums,]
  } else {
    message("No proportions found")
    return(NULL)
  }
  
  # 4. Finally, calculate the real proportion
  if (nrow(PreTest) > 1) {
    PreTest <- PreTest %>% 
      dplyr::mutate(
        Proportion = Pre - dplyr::lag(Pre),
        TrueLowerCI = LowerCI - dplyr::lag(Pre),
        TrueUpperCI = UpperCI - dplyr::lag(Pre)
      )
    PreTest$Proportion[1] <- PreTest$Pre[1]
    PreTest$TrueLowerCI[1] <- PreTest$LowerCI[1]
    PreTest$TrueUpperCI[1] <- PreTest$UpperCI[1]
    PreTest <- PreTest %>% dplyr::mutate(
      Proportion = ifelse(Proportion <= 0, NA, Proportion),
      TrueLowerCI = ifelse(TrueLowerCI < 0, 0, TrueLowerCI),
      TrueUpperCI = ifelse(TrueUpperCI > 1, 1, TrueUpperCI)
    )
  } else if (nrow(PreTest) == 1) {
    PreTest$Proportion <- PreTest$Pre
    PreTest$TrueLowerCI <- PreTest$LowerCI
    PreTest$TrueUpperCI <- PreTest$UpperCI
  }
  
  # Make proportions data.frame
  Proportions <- PreTest %>%
    dplyr::select(Name, Proportion, TrueLowerCI, TrueUpperCI) %>%
    dplyr::rename(LowerCI = TrueLowerCI, UpperCI = TrueUpperCI)
  
  # 5. Order and rename
  Proportions$Name <- factor(Proportions$Name, levels = colnames(AbundanceMatrix)[2:ncol(AbundanceMatrix)])
  Proportions <- Proportions %>% dplyr::arrange(Name)
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
    
    Plot <- ggplot2::ggplot(PercentPlot, ggplot2::aes(x = x, y = Proportion, fill = Modification)) +
      ggplot2::geom_bar(stat = "identity", position = "dodge") + ggplot2::theme_bw() +
      ggplot2::geom_errorbar(ggplot2::aes(ymin = LowerCI, ymax = UpperCI), width = 0.2, position = ggplot2::position_dodge(.9)) +
      ggplot2::theme(axis.title.x = ggplot2::element_blank(), axis.text.x = ggplot2::element_blank()) 
    
    return(list(Proportions, Plot))
    
  } else {return(Proportions)}
  
}