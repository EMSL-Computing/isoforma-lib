#' Select ions based on dissociation method
#' 
#' @param ActivationMethod A dissociation method of "HCD", "CID", or "ETD". Otherwise,
#'     all ions will be returned. Required.
#'
#' @export
.determine_ions <- function(activation_method) {
  
  Ions <- c("a", "b", "c", "x", "y", "z")
  
  # Select specific ions if the activation method is HCD, CID, ETD
  if (activation_method == "HCD") {Ions <- c("b", "y")} else 
  if (activation_method == "CID") {Ions <- c("a", "b", "y")} else
  if (activation_method == "ETD") {Ions <- c("c", "z")} 
  
  return(Ions)
}