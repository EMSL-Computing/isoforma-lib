Glossary <- data.table::fread(
  system.file("extdata", "Unimod_v20220602.csv", package = "ProteoMatch")
)

#' A function to generate a list of properly formatted proforma strings 
#' 
#' @param Sequence A valid protein sequence. It may start with "M." which gets removed. Required.
#' @param Modification An IsoForma modifications annotation written as 
#'     "PTM,Residue(Positions)[Number of Modifications]". Required. 
#'     
#' @examples
#' \dontrun{
#' 
#' # Single modification of a single PTM examples 
#' create_modifications_object("TRICITIES", "Methyl,I(3)[1]")
#' create_modifications_object("TRICITIES", "Methyl,I(3,5,7)[1]")
#' 
#' # Multiple modifications of a single PTM example
#' create_modifications_object("TRICITIES", "Methyl,I(3,5,7)[2]") 
#' 
#' # Multiple modifications with a fixed position
#' create_modifications_object("TRICITIES", "Methyl,I(3^,5,7)[2]") 
#' 
#' # Multiple modifications with two fixed positions 
#' create_modifications_object("TRICITIES", "Methyl,I(3^,5,7^)[2]")
#' 
#' # Multiple modifications of a single PTM with any "X" residue
#' create_modifications_object("TRICITIES", "Methyl,X(1,2,3,4,5,6,7,8,9)[2]")
#' 
#' # Multiple modifications with multiple PTMs examples 
#' create_modifications_object("TRICITIES", "Methyl,T(1)[1];Acetyl,X(2,4,9)[1]")
#' create_modifications_object("TRICITIES", "Methyl,T(1,2,3,4,5,6,7,8,9)[1];Acetyl,X(2,4,9)[1]")
#' create_modifications_object("TRICITIES", "Methyl,T(1^,2,3,4,5,6,7^,8,9)[3];Acetyl,X(2,4,9)[1]") 
#' 
#' }
#' @export
create_modifications_object <- function(Sequence, Modification) {
  
  ##################
  ## CHECK INPUTS ##
  ##################
  
  # Check that sequence is a character string first
  if (!is.character(Sequence)) {
    stop("Sequence must be a string.")
  }
  
  # Remove M. and any "." from the sequence 
  Sequence <- Sequence %>% 
    gsub(pattern = "M.", replacement = "", fixed = TRUE) %>%
    gsub(pattern = ".", replacement = "", fixed = TRUE)
  
  # Use pspecterlib's sequence check
  if (!pspecterlib::is_sequence(Sequence)) {
    stop("Sequence is not an acceptable peptide/protein sequence. See ?pspecterlib::is_sequence for more details.")
  }
  
  # Check that Modification is a string
  if (!is.character(Modification)) {
    stop("Modification must be a string.")
  }
  
  # The modification must have a comma, a start parenthesis, end parenthesis, 
  # start and end bracket
  checks <- c(",", "\\(", "\\)", "\\[", "\\]")
  if (all(unlist(lapply(checks, function(x) {grepl(x, Modification)}))) == FALSE) {
    stop("The proper Modification format is PTM,Residue(Positions)Number of Modifications, separated by semicolon.")
  }
  
  #####################################
  ## GENERATE LIST OF MOD DATAFRAMES ##
  #####################################
  
  # Split into dataframes to iterate through
  Mods <- lapply(unlist(strsplit(Modification, ";")), function(x) {
    
    # Get the modifications
    Mods <- dplyr::bind_rows(
      PTM  = strsplit(x, ",") %>% unlist() %>% head(1),
      Residues = x %>% strsplit(",|\\(") %>% unlist() %>% head(2) %>% tail(1),
      Positions = x %>% strsplit("\\(|\\)") %>% unlist() %>% head(2) %>% 
        tail(1) %>% strsplit(",") %>% unlist(),
      ModNum = x %>% strsplit("\\[") %>% unlist() %>% tail(1) %>% 
        gsub(pattern = "\\]", replacement = "") %>% as.numeric()
    )
    
    # If ModNum is greater than 1, generate the results with the commas in them
    if (unique(Mods$ModNum) > 1) {
      
        Mods <- tibble::tibble(
          PTM = Mods$PTM[1],
          Residues = Mods$Residues[1],
          Positions = utils::combn(Mods$Positions, Mods$ModNum[1]) %>% apply(2, function(x) paste(x, collapse = ",")),
          ModNum = Mods$ModNum[1]
        )
    }
      
    # Filter by fixed positions
    if (any(grepl("^", Mods$Positions, fixed = T))) {
      
      # Count the number of ^ in the modification call
      mod_fixed_count <- stringr::str_count(x, stringr::fixed("^"))
      
      # Subset down to appropriate number of fixed positions
      fixed_pos <- Mods$Positions[stringr::str_count(Mods$Positions, stringr::fixed("^")) == mod_fixed_count]
      
      # Now, remove the symbols and make the Mods tibble
      Mods <- tibble::tibble(
        PTM = Mods$PTM[1],
        Residues = Mods$Residues[1],
        Positions = gsub("^", "", fixed_pos, fixed = T),
        ModNum = Mods$ModNum[1]
      )
      
    } 
    
    return(Mods)
    
  })
  
  ################
  ## CHECK MODS ##
  ################
  
  # Create a dataframe that is easier to check
  checkMod <- do.call(rbind, Mods)
  
  # Check that the modification is in the glossary
  if (all(unique(checkMod$PTM) %in% Glossary$Modification) == FALSE) {
    stop(paste("Modification", checkMod$PTM, "is not in the backend glossary."))
  }
  
  # Check that the residues are acceptable options
  if (any(unique(checkMod$Residues) %in% c("B", "J", "O", "U", "Z"))) {
    stop("B, J, O, U, and Z are not acceptable inputs for residues.")
  }
  
  #####################
  ## COMBINE OPTIONS ##
  #####################
  
  # Generate all combinations of PTMs
  allCombs <- expand.grid(lapply(Mods, function(x) {1:nrow(x)}))
  
  # If there are more than 300 combinations to test, error out 
  if (nrow(allCombs) > 300) {
    stop("Isoforma currently only supports 300 modifications or less at a time.")
  }
  
  ###############################
  ## GENERATE PROFORMA STRINGS ##
  ###############################
  
  # Split the sequence
  SeqSplit <- strsplit(Sequence, "") %>% unlist()
  
  # Generate the proforma strings 
  Proforma_Strings <- lapply(1:nrow(allCombs), function(row) {
    
    # Pull entries to query
    accessList <- unlist(allCombs[row,])
    
    # Get the correct query 
    ModsToApply <- do.call(rbind, lapply(1:length(accessList), function(el) {
      Mods[[el]][accessList[el],]
    }))
    
    # Format PTM and position
    PTM_Names <- lapply(1:nrow(ModsToApply), function(theX) {rep(ModsToApply$PTM[theX], ModsToApply$ModNum[theX])}) %>% unlist()
    Position <- lapply(ModsToApply$Positions, function(val) {strsplit(val, ",") %>% unlist() %>% as.numeric()}) %>% unlist()
    
    # Use a fresh SplitSeq 
    SeqS <- SeqSplit
    
    # Expand the split sequence. If a duplicate position shows up, replace it. 
    out <- lapply(1:length(Position), function(el) {
      SeqS[Position[el]] <<- paste0(SeqSplit[Position[el]], "[", PTM_Names[el], "]")
    })
    
    return(SeqS %>% paste0(collapse = ""))
    
  }) %>% unlist()
  
  return(Proforma_Strings)
  
}
