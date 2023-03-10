% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pull_scan_numbers.R
\name{pull_scan_numbers}
\alias{pull_scan_numbers}
\title{Return a vector of scan numbers based on a Precursor Mass, mass tolerance, and
isolation window filtering}
\usage{
pull_scan_numbers(
  Sequence,
  Modification,
  ScanMetadata,
  RTStart,
  RTEnd,
  MassWindow = 10,
  MinMS1Matches = 3,
  ActivationMethod = "ETD",
  AsDataframe = FALSE,
  MinAbundance = 5
)
}
\arguments{
\item{Sequence}{A valid amino acid sequence. Any non-traditional symbols are
removed and "M." is ignored. Required.}

\item{Modification}{An IsoFormA modifications annotation written as
"PTM,Residue(Positions)\link{Number of Modifications}". Required.}

\item{ScanMetadata}{An object of the scan_metadata class from pspecterlib
get_scan_metadata. Required.}

\item{MassWindow}{A m/z window for acceptable precursor masses. Default is 5.}

\item{MinMS1Matches}{Minimum number of MS1 peaks from the isotopic distribution
of the full peptide sequence that fall within the MS2 window. Default is 3.}

\item{ActivationMethod}{The method used by the MS instrument to dissociate fragments.
Default is ETD.}

\item{AsDataframe}{A boolean to indicate whether you would like the values used
to determine the ScanNumbers included or not. Default is FALSE.}
}
\description{
Return a vector of scan numbers based on a Precursor Mass, mass tolerance, and
isolation window filtering
}
\examples{
\dontrun{

# Generate the ScanMetadata object 
ScanMetadata <- pspecterlib::get_scan_metadata(
    MSPath = "/Users/degn400/Desktop/IsoForma_Test/Sorghum-Histone0622162L11.mzML"
)

# Run function
pull_scan_numbers(Sequence = "M.SGRGKGGKGLGKGGAKRHRKVLRDNIQGITKPAIRRLARRGGVKRISGLIYEETRGVLKIFLENVIRDAVTYTEHARRKTVTAMDVVYALKRQGRTLYGFGG",
                  Modification = "Acetyl,X(1*,5,8,12,16,20)[2]",
                  ScanMetdata = ScanMetadata,
                  PrecursorMass = 11355.38,
                  RTStart = 30,
                  RTEnd = 60)

}
}
