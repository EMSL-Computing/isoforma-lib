# isoforma-lib

R package for the relative quantification of proteoform positional isomers (QPI) from MS2 scan data

<img src="IsoForma_Logo.png" alt="IsoForma_Logo" style="width:250px"/>

# How to install 
`devtools::install_github("EMSL-Computing/isoforma-lib")`

To get started, read our [vignette](https://emsl-computing.github.io/isoforma-lib/)

You may need to install pspecterlib separately, though it should auto-install
with the above command. If not, try the command below. Check that both packages 
install with `library(isoforma-lib)` and `library(pspecterlib)`. If one of the 
packages does not install, 

`devtools::install_github("EMSL-Computing/pspecterlib")`

# Installation Checks

The package has been succesfully installed on an Apple M1 Max using 13.6.1 and a Windows PC running Windows 11. **We encourage users to please report any installation issues on the Github Issues page, and we are more than happy to help!**

# Usage (Windows)

1. Download and install R: https://ftp.osuosl.org/pub/cran/
2. Download, install, and open RStudio (Free version): https://www.rstudio.com/products/rstudio/download/
3. If re-installing, remove any previous versions of pspecter
	remove.packages("pspecterlib")
4. Clone/download IsoForma-paper repository

# Reference

If you use IsoForma or any portions of this code please cite: Degnan et al. "IsoForma: An R Package for Quantifying and Visualizing Positional Isomers in Top-Down LC-MS/MS Data". Submitted to Journal of Proteome Research.

# Modifications

For a list of acceptable modifications, see the "Modification" column of [Glossary](https://github.com/EMSL-Computing/isoforma-lib/blob/david_develop/inst/extdata/Unimod_v20220602.csv). If there is a modification needed that is not in this list, please report the modification to our issues page. 
