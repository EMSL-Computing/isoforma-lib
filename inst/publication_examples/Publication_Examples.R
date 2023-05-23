## Brunner Benchmark ## 

# Designate modifications
MultipleMods <- pspecterlib::multiple_modifications(
  Sequence = "LQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG",
  Modification = "6.018427,V(17,26,70)[1]",
  ReturnUnmodified = TRUE
)

ion <- "c"

# 1:1:1
PeakSum111 <- sum_ms2_spectra(
  PeakDataList = list(
    readRDS(system.file("extdata", "PeakData_1to1to1_1.RDS", package = "isoforma")),
    readRDS(system.file("extdata", "PeakData_1to1to1_2.RDS", package = "isoforma")),
    readRDS(system.file("extdata", "PeakData_1to1to1_3.RDS", package = "isoforma"))
  )
)
IsoForma111 <- isoforma_pipeline(Sequences = MultipleMods, SummedSpectra = PeakSum111,
  PrecursorCharge = 10, ActivationMethod = "ETD", IonGroup = ion, CorrelationScore = 0,
  IsotopeAlgorithm = "Rdisop", Messages = TRUE)

# 1:4:8
PeakSum148 <- sum_ms2_spectra(
  PeakDataList = list(
    readRDS(system.file("extdata", "PeakData_1to4to8_1.RDS", package = "isoforma")),
    readRDS(system.file("extdata", "PeakData_1to4to8_2.RDS", package = "isoforma")),
    readRDS(system.file("extdata", "PeakData_1to4to8_3.RDS", package = "isoforma"))
  )
)
IsoForma148 <- isoforma_pipeline(Sequences = MultipleMods, SummedSpectra = PeakSum148,
                                 PrecursorCharge = 11, ActivationMethod = "ETD", IonGroup = ion, CorrelationScore = 0,
                                 IsotopeAlgorithm = "Rdisop", Messages = TRUE)

# 4:1:8
PeakSum418 <- sum_ms2_spectra(
  PeakDataList = list(
    readRDS(system.file("extdata", "PeakData_4to1to8_1.RDS", package = "isoforma")),
    readRDS(system.file("extdata", "PeakData_4to1to8_2.RDS", package = "isoforma")),
    readRDS(system.file("extdata", "PeakData_4to1to8_3.RDS", package = "isoforma"))
  )
)
IsoForma418 <- isoforma_pipeline(Sequences = MultipleMods, SummedSpectra = PeakSum418,
                                 PrecursorCharge = 11, ActivationMethod = "ETD", IonGroup = ion, CorrelationScore = 0,
                                 IsotopeAlgorithm = "Rdisop", Messages = TRUE)

# 8:1:4
PeakSum814 <- sum_ms2_spectra(
  PeakDataList = list(
    readRDS(system.file("extdata", "PeakData_8to1to4_1.RDS", package = "isoforma")),
    readRDS(system.file("extdata", "PeakData_8to1to4_2.RDS", package = "isoforma")),
    readRDS(system.file("extdata", "PeakData_8to1to4_3.RDS", package = "isoforma"))
  )
)
IsoForma814 <- isoforma_pipeline(Sequences = MultipleMods, SummedSpectra = PeakSum814,
                                 PrecursorCharge = 11, ActivationMethod = "ETD", IonGroup = ion, CorrelationScore = 0,
                                 IsotopeAlgorithm = "Rdisop", Messages = TRUE)

# Quick Visualize Results
library(patchwork)
library(ggplot2)
(IsoForma111[[5]] + ggtitle("1:1:1")) +
(IsoForma148[[5]] + ggtitle("1:4:8")) +
(IsoForma418[[5]] + ggtitle("4:1:8")) +
(IsoForma814[[5]] + ggtitle("8:1:4")) + patchwork::plot_annotation(title = paste("Ion Type:", ion))




