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

rbind(
  IsoForma111[[4]] %>% dplyr::mutate(Type = "IsoForma", Ratio = "1:1:1"),
  IsoForma111[[4]] %>% dplyr::mutate(Type = "Truth", Ratio = "1:1:1", Proportion = c(1/3, 1/3, 1/3), LowerCI = NA, UpperCI = NA),
  IsoForma148[[4]] %>% dplyr::mutate(Type = "IsoForma", Ratio = "1:4:8"),
  IsoForma148[[4]] %>% dplyr::mutate(Type = "Truth", Ratio = "1:4:8", Proportion = c(1/13, 4/13, 8/13), LowerCI = NA, UpperCI = NA),
  IsoForma418[[4]] %>% dplyr::mutate(Type = "IsoForma", Ratio = "4:1:8"),
  IsoForma418[[4]] %>% dplyr::mutate(Type = "Truth", Ratio = "4:1:8", Proportion = c(4/13, 1/13, 8/13), LowerCI = NA, UpperCI = NA),
  IsoForma814[[4]] %>% dplyr::mutate(Type = "IsoForma", Ratio = "8:1:4"),
  IsoForma814[[4]] %>% dplyr::mutate(Type = "Truth", Ratio = "8:1:4", Proportion = c(8/13, 1/13, 4/13), LowerCI = NA, UpperCI = NA)
) %>%
  dplyr::mutate(Ratio = factor(Ratio, levels = c("1:1:1", "1:4:8", "4:1:8", "8:1:4"))) %>%
  dplyr::mutate(Type = factor(Type, levels = c("Brunner et. al", "IsoForma", "Truth"))) %>%
  ggplot(aes(x = Modification, y = Proportion, fill = Type)) + 
   geom_bar(stat = "identity", position = "dodge", color = "black") + 
   scale_fill_manual(values = c("steelblue", "forestgreen")) + theme_bw() +
  ggplot2::geom_errorbar(ggplot2::aes(ymin = LowerCI, ymax = UpperCI), width = 0.2, position = ggplot2::position_dodge(.9)) +
  facet_wrap(.~Ratio)




