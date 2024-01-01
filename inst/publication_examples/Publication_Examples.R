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
  PrecursorCharge = 10, ActivationMethod = "ETD", IonGroup = ion, 
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
                                 PrecursorCharge = 11, ActivationMethod = "ETD", IonGroup = ion, 
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
                                 PrecursorCharge = 11, ActivationMethod = "ETD", IonGroup = ion, 
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
                                 PrecursorCharge = 11, ActivationMethod = "ETD", IonGroup = ion,
                                 IsotopeAlgorithm = "Rdisop", Messages = TRUE)

# Quick Visualize Results
library(patchwork)
library(ggplot2)

Brunner <- rbind(
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
  dplyr::mutate(Modification = ifelse(Modification == "6.018427@V17", "V17",
                                      ifelse(Modification == "6.018427@V26", "V26", "V70"))) %>%
  ggplot(aes(x = Modification, y = Proportion, fill = Type)) + ylim(c(0,0.7)) +
   geom_bar(stat = "identity", position = "dodge", color = "black") +
   scale_fill_manual(values = c("steelblue", "forestgreen")) + theme_bw() +
  ggplot2::geom_errorbar(ggplot2::aes(ymin = LowerCI, ymax = UpperCI), width = 0.2, position = ggplot2::position_dodge(.9)) +
  theme(legend.position = "bottom", plot.title = element_text(hjust = 0.5)) +
  ggtitle("Valine Dataset") +
  xlab("") +
  theme(plot.title = element_text(size = 22), axis.title.y = element_text(size = 20), 
        axis.text.x = element_text(size = 18), axis.text.y = element_text(size = 18),
        legend.text = element_text(size = 18), strip.text = element_text(size = 18)) + 
  facet_wrap(.~Ratio)

Brunner

## Pasavento Benchmark ## 

# Load raw mzML data 
xml_data <- pspecterlib::get_scan_metadata(MSPath = system.file("extdata", "Example.mzML", package = "isoforma"))

# Pull scan numbers
Sequence <- "SGRGKGGKGLGKGGAKRHRKVLRDNIQGITKPAIRRLARRGGVKRISGLIYEETRTVLKTFLENVIRDSVTYTEHARRKTVTAMDVVYALKRQGRTLYGFGG"
Modifications <- "Acetyl,X(1^,5,8,12,16)[2];Methyl,K(79)[1];Oxidation,M(84)[1]"
Modified_Sequences <- pspecterlib::multiple_modifications(Sequence, Modifications, ReturnUnmodified = TRUE)
Scan_Numbers <- pull_scan_numbers(Sequence = Modified_Sequences[5], ScanMetadata = xml_data, RTStart = 100, RTEnd = 110)

# Sum Peaks
Summed_Peaks <- sum_ms2_spectra(ScanMetadata = xml_data, ScanNumbers = Scan_Numbers)

# Run the pipeline
IsoForma2 <- isoforma_pipeline(
    Sequences = Modified_Sequences,
    SummedSpectra = Summed_Peaks,
    PrecursorCharge = 16, 
    ActivationMethod = "ECD",
    IonGroup = "c",
    IsotopeAlgorithm = "Rdisop" # Rdisop is preferred, is faster, and is more accurate, but it tends to crash on Windows
)

Pasavento <- IsoForma2[[4]] %>%
  dplyr::mutate(
    Residue = factor(c("K5", "K8", "K12", "K16"), levels = c("K5", "K8", "K12", "K16")),
    Expected = c(0.38, 0.08, 0.42, 0.12)
  ) %>%
  dplyr::select(Proportion, LowerCI, UpperCI, Residue, Expected) %>%
  tidyr::pivot_longer(c(Proportion, Expected)) %>%
  dplyr::rename(Type = name, Proportion = value) %>%
  dplyr::mutate(
    Type = ifelse(Type == "Proportion", "IsoForma", "Truth"),
    LowerCI = ifelse(Type == "Truth", NA, LowerCI),
    UpperCI = ifelse(Type == "Truth", NA, UpperCI)
  ) %>%
  ggplot(aes(x = Residue, y = Proportion, fill = Type)) + ylim(c(0,0.7)) +
   geom_bar(stat = "identity", position = "dodge", color = "black") +
   scale_fill_manual(values = c("steelblue", "forestgreen")) + theme_bw() +
   ggplot2::geom_errorbar(ggplot2::aes(ymin = LowerCI, ymax = UpperCI), width = 0.2, position = ggplot2::position_dodge(.9)) +
   theme(legend.position = "bottom", plot.title = element_text(hjust = 0.5)) +
    xlab("") +
  theme(plot.title = element_text(size = 22), axis.title.y = element_text(size = 20), 
        axis.text.x = element_text(size = 18), axis.text.y = element_text(size = 18),
        legend.text = element_text(size = 18), strip.text = element_text(size = 18)) + 
   ggtitle("Histone Dataset")

Fig1 <- Brunner + Pasavento + patchwork::plot_annotation(tag_levels = "A")
Fig1

ggsave("~/Downloads/IsoForma_Fig1_300dpi.png", dpi = 300)




