# Plots illustrating QTL effects

rm(list = ls())

setwd("C:/Users/msmg/Desktop/BMWpop/Interval_mapping/")

library(tidyverse)

load("BMWpop.Rdata")

IM_table <- read.table("Results/IM_table.csv", sep = ";", header = TRUE)



# Create long table for effect plots --------------------
IM_table_long <- gather(data = IM_table, key = key, value = value, Mean_1:Tukey_8) %>% 
  extract(col = key, into = c("Statistic", "Allele"), regex = "(.*)_([[:digit:]])") %>%
  spread(key = Statistic, value = value) 

IM_table_long$N <- rep(NA, nrow(IM_table_long))
for (i in 1:nrow(IM_table_long)) {
  estfnd <- BMWpop$estfnd[[as.character(IM_table_long$Chr[i])]]
  IM_table_long$N[i] <- sum(estfnd[, colnames(estfnd) %in% IM_table_long$QTL[i]] 
                            %in% IM_table_long$Allele[i])
}

IM_table_long <- IM_table_long %>%
  mutate(Allele, Allele = case_when(Allele == 1 ~ "Event", 
                                    Allele == 2 ~ "BAYP4535",
                                    Allele == 3 ~ "Ambition",
                                    Allele == 4 ~ "Firl3565",
                                    Allele == 5 ~ "Format",
                                    Allele == 6 ~ "Potenzial",
                                    Allele == 7 ~ "Bussard",
                                    Allele == 8 ~ "Julius")) %>% 
  mutate(Trait = factor(Trait, levels = unique(Trait))) %>%
  mutate(QTL = factor(QTL, levels = unique(QTL))) %>%
  mutate(Allele = factor(Allele, levels = unique(Allele))) %>%
  group_by(Trait, QTL) %>%
  mutate(Mean = as.numeric(Mean), Mean_CT = as.numeric(Mean_CT), Mean_CT_SD = as.numeric(Mean_CT_SD), SE = as.numeric(SE)) %>%
  mutate(Cols_bin = case_when(Mean_CT < 0 ~ 0, Mean_CT > 0 ~ 1)) %>% 
  mutate(Cols_cont = (Mean - min(Mean, na.rm = TRUE)) / diff(range(Mean, na.rm = TRUE))) %>%
  arrange(.by_group = TRUE)


# QTL-wise barplots for comparison of means --------------------
unlink("Results/Effects/QTL-wise/*")
QTL_unique <- unique(IM_table$QTL_unique[!is.na(IM_table$QTL_unique)])
for (i in 1:length(QTL_unique)) {
  IM_table_long_sel <- IM_table_long[as.character(IM_table_long$QTL_unique) %in% QTL_unique[i], ]
  if (any(!is.na(IM_table_long_sel$Mean))) {
    tiff(paste0("Results/Effects/QTL-wise/Effects_", QTL_unique[i], ".tif"), width = 2200,
         height = length(unique(IM_table_long_sel$Trait)) * 700, res = 350, compression = "lzw")
    print(ggplot(IM_table_long_sel, aes(Allele, Mean_CT, fill = Cols_bin, label = Tukey)) +
            geom_bar(stat="identity", width = 0.6) +
            facet_grid(rows = vars(Trait, QTL), scales = "free", switch = "y") +
            scale_fill_gradientn(colors = c("magenta", "green")) +
            geom_errorbar(aes(ymin = Mean_CT - 1.96 * SE, ymax = Mean_CT + 1.96 * SE), width = 0.2) +
            geom_text(hjust = - 0.5) +
            annotate("text", x = 0.8, y = Inf, label = 'italic("N =")', parse = TRUE, vjust = "inward", size = 2) +
            geom_text(aes(Allele, Inf , label = N), vjust = "inward", size = 2) +
            geom_hline(yintercept = 0) + 
            theme(axis.text.x = element_text(angle = 50, hjust = 1)) +
            theme(legend.position = "none") +
            ylab(""))
    dev.off()
  }
  
}



# Trait-wise barplots for comparison of means --------------------
unlink("Results/Effects/Trait-wise/*")
trait_sel <- unique(gsub("_[[:alnum:]]*", "",  IM_table$Trait))
for (i in 1:length(trait_sel)) {
  IM_table_long_sel <- IM_table_long[grepl(paste0("^", trait_sel[i]), IM_table_long$Trait), ]
  tiff(paste0("Results/Effects/Trait-wise/Effects_", trait_sel[i], ".tif"), width = 2200,
       height = length(unique(paste(IM_table_long_sel$Trait, IM_table_long_sel$QTL))) * 700, res = 350, compression = "lzw")
  print(ggplot(IM_table_long_sel, aes(Allele, Mean_CT, fill = Cols_bin, label = Tukey)) +
          geom_bar(stat="identity", width = 0.6) +
          facet_grid(rows = vars(Trait, QTL, QTL_unique), switch = "y") +
          scale_fill_gradientn(colors = c("magenta", "green")) +
          geom_errorbar(aes(ymin = Mean_CT - 1.96 * SE, ymax = Mean_CT + 1.96 * SE), width = 0.2) +
          geom_text(hjust = - 0.5) +
          annotate("text", x = 0.8, y = Inf, label = 'italic("N =")', parse = TRUE, vjust = "inward", size = 2) +
          geom_text(aes(Allele, Inf , label = N), vjust = "inward", size = 2) +
          geom_hline(yintercept = 0) + 
          theme(axis.text.x = element_text(angle = 50, hjust = 1)) +
          theme(legend.position = "none") +
          ylab(""))
  dev.off()
}
