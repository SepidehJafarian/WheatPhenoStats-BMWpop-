# Composite interval mapping using mpMap

rm(list = ls())

setwd("~/Desktop/bmwpopulation-2019/Interval_mapping/")

library(mpMap)
library(lme4)
library(aods3)
library(qtl)
library(VPdtw)
library(tidyverse)
library(data.table)
library(multcompView)

source("Functions/mpIM.R")
source("Functions/plot.mpqtl.R")
source("Functions/fit.mpqtl.R")
 
load("BMWpop.Rdata")

# Loop for composite interval mapping over all traits --------------------
# choosing only traits using our population.
# applying threshhold.
#summary: to store all the effects on it

head(BMWpop$pheno)
traits <- colnames(BMWpop$pheno)
View(traits)
thr <- 4.256106
IM_results <- list()
IM_summary <- list()

#(SIM)-sinmple interval mapping which is pre-requisition for composit interval mapping.
#(CIM)- composite interval mapping.
#cat(): To returne character string in a readable format.
#attr():To access and extract single attribute of our data.
#findqtl:To detect QTL peakes.given output from scan of a chromosom.locates local maxima exceeding a #significance threshold in a QTL profile.
#mpIM: obj: object of class mpcross
#ncove: number of markers covariates to search. as many as possible=default(forward/backward selection)
#responsname: optinal input of response name to look for object$pheno
#step:step size at which to compute QTL profi:
#dwindow:Window over which to smooth p-values - default is five markers

 for (i in 1:length(traits)){
  cat("------ QTL analysis for", traits[i], "----------\n")
  # First perform simple interval mapping
  QTL <- mpIM(object = BMWpop, ncov = 0, responsename = traits[i], step = 1)
  QTL <- findqtl(QTL, dwindow = 100, threshold = thr)
  nqtl <- attr(QTL$QTLresults$qtl, "nqtl")
  View(QTL)
  View(QTL$QTLresults$qtl)
  View(nqtl)

  # Then perform composite interval mapping using the number of QTL from SIM as cofactors
  # Then perform composite interval mapping using the number of QTL from SIM as cofactors
  #The difference between them is that [[ ]] is used to access a component in a
  #list or matrix whereas [ ] is used to access a single element in a matrix or
  #array
  QTL <- mpIM(object = BMWpop, ncov = nqtl, responsename = traits[i], step = 1)
  QTL <- findqtl(QTL, dwindow = 100, threshold = thr)
  IM_results[[i]] <- unclass(QTL)
  head(QTL)
  View(QTL)
  ##The difference between them is that [[ ]] is used to access a component in a
  #list or matrix whereas [ ] is used to access a single element in a matrix or
  #array
  names(IM_results)[[i]] <- traits[i]
  head(IM_results)
  IM_summary[[i]] <- fit.mpqtl(QTL)
  names(IM_summary)[[i]] <- traits[i]
 }






  # Create table summarizing results for all detected QTL --------------------

IM_table <- sapply(IM_summary, function(x) x[1])
names(IM_table) <- traits

IM_table <- bind_rows(IM_table, .id = "Trait")
View(IM_table)
IM_table <- mutate_if(IM_table, is.factor, as.character) 
View(chr)
chr <- names(BMWpop$map) # unique QTL identification
IM_table_QTL <- vector("list", length(chr))
View(IM_table_QTL)
for (i in 1:length(chr)) {
  IM_table_chr <- IM_table[IM_table$Chr %in% chr[i] &
                             !IM_table$QTL %in% "No QTL detected" , ]
  setDT(IM_table_chr)
  setkey(IM_table_chr, SI_lower, SI_upper)
  IM_table_chr[, row_id :=1:nrow(IM_table_chr)]
  temp <- foverlaps(IM_table_chr, IM_table_chr)
  temp[, `:=`(c("SI_lower", "SI_upper"), list(min(SI_lower, i.SI_lower), max(SI_upper, i.SI_upper))), by = row_id]
  temp[, `:=`(c("SI_lower", "SI_upper"), list(min(SI_lower, i.SI_lower), max(SI_upper, i.SI_upper))), by = i.row_id]
  temp <- temp[, list(QTL_unique = .GRP, row_id = unique(c(row_id, i.row_id))), by = .(SI_lower, SI_upper)][, .(row_id, QTL_unique)]
  setkey(IM_table_chr, row_id)
  setkey(temp, row_id)
  IM_table_QTL[[i]] <- temp[IM_table_chr]  
}

IM_table_QTL <- do.call("rbind", IM_table_QTL)
IM_table_QTL <- IM_table_QTL[!duplicated(paste(IM_table_QTL$Chr, IM_table_QTL$row_id))]
IM_table_QTL$QTL_unique <- paste0("Q_", IM_table_QTL$Chr, "_", IM_table_QTL$QTL_unique)
IM_table <- merge(IM_table, IM_table_QTL[, c("Trait", "QTL", "QTL_unique", "Chr")], by = c("Trait", "QTL", "Chr"), all.x = TRUE, sort = FALSE)
IM_table <- IM_table[, c(1, 2, ncol(IM_table), 3:(ncol(IM_table) - 1))]


mean_raw <- IM_table[, seq(13, 27, 2)] # add intercept to all effects
for (i in 1:nrow(mean_raw)) {
  if (any(!is.na(mean_raw[i, ]))) {
    ref <- which(!is.na(mean_raw[i, ]))[1]
    mean_raw[i, - ref] <- mean_raw[i, - ref] + mean_raw[i, ref]
  }
}
colnames(mean_raw) <- sub("Effect", "Mean", colnames(mean_raw))
head(mean_raw)

mean_ct <- mean_raw # centering raw means by substracting mean of raw means 
for (i in 1:nrow(mean_ct)) {
  mean_ct[i, ] <- mean_ct[i, ] - mean(as.numeric(mean_ct[i, ]), na.rm = TRUE)
}
colnames(mean_ct) <- sub("Mean", "Mean_CT", colnames(mean_ct))
head(mean_ct)
mean_ct_sd <- mean_ct # divide centered means by standard deviation
for (i in 1:nrow(mean_ct_sd)) {
  mean_ct_sd[i, ] <- mean_ct_sd[i, ] / sd(BMWpop$pheno[, IM_table$Trait[i]], na.rm = TRUE)
}
colnames(mean_ct_sd) <- sub("Mean_CT", "Mean_CT_SD", colnames(mean_ct_sd))
head(mean_ct_sd)
Effect_MAX <- apply(mean_ct, 1, max, na.rm = TRUE) - apply(mean_ct, 1, min, na.rm = TRUE) # maximum effect per QTL
Effect_MAX_SD <- apply(mean_ct_sd, 1, max, na.rm = TRUE) - apply(mean_ct_sd, 1, min, na.rm = TRUE)

tuk <- matrix(nrow = nrow(IM_table), ncol = 8) # pairwise comparison of means
colnames(tuk) <- paste("Tukey", 1:8, sep = "_")
for (i in 1:nrow(IM_table)) {
  if (!IM_table[i, "QTL"] %in% "No QTL detected" &
      !is.null(IM_summary[[IM_table[i, "Trait"]]][[IM_table[i, "QTL"]]])) {
    mod <- aov(IM_summary[[IM_table[i, "Trait"]]][[IM_table[i, "QTL"]]])
    tuk_let <- multcompLetters4(mod, TukeyHSD(mod))[[1]][["Letters"]]
    tuk[i, as.numeric(names(tuk_let))] <- tuk_let
  }
}

IM_table <- cbind(IM_table[, 1:12], Effect_MAX, Effect_MAX_SD, mean_raw, mean_ct, mean_ct_sd, IM_table[, seq(14, 28, 2)], tuk)

write.table(IM_table, file = "Results/IM_table.csv", sep = ";", row.names = FALSE)


# p-value curves --------------------

# Combine p-value curves for each unique QTL
unlink("Results/p_values/QTL-wise/*")
QTL_unique <- unique(IM_table$QTL_unique[!is.na(IM_table$QTL_unique)])
for (i in 1:length(QTL_unique)) {
  trait_sel <- IM_table[IM_table$QTL_unique %in% QTL_unique[i], "Trait"]
  IM_results_sel <- IM_results[names(IM_results) %in% trait_sel]
  height <- (length(IM_results_sel) + 1) * 100
  tiff(paste0("Results/p_values/QTL-wise/p_values_", QTL_unique[i], ".tif"), width = 1500, height = height, compression = "lzw")
  par(mfrow = c(length(IM_results_sel) + 1, 1), xaxt = 'n', mar = c(0, 4, 0, 10), cex = 1)
  for (i in 1:length(IM_results_sel)) {
    plot.mpqtl(IM_results_sel[[i]], incl.markers = FALSE, xlab = "", ylab = "", axes = FALSE)
    axis(side = 2, las = 1, 
         at = seq(0, par("yaxp")[2] - 1, 2),
         labels = as.character(seq(0, par("yaxp")[2] - 1, 2)))
    mtext(text = names(IM_results_sel[i]), side = 4, las = 1, line = 1, cex = 2)
    abline(h = thr, col = "red")
    box(lwd = 3)
  }
  par(mar = c(6.9, 4, 0, 10), xaxt = 's', yaxt = 'n', cex.axis = 2)
  plot.mpqtl(IM_results_sel[[i]], incl.markers = FALSE, xlab = "", ylab = "")
  dev.off()
}

# Combine p-value curves for each trait
unlink("Results/p_values/Trait-wise/*")
trait_sel <- unique(gsub("_[[:alnum:]]*", "", names(IM_results)))
for (i in 1:length(trait_sel)) {
  IM_results_sel <- IM_results[grepl(paste0("^", trait_sel[i]), names(IM_results))]
  height <- (length(IM_results_sel) + 1) * 100
  tiff(paste0("Results/p_values/Trait-wise/p_values_", trait_sel[i], ".tif"), width = 1500, height = height, compression = "lzw")
  par(mfrow = c(length(IM_results_sel) + 1, 1), xaxt = 'n', mar = c(0, 4, 0, 10), cex = 1)
  for (i in 1:length(IM_results_sel)) {
    plot.mpqtl(IM_results_sel[[i]], incl.markers = FALSE, xlab = "", ylab = "", axes = FALSE)
    axis(side = 2, las = 1, 
         at = seq(0, par("yaxp")[2] - 1, 2),
         labels = as.character(seq(0, par("yaxp")[2] - 1, 2)))
    mtext(text = names(IM_results_sel[i]), side = 4, las = 1, line = 1, cex = 2)
    abline(h = thr, col = "red")
    box(lwd = 3)
  }
  par(mar = c(6.9, 4, 0, 10), xaxt = 's', yaxt = 'n', cex.axis = 2)
  plot.mpqtl(IM_results_sel[[i]], incl.markers = FALSE, xlab = "", ylab = "")
  dev.off()
}


# Table with descriptive trait statistics --------------------
trait_stat <- data.frame(Trait = traits,
                         NQTL = rep(NA, length(traits)),
                         R2_min = rep(NA, length(traits)),
                         R2_max = rep(NA, length(traits)),
                         Effect_MAX = rep(NA, length(traits)))
for (i in 1:length(traits)) {
  trait_stat$NQTL[i] <- sum(IM_table$Trait %in% traits[i] & !IM_table$QTL %in% "No QTL detected")
  if (trait_stat$NQTL[i] > 0 & !all(is.na(IM_table[IM_table$Trait %in% traits[i], "R2"]))) {
    trait_stat$R2_min[i] <- min(IM_table[IM_table$Trait %in% traits[i], "R2"], na.rm = TRUE)
    trait_stat$R2_max[i] <- max(IM_table[IM_table$Trait %in% traits[i], "R2"], na.rm = TRUE)
    trait_stat$Effect_MAX[i] <- max(IM_table[IM_table$Trait %in% traits[i], "Effect_MAX"], na.rm = TRUE)
  }
}
write.table(trait_stat, file = "Results/Trait_statistics.csv", sep = ";", row.names = FALSE)
