# LD analysis in QTL support intervals

rm(list = ls())

setwd("C:/Users/msmg/Desktop/BMWpop/Interval_mapping/")

library(tidyverse)
library(readxl)
library(genetics)
library(LDheatmap)
library(grid)


# Load QTL peak markers and support intervals --------------------
IM_table <- read.table("Results/IM_table.csv", sep = ";", header = TRUE, stringsAsFactors = FALSE)
IM_table$QTL[grep("TG0010", IM_table$QTL)] <- "TG0010a"
IM_table$QTL[grep("TG0011", IM_table$QTL)] <- "TG0011a"
IM_table$NextMrk[grep("TG0010", IM_table$NextMrk)] <- "TG0010a"
IM_table$NextMrk[grep("TG0011", IM_table$NextMrk)] <- "TG0011a"
IM_table <- IM_table[!IM_table$QTL %in% "No QTL detected", c(1:7, 10)]
IM_table$Pos <- trunc(IM_table$Pos * 100) / 100


# Load genetic and physical map --------------------
map <- read_excel("../Linkage_map/Genetic_map_5436.XLSX", "BMWpop_Map", skip = 1)
map$Marker[grep("TG0010", map$Marker)] <- "TG0010a"
map$Marker[grep("TG0011", map$Marker)] <- "TG0011a"
blast <- read.table("../Genotypes/TG20k/TAS16024_20k_blastn_IWGSC_WGA_v1.0.csv",
                    sep = ";", skip = 1, header = TRUE, stringsAsFactors = FALSE)
blast <- blast[, c("Query_id", "Subject_id", "pos")]
colnames(blast) <- c("marker", "chr", "pos")
blast$marker <- sub("-", ".", blast$marker)
blast <- blast[paste(blast$marker, sub("chr", "", blast$chr)) %in%
                 paste(map$Marker, map$Chromosome), ]
map <- left_join(map, blast[, c("marker", "pos")], by = c("Marker" = "marker"))
colnames(map)[colnames(map) %in% c("Position", "pos")] <- c("Pos_genetic", "Pos_physical")
map <- as.data.frame(map)


# Load genotypic data --------------------
# SNPs were filtered with nmiss = 0.1, maf = 0.01 and mapped markers were kept
geno <- read.table("../Genotypes/TG20k/Genotypic_data_394_lines_10842_SNPs.csv", sep = ";", header = TRUE, row.names = 1)
dim(geno) # 394 10842


# LD heatmaps for all QTL --------------------
unlink("Results/LD/*")
for (i in 1:nrow(IM_table)) {
  trait <- IM_table[i, "Trait"]
  qtl <- IM_table[i, "QTL"]
  chr <- IM_table[i, "Chr"]
  peakmrk <- map[map$Pos_genetic %in% IM_table$Pos[i] &
                   map$Chromosome %in% chr, "Marker"]
  si_lower <- IM_table[i, "SI_lower"]
  si_upper <- IM_table[i, "SI_upper"]
  
  map_si <- map[map$Chromosome %in% chr &
                  map$Pos_genetic >= si_lower &
                  map$Pos_genetic <= si_upper, ]
  map_si <- map_si[order(map_si$Pos_physical), ]
  miss_pos <- sum(is.na(map_si$Pos_physical))
  map_si <- map_si[!is.na(map_si$Pos_physical), ]
  si_size <- max(map_si$Pos_physical) - min(map_si$Pos_physical)
  if(si_size > 0) {
    
    geno_si <- geno[, map_si$Marker]
    geno_si <- makeGenotypes(geno_si)
    
    ld <- LD(geno_si)[["R^2"]]
    ld <- LDheatmap(ld, genetic.distance = map_si$Pos_physical, SNP.name = peakmrk, geneMapLabelX = 2)
    ld <- editGrob(ld$LDheatmapGrob, gPath("geneMap", "SNPnames"), gp = gpar(cex = 0.5))
    pdf(paste0("Results/LD/LD_", trait, "_", chr, "_", qtl, ".pdf"))
    grid.draw(ld)
    grid.text(paste("Interval size:", round(si_size / 1e6, 2) , "Mb"),
              x = 0.8, y = 0.25, gp = gpar(fontsize = 16))
    if(miss_pos > 0) {
      grid.text(paste(miss_pos, "marker(s) not shown"), x = 0.8, y = 0.2, gp = gpar(fontsize = 16))
    }
    dev.off()
  }
}
