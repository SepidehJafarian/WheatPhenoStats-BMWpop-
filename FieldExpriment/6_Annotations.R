# Analysis of annotated genes at peak marker positions

rm(list = ls())

setwd("C:/Users/SepidehJafarian/Desktop/bmwpopulation-mudified-2019 (copy)/Interval_mapping")

library(tidyverse)
library(readxl)
library(biomaRt)
library(dplyr)

# Load QTL peak markers and support intervals --------------------
IM_table <- read.table("Results/IM_table.csv", sep = ",", header = TRUE, stringsAsFactors = FALSE)
IM_table$QTL[grep("TG0010", IM_table$QTL)] <- "TG0010a"
IM_table$QTL[grep("TG0011", IM_table$QTL)] <- "TG0011a"
IM_table$NextMrk[grep("TG0010", IM_table$NextMrk)] <- "TG0010a"
IM_table$NextMrk[grep("TG0011", IM_table$NextMrk)] <- "TG0011a"
IM_table$Pos <- trunc(IM_table$Pos * 100) / 100
IM_table <- IM_table[!IM_table$QTL %in% "No QTL detected", c(1:7, 10)]


# Load genetic and physical map --------------------
map <- read_excel("../Linkage_map/Genetic_map_5436.XLSX", "BMWpop_Map", skip = 1)
map$Marker[grep("TG0010", map$Marker)] <- "TG0010a"
map$Marker[grep("TG0011", map$Marker)] <- "TG0011a"
blast <- read.table("../Genotypes/TG20k/TAS16024_20k_blastn_IWGSC_WGA_v1.0.csv",
                    sep = ";", skip = 1, header = TRUE, stringsAsFactors = FALSE)
blast <- blast[, c("Query_id", "Subject_id", "pos", "e.value")]
colnames(blast) <- c("marker", "chr", "pos", "e")
blast$marker <- sub("-", ".", blast$marker)
blast_mapped <- blast[paste(blast$marker, sub("chr", "", blast$chr)) %in% # blast data to anchor mapped markers
                        paste(map$Marker, map$Chromosome), ]
map <- left_join(map, blast_mapped[, c("marker", "pos")], by = c("Marker" = "marker"))
colnames(map)[colnames(map) %in% c("Position", "pos")] <- c("Pos_genetic", "Pos_physical")
map <- as.data.frame(map)
blast <- blast %>% group_by(marker) %>% top_n(n = -1, wt = e) %>% as.data.frame() # blast data to enrich mapped markers


# Load automatically annotated genes --------------------
# Note: You should use the current RefSeq version, not v1.1.
genes <- read.table("C:/Users/SepidehJafarian/Desktop/bmwpopulation-mudified-2019 (copy)/Interval_mapping/IWGSC_v1.1_HC_20170706.gff3", sep = "\t")
colnames(genes) <- c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute")
genes <- genes %>% mutate_if(is.factor, as.character)
genes <- genes[genes$feature %in% "gene", ]


# Create tables with annotated genes --------------------
unlink("Results/Annotations/*")
for (i in 1:nrow(IM_table)) {
  trait <- IM_table$Trait[i]
  qtl <- IM_table$QTL[i]
  qtl_unique <- IM_table$QTL_unique[i]
  chr <- IM_table$Chr[i]
  pos <- IM_table$Pos[i]
  
  # Define mapped peak markers
  map_temp <- map[map$Chromosome %in% chr, ]
  map_temp$dist <- map_temp$Pos_genetic - pos
  map_temp <- map_temp[order(abs(map_temp$dist)), ]
  
  if (any(map_temp$dist == 0)) { 
    nextmrk <- map_temp[map_temp$dist == 0, "Marker"] 
  } else {
    if(abs(map_temp$dist)[1] <= 5) {
      nextpos <- abs(map_temp$dist)[1]
      nextmrk <- map_temp[abs(map_temp$dist) == nextpos, "Marker"]
    } else {
      nextmrk <- NA
    }
  }
  
  # Search and characterize genes at peak markers +/- 5 Mb 
  pos <- map[map$Marker %in% nextmrk, "Pos_physical"]
  pos_min <- min(pos, na.rm = TRUE) - 5e+06
  pos_max <- max(pos, na.rm = TRUE) + 5e+06
  genes_qtl <- genes[genes$seqname %in% paste0("chr", chr) &
                   genes$end >= pos_min &
                   genes$start <= pos_max, ]
  if (nrow(genes_qtl) >= 1) { 
    genes_qtl$gene_id <- str_extract(genes_qtl$attribute, "(?<=ID=)[[:alnum:]]*")
    # Add gene description and GO terms to candidate genes using Ensembl BioMart
   # Sys.setenv(http_proxy = "www.proxy.bybn.de:80") # you may omit or change this proxy setting
    # listAttributes(mart = useDataset("taestivum_eg_gene", useMart("plants_mart", host = "plants.ensembl.org")))
    bm_query <- getBM(attributes = c("ensembl_gene_id", "description", "name_1006"),
                      filters = "ensembl_gene_id",
                      values = genes_qtl$gene_id,
                      uniqueRows = TRUE,
                      mart = useDataset("taestivum_eg_gene", useMart("plants_mart", host = "https://plants.ensembl.org")))
    bm_query <- bm_query %>%
      mutate(ensembl_gene_id = factor(ensembl_gene_id, levels = unique(ensembl_gene_id)),
             go_term = na_if(name_1006, "")) %>%
      group_by(ensembl_gene_id) %>%
      summarise(description = unique(description),
                go_term = paste(na.omit(unique(go_term)), collapse="; ")) %>%
      mutate(ensembl_gene_id = as.character(ensembl_gene_id))
    genes_qtl <- left_join(genes_qtl, bm_query, by = c("gene_id" = "ensembl_gene_id"))
    write.table(genes_qtl, paste0("Results/Annotations/Annot_", trait, "_", chr, "_", qtl, ".csv"), sep = ";", row.names = FALSE, na = "")
  }
}





