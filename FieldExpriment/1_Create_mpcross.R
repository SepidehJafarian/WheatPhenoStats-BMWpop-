# Create mpcross object

rm(list = ls())

setwd("C:/Users/SepidehJafarian/Desktop/bmwpopulation-mudified-2019 (copy)/Interval_mapping")


library(mpMap) 1.14 -> 2.0.2
library(lme4) 1.1-35.5
library(aods3) 0.4-1.2
library(qtl)
library(VPdtw) 2.2.1


MR19 <- read.table("../Phenotypes/Analysis/2018_2019/MR19_means.csv", sep = ";", header = TRUE, row.names = 1)
colnames(MR19) <- paste(colnames(MR19), "MR19", sep = "_")
pheno <- MR19

IDs <- read.table("../IDs_BMWpop.txt", sep = ";", header = TRUE)
View (IDs)
pheno <- merge(pheno, IDs, by.x = 0, by.y = "Line", all.y = TRUE)
pheno <- pheno[, c(ncol(pheno), 2:(ncol(pheno) - 1))]
rownames(pheno) <- pheno$ID
pheno <- pheno[, - 1]
write.table(pheno, "../Phenotypes/Analysis/Phenotypes_mpmap.txt", sep = "\t", col.names = NA, quote = FALSE)

BMWpop <- read.mpcross(founderfile = "../Genotypes/Genotypes_founders_2804.txt",
                       finalfile = "../Genotypes/Genotypes_BMWpop_2804.txt",
                       pedfile = "C:/Users/SepidehJafarian/Desktop/bmwpopulation-mudified-2019 (copy)/Pedigree.txt",
                       mapfile =  "../Linkage_map/Genetic_map_2804.txt",
                       phenofile = "../Phenotypes/Analysis/Phenotypes_mpmap.txt")

identical(rownames(BMWpop$pheno), rownames(BMWpop$finals))
rownames(BMWpop$pheno)
rownames(BMWpop$finals)
print(dim(BMWpop$founders))
print(dim(BMWpop$finals))
print(dim(BMWpop$pheno))
print(dim(BMWpop$pedigree))

BMWpop <- mpprob(BMWpop, program = "qtl", mapfx = "haldane", threshold = 0.5, step = 1)
save(BMWpop, file = "BMWpop.Rdata")
