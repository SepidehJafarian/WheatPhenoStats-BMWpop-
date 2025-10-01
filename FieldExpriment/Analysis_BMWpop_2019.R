# Phenotypic analysis of trial MR19


setwd("~/Desktop/bmwpopulation-2019/Phenotypes/Analysis/2018_2019")
rm(list = ls())
library(tidyverse)
library(readxl)
library(lme4)
library(emmeans)
library(outliers)


source("../grubbs_loop.R")


# Loading and formatting --------------------

# Field book
MR19 <- read.table("../../Field_books/2018_2019/Field_book_BMWpop_MR19.csv", header = TRUE, sep = ";")
View(MR19)
# Heading date
MR19_HD <- read_excel("MR19_raw_data_agronomic.xlsx", "Boniturtabelle")
MR19_HD <- MR19_HD[1:860, c(3, 6)]
colnames(MR19_HD) <- c("id", "HD")
MR19_HD$plot <- as.numeric(sub("[[:space:]][[:alnum:]]*[[:punct:]]?[[:alnum:]]*", "", MR19_HD$id))
MR19_HD <- MR19_HD[, c("plot", "HD")]
for (i in 1:nrow(MR19_HD)){
  if(!is.na(MR19_HD$HD[i])) {
    if (MR19_HD$HD[i] < 29) {
      MR19_HD$HD[i] <- MR19_HD$HD[i] + 31
    }
  }
}
View(MR19_HD)
# Grain yield
MR19_GY <- read_excel("MR19_raw_data_GY.xlsx", "Ernte")
MR19_GY <- MR19_GY[, c(1, 2, 6)]
colnames(MR19_GY) <- c("RND", "plot", "GY")
head(MR19_GY)
View(MR19_GY)
# Thousand-grain weight, grain width, grain length, grain area and grain protein content
MR19_TGW_GPC <- read_excel("MR19_raw_data_TGW_GPC.xlsx")
print(MR19_TGW_GPC )
MR19_TGW_GPC <- MR19_TGW_GPC[, c(1, 2, 4, 5, 3, 7)]
print(MR19_TGW_GPC)
colnames(MR19_TGW_GPC) <- c("RND", "TGW", "GW", "GL", "GA", "GPC")
MR19_TGW_GPC <- MR19_TGW_GPC[-1, ]

MR19_TGW_GPC$RND[MR19_TGW_GPC$RND %in% "b3083"] <- 23083
head(MR19_TGW_GPC)
MR19_TGW_GPC$RND[MR19_TGW_GPC$RND %in% "b3237"] <- 23237
MR19_TGW_GPC$RND <- as.numeric(MR19_TGW_GPC$RND)
View(MR19_TGW_GPC)
head(MR19_TGW_GPC)

# Isotope discrimination measurements could be loaded here 
MR19_ISD <- read_excel("Morgenrot_2019_ISD.xlsx")
View(MR19_ISD)
MR19_ISD<- MR19_ISD[, c(1,9,16)]
view(MR19_ISD)
colnames(MR19_ISD)<- c("plot", "N_take", "Isotope")
head(MR19_ISD)
# Merge phenotypic data

MR19 <- left_join(MR19, MR19_HD, by = c("plots" = "plot"))
MR19 <- left_join(MR19, MR19_GY, by = c("plots" = "plot"))
MR19 <- left_join(MR19, MR19_TGW_GPC, by = c("RND" = "RND"))
MR19<-  left_join(MR19, MR19_ISD, by = c("plots"="plot"))
head(MR19)
str(MR19)
View(MR19)

MR19 <- MR19[, c(4, 5, 3, 7, 9:16)]
colnames(MR19)[1:2] <- c("line", "rep")
View(MR19)
MR19$line <- as.character(MR19$line)
MR19$rep <- as.factor(MR19$rep)
MR19$block <- as.factor(MR19$block)
MR19$GW <- as.numeric(MR19$GW)
MR19$GL <- as.numeric(MR19$GL)
MR19$GA <- as.numeric(MR19$GA)

str(MR19)
view(MR19)

# Rename parents and standards (some names were chosen poorly)
unique(MR19$line)
view(MR19$line)
MR19$line[grep("Ambition", MR19$line)] <- "Ambition"
MR19$line[grep("Me/De", MR19$line)] <- "BAYP4535"
MR19$line[grep("Bussard", MR19$line)] <- "Bussard"
MR19$line[grep("Event", MR19$line)] <- "Event"
MR19$line[grep("Amigo", MR19$line)] <- "Firl3565"
MR19$line[grep("Format", MR19$line)] <- "Format"
MR19$line[grep("Julius", MR19$line)] <- "Julius"
MR19$line[grep("Potenzial", MR19$line)] <- "Potenzial"
MR19$line[grep("Reform", MR19$line)] <- "Reform"

MR19$line <- as.factor(MR19$line) 
View(MR19)


# Discarding outliers --------------------

# pivot_longer(MR19, HD:GPC, names_to = "trait", values_to = "value")
#ggplot(., aes(value)) +
#geom_histogram() +
#facet_wrap(~trait, scale = "free")


MR19 %>%
  pivot_longer(HD:Isotope, names_to = "trait", values_to = "value") %>%
  ggplot(., aes(value)) +
  geom_histogram() +
  facet_wrap(~trait, scale = "free")


MR19[, 4:ncol(MR19)] <- apply(MR19[, 4:ncol(MR19)], 2, grubbs_loop)

MR19 %>%
  pivot_longer(HD:Isotope, names_to = "trait", values_to = "value") %>%
  ggplot(., aes(value)) +
  geom_histogram() +
  facet_wrap(~trait, scale = "free")
# Derive GPY, GPD and GYD --------------------
MR19$GPY <- MR19$GY * MR19$GPC / 100
MR19$GPD <- resid(lm(GPC ~ GY, data = MR19, na.action=na.exclude))
MR19$GYD <- resid(lm(GY ~ GPC, data = MR19, na.action=na.exclude))
head(MR19)
# Discard outliers for newly derived traits
MR19 %>%
  pivot_longer(HD:GYD, names_to = "trait", values_to = "value") %>%
  ggplot(., aes(value)) +
  geom_histogram() +
  facet_wrap(~trait, scale = "free")

MR19[, 4:ncol(MR19)] <- apply(MR19[, 4:ncol(MR19)], 2, grubbs_loop)

View(MR19)  


# Calculate repeatability (genetic contribution to observed variance, similar to heritability)  --------------------
env <- MR19
view(env)
traits <- colnames(env)[4:ncol(env)]
var <- matrix(nrow = 3, ncol = length(traits))
rownames(var) <- c("varG", "varR", "rep")
colnames(var) <- traits

for (i in 1:length(traits)) {
  m.var <- lmer(get(traits[i]) ~ (1|line) + rep + (1|block), data = env)
  print(m.var)
  var["varG", i] <- VarCorr(m.var)$line[[1]]  
  var["varR", i] <- attr(VarCorr(m.var), "sc")^2    
  var["rep", i] <- var["varG", i] / (var["varG", i] + var["varR", i] / length(unique(env$rep))) 
   

}
head(var)
finalvar<- round(var, digits = 2)
write.table (finalvar, "repeatability.csv")

#var mit Isotope by sepideh
#       HD    GY   TGW   GW   GL   GA  GPC N_take  Isotope  GPY   GPD  GYD
#varG 3.51 16.10 11.87 0.01 0.06 0.80 0.49   0.01    0.05  0.11  0.18  7.14
#varR 0.56  2.71  0.76 0.00 0.00 0.08 0.15   0.01    0.01  0.06  0.12  2.27
#rep  0.93  0.92  0.97 0.95 0.97 0.95 0.86   0.69    0.89  0.78  0.75  0.86


# round(var, 2)
#        HD    GY   TGW   GW   GL   GA  GPC  GPY  GPD  GYD
# varG 3.51 16.10 11.87 0.01 0.06 0.80 0.49 0.11 0.18 7.14
# varR 0.56  2.71  0.76 0.00 0.00 0.08 0.15 0.06 0.12 2.27
# rep  0.93  0.92  0.97 0.95 0.97 0.95 0.86 0.78 0.75 0.86

# Calculate adjusted mean --------------------

env <- MR19
View(MR19)
traits <- colnames(env)[4:ncol(env)]
view(traits)
means <- matrix(nrow = length(unique(env$line)), ncol = length(traits))
view(means)
rownames(means) <- unique(env$line)
colnames(means) <- traits
view(means)

for (i in 1:length(traits)) {
  m <- lmer(get(traits[i]) ~ line + rep + (1|block), data = env)
  means_trait <- summary(emmeans(m, "line", df = 1))
  line_names <- means_trait$line
  
    for (j in 1:nrow(means)) {
    line_mean <- means_trait$emmean[line_names %in% rownames(means)[j]]
    if (length(line_mean) == 1) {
      means[j, i] <- line_mean
    }
  }
}
?emmeans
means <- means[order(rownames(means)), ]
View(means)
write.table(means, "MR19_means.csv", sep = ";", row.names = TRUE, col.names = NA)

# Summary statistics for BMWpop only --------------------
BMW <- read.table("../../../IDs_BMWpop.txt", header = TRUE, sep = ";")$Line
MR19_BMW <- means[rownames(means) %in% BMW, ]
MR19_BMW <- as.data.frame(MR19_BMW)
view(MR19_BMW)

summary_tbl <- summarise_all(.tbl = MR19_BMW, .funs = list(min, mean, max, sd), na.rm = TRUE)
summary_tbl <- summary_tbl[order(names(summary_tbl))]
summary_tbl <- as.data.frame(matrix(unlist(summary_tbl), byrow = TRUE, ncol = 4))
rownames(summary_tbl) <- colnames(MR19_BMW)[order(colnames(MR19_BMW))]
colnames(summary_tbl) <- c("min", "mean", "max", "SD")
write.table(summary_tbl, "MR19_summary.csv", sep = ";", col.names = NA)














#################################################

#choosing exteremes disregarding allel

x<- read.csv("MR19_means.csv", header = TRUE, sep=";")
View(x)
colnames(x)[1] <-"name"
x <- x[ , c(1,10)]
x<- x[c(3:402),]
length(x$Isotope)
length(x$name)

#mean and median calculation

#min<- min(x$Isotope)
#max<- max(x$Isotope)
#mean<- mean(x$Isotope)
#sd(x$Isotope)
#median<- median(x$Isotope)
#distance.top <- mean(x$Isotope)-(top5$Isotope)
#distance.Low <- mean(x$Isotope)-(Low5$Isotope)
##head(distance.top)
#head(distance.Low)

# loading alleles file
data_Rht <- read.csv("BMWpop_Rht_alleles.csv")
View(data_Rht)
data_Rht <- data_Rht[c(1:400),]
length(data_Rht$name)

#merging files


library(tidyverse)
data_I_A<- left_join(data_Rht, x , by=c("name"="name"))
View(data_I_A)
colnames(data_I_A)[2:3]<-c("RB1","RD1")
View(data_I_A)
data_I_A_G<- data_I_A[data_I_A$RD1=="G",]
View(data_I_A_G)

data_I_A_T<- data_I_A[data_I_A$RD1=="T",]
View(data_I_A_T)
#............sorting data. choosing top and low5.........
data_ordered_I_A_G<-arrange(data_I_A_G,desc(Isotope))
View(data_ordered_I_A_G)
top5<- data_ordered_I_A_G%>% top_n(5)
View(top5)
mean(top5$ Isotope)
Low5<- data_ordered_I_A_G %>% top_n(-5)
View(Low5)
mean(Low5$ Isotope)
final10lines<- rbind(top5,Low5)
View(final10lines)
#.......................................
data_ordered_I_A_T<-arrange(data_I_A_T,desc(Isotope))
View(data_ordered_I_A_T)
top10T<- data_ordered_I_A_T%>% top_n(10)
View(top10T)
mean(top5$ Isotope)
Low10T<- data_ordered_I_A_T %>% top_n(-10)
View(Low10T)
mean(Low5$ Isotope)
final10linesT<- rbind(top10T,Low10T)
View(final10linesT)



top10<- data_ordered_I_A_G%>% top_n(10)
View(top10)
mean(top5$ Isotope)
Low10<- data_ordered_I_A_G %>% top_n(-10)
View(Low10)
mean(Low5$ Isotope)
final20lines<- rbind(top10,Low10)
View(final20lines)


#............sorting data. choosing tot and low50.........
top50<- data_ordered_I_A %>% top_n(50)
View(top50)
low50<- data_ordered_I_A %>% top_n(-50)
View(low50)
low.top.50<- rbind(top50,low50)
View(low.top.50)
write.table(low.top.50, "extreme 100lines.csv", sep = ";", row.names = TRUE, col.names = NA)
#....................choosing top and low25lines..............
top25<- data_ordered_I_A %>% top_n(25)
View(top25)
low25<- data_ordered_I_A %>% top_n(-25)
View(low25)
low.top.25<- rbind(top25,low25)
View(low.top.25)
length(low.top.25$Isotope)
write.table(low.top.25, "extreme 50lines.csv", sep = ";", row.names = TRUE, col.names = NA)



#................heatmap of 100-extreme lines..............

library(Nbclust)
library(cluster)
library(ggplot2)
library(heatmaply)
heatmap_extremelines<- read.csv("Genetic_distance_100lines_extremes.csv", header=TRUE, sep=";")
View(heatmap_extremelines)
length(heatmap_extremelines)

heatmap(values1, scale="none")
heatmaply(values1)

values1 <- heatmap_extremelines[,2:length(heatmap_extremelines)] 
View(values1)
allColNames <- colnames(heatmap_extremelines)
rownames(values1) <- allColNames[2:length(allColNames)]
View(values1)
class(values1)
values1<- as.matrix(values1)
library(heatmaply)

png("mutantheatmap.png")

heatmaply_cor(cor(values1), xlab="lines", ylab="lines", k_col=4, k_row=4, limits =c(0.5,0.7),
              fontsize_row=5, fontsize_col=5,row_text_angle = 0,column_text_angle = 90)
dev.off()


dev.off()
#..................heatmap of 50-extreme lines.........
 
heatmap_extremelines2<- read.csv("Genetic_distance_50lines_extreme.csv", header=TRUE, sep=";")
View(heatmap_extremelines2)
str(heatmap_extremelines2)
length(heatmap_extremelines2)
class(heatmap_extremelines2)

values <- heatmap_extremelines2[,2:length(heatmap_extremelines2)] 
View(values)
allColNames <- colnames(heatmap_extremelines2)
View(allColNames)
rownames(values) <- allColNames[2:length(allColNames)]
View(values)
values<- as.matrix(values)
class(values)
library(RColorBrewer)
col<- colorRampPalette(brewer.pal(10, "RdYlBu"))(256)
heatmap(values, scale="none", col=col,
        RowSideColors = rep(c("blu","pink"), each=16),
        ColSideColors = c(rep("purple",5), rep("orange",16)))
heatmap(values)
heatmaply(values)
#.................................
BMW_genetic_distance<- read.csv("BMWpop_genetic_distance.csv", header=TRUE, sep=";")
View(BMW_genetic_distance)
class(BMW_genetic_distance)
BMW_genetic_distance<- as.matrix(BMW_genetic_distance)
class(BMW_genetic_distance)
class(values)
dim(values)
df <- as.data.frame(values)
class(df)
a = colnames(df)
a
dd<- BMW_genetic_distance[a, ]
View(dd)
