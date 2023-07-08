#########################################
#SUPPLEMENTAL FIGURES
#########################################
library(MetBrewer)
library("devtools")
library(ggplot2)
library("ggrepel")
library("FactoMineR")
library("factoextra")
library(dartR)
library(adegenet)
library(vcfR)
library(SNPRelate)
library(remotes)
library(gwscaR)
library("readxl")
library(grid)
library(gridExtra)

#########################################
#SUPPLEMENTARY FIGURE 1
#########################################

wild9 <- read.table("~/R/Cannabis_Genetics_Manuscript/pi/pi_all_datasets/lw_n498_10chr_strictfilter_10k.windowed.pi", header=T)
head(wild9)
wild9$Name<-as.factor("LeafWorks")

wild8 <- read.table("~/R/Cannabis_Genetics_Manuscript/pi/pi_all_datasets/Phylos_n845_10chr_strictfilter_pi_10k.windowed.pi", header=T)
head(wild8)
wild8$Name<-as.factor("Phylos Biosciences (n=845)")

wild7 <- read.table("~/R/Cannabis_Genetics_Manuscript/pi/pi_all_datasets/Phylos_n1378_10chr_strictfilter_pi_10k.windowed.pi", header=T)
head(wild7)
wild7$Name<-as.factor("Phylos Biosciences (n=1378)")

wild6 <- read.table("~/R/Cannabis_Genetics_Manuscript/pi/pi_all_datasets/Sunrise_n25_10chr_strictfilter_pi_10k.windowed.pi", header=T)
head(wild6)
wild6$Name<-as.factor("Sunrise Genetics")

wild5 <- read.table("~/R/Cannabis_Genetics_Manuscript/pi/pi_all_datasets/Soorni_n94_10chr_strictfilter_pi_10k.windowed.pi", header=T)
head(wild5)
wild5$Name<-as.factor("Soorni")

wild4 <- read.table("~/R/Cannabis_Genetics_Manuscript/pi/pi_all_datasets/Colorado_n162_10chr_strictfilter_pi_10k.windowed.pi", header=T)
wild4$Name<-as.factor("Colorado")

wild3 <- read.table("~/R/Cannabis_Genetics_Manuscript/pi/pi_all_datasets/Courtagen_n58_10chr_strictfilter_pi_10k.windowed.pi", header=T)
head(wild3)
wild3$Name<-as.factor("Courtagen")

wild2 <- read.table("~/R/Cannabis_Genetics_Manuscript/pi/pi_all_datasets/Kannapedia1_n61_10chr_strictfilter_pi_10k.windowed.pi", header=T)
wild2$Name<-as.factor("Medicinal Genomics (n=61)")

wild1 <- read.table("~/R/Cannabis_Genetics_Manuscript/pi/pi_all_datasets/MG_V1_n289_10chr_strictfilter_pi_10k.windowed.pi", header=T)
wild1$Name<-as.factor("Medicinal Genomics StrainSEEK V1 (n=289)")


wild_all<-rbind(wild1,wild2,wild3,wild4,wild5,wild6,wild7,wild8,wild9)
head(wild_all)
str(wild_all)
write.csv(wild_all, "~/R/Cannabis_Genetics_Manuscript/pi/pi_all_datasets/all_pi_9_datasets.csv")

wild_all<-read.csv("~/R/Cannabis_Genetics_Manuscript/pi/pi_all_datasets/all_pi_9_datasets.csv")
head(wild_all)

wild_all$CHROM<-as.factor(wild_all$CHROM)
wild_all$Name<-as.factor(wild_all$Name)
wild_all$CHROM<- factor(wild_all$CHROM,levels=unique(wild_all$CHROM))
wild_all$BIN_START<-as.numeric(wild_all$BIN_START)
str(wild_all)

#Supplementary Figure 1B
#Plot all groups by Chromosome side by side 
p1<- ggplot(data= wild_all, aes(x=BIN_START, y = PI, colour = Name, group = Name)) +
  geom_point(size=0.5) +
  #geom_smooth(se = T, method = loess) +
  facet_wrap(~CHROM) +
  theme_bw() + theme(legend.text=element_text(size=20)) +
  scale_color_manual(values=met.brewer("Hokusai1", 9))
p1
p1 + theme(axis.text = element_text(size = 8))  + theme(text = element_text(size = 20))


met.brewer("Hokusai3", 7)
met.brewer("Hokusai1", 7)
met.brewer("Hokusai1", 10)


#SUPPLEMENTARY FIGURE 1A
#Box Plot of Chromosomes
p2 <- ggplot(wild_all, aes(x = Name, y = PI, fill = CHROM)) +
  geom_boxplot(outlier.size = 0) +
  theme_bw() + theme(legend.text=element_text(size=20)) +
  scale_fill_manual(values=met.brewer("Hokusai1", 10))
p2




#########################################
#SUPPLEMENTARY FIGURE 2
#########################################
#Phylos n1378

#SNPRelate 
setwd("~/R/Cannabis_Genetics_Manuscript/Filtered_VCF_Files/")
vcf.fn <- "Phylos_n1378_10chr_strictfilter.vcf"

#Imports variants & convert to Genlight object
snpgdsVCF2GDS(vcf.fn, "phylos1378.gds", method="copy.num.of.ref")
snpgdsSummary("phylos1378.gds")
genofile <- snpgdsOpen("phylos1378.gds")
set.seed(1000)
snpset <- snpgdsLDpruning(genofile, ld.threshold=0.2,autosome.only = F) 
names(snpset)
snpset.id <- unlist(snpset)
pca <- snpgdsPCA(genofile, snp.id=snpset.id, num.thread=2, autosome.only = F)

#For PCA 
pc.percent <- pca$varprop*100
head(round(pc.percent, 2))

#make a data frame
tab <- data.frame(sample.id = pca$sample.id,
                  EV1 = pca$eigenvect[,1],
                  EV2 = pca$eigenvect[,2],
                  EV3 = pca$eigenvect[,3],
                  EV4 = pca$eigenvect[,4],
                  stringsAsFactors = FALSE)
head(tab)
write.csv(tab, "~/R/Cannabis_Genetics_Manuscript/csv_files/phylos_1378_PRJNA510566.csv")

#########################################
#Plot 2D PCA 
pca <- read.csv("~/R/Cannabis_Genetics_Manuscript/csv_files/phylos_1378_PRJNA510566_typeknown.csv")
head(pca)
p<-ggplot(pca,aes(x=EV1,y=EV2))
p<-p+geom_point()
p

pca <- read.csv("~/R/Cannabis_Genetics_Manuscript/csv_files/phylos_1378_PRJNA510566_typeknown.csv")
p<-ggplot(data=pca, aes(x=EV1, y=EV2,group=Type,shape=Type,color=Type))+
  geom_point() +
  scale_shape_manual(values=seq(0,15)) +
  scale_color_manual(values = c("turquoise","blue","red", "green", "grey","orange")) +
  stat_ellipse() +
  theme_classic() +
  xlab("PC 1 (8.14%)") + 
  ylab("PC 2 (4.61%)")
p
p + theme(axis.text = element_text(size = 10))  + theme(text = element_text(size = 30)) 


#########################################
#Plot Hierarchical Cluster
pca <- read.csv("~/R/Cannabis_Genetics_Manuscript/csv_files/phylos_1378_PRJNA510566_typeknown.csv")
head(pca)
pca2 <- pca[3:6]
head(pca2)
row.names(pca2) <- pca$Cultivar_full


res.pca <- PCA(pca2, graph = FALSE)
res.hcpc <- HCPC(res.pca, graph = FALSE)


f<- fviz_dend(res.hcpc, 
              cex = 0.03,                     # Label size
              palette = "jco",               # Color palette see ?ggpubr::ggpar
              labels_track_height = 1.0,     # Augment the room for labels
              labels_cols = TRUE, pca$Cultivar_full)
f + scale_color_manual(values=met.brewer("Hokusai1", 8))


#########################################
#SUPPLEMENTARY FIGURE 3
#########################################
#Sunrise 

#SNPRelate
setwd("/Users/annamccormick/R/VCF_analysis/nuclear3/")

setwd("~/R/Cannabis_Genetics_Manuscript/Filtered_VCF_Files/")
vcf.fn <- "Sunrise_n25_10chr_strictfilter.vcf"

#Imports variants & convert to Genlight object
snpgdsVCF2GDS(vcf.fn, "sunrise_n25.gds", method="copy.num.of.ref")
snpgdsSummary("sunrise_n25.gds")
genofile <- snpgdsOpen("sunrise_n25.gds")
set.seed(1000)
snpset <- snpgdsLDpruning(genofile, ld.threshold=0.2,autosome.only = F) 
names(snpset)
snpset.id <- unlist(snpset)
pca <- snpgdsPCA(genofile, snp.id=snpset.id, num.thread=2, autosome.only = F)

#For PCA 
pc.percent <- pca$varprop*100
head(round(pc.percent, 2))

#make a data frame
tab <- data.frame(sample.id = pca$sample.id,
                  EV1 = pca$eigenvect[,1],
                  EV2 = pca$eigenvect[,2],
                  EV3 = pca$eigenvect[,3],
                  EV4 = pca$eigenvect[,4],
                  stringsAsFactors = FALSE)
head(tab)
write.csv(tab, "~/R/Cannabis_Genetics_Manuscript/csv_files/sunrise_n25.csv")

#########################################
#Plot 2D PCA 
library("ggrepel")
pca <- read.csv("~/R/Cannabis_Genetics_Manuscript/csv_files/sunrise_n25.csv")
head(pca)
p<-ggplot(pca,aes(x=EV1,y=EV2))
p<-p+geom_point()
p

pca <- read.csv("~/R/Cannabis_Genetics_Manuscript/csv_files/sunrise_n25.csv")
p<-ggplot(data=pca, aes(x=EV1, y=EV2,group=Type,shape=Type,color=Type))+
  geom_point() +
  scale_shape_manual(values=seq(0,15)) +
  scale_color_manual(values = c("red", "grey")) +
  stat_ellipse() +
  theme_classic() +
  xlab("PC 1 (15.49%)") + 
  ylab("PC 2 (8.49%)")
p
p + theme(axis.text = element_text(size = 10))  + theme(text = element_text(size = 30)) 


#########################################
#Plot Hierarchical Cluster
library("FactoMineR")
library("factoextra")
pca <- read.csv("~/R/Cannabis_Genetics_Manuscript/csv_files/sunrise_n25.csv")
head(pca)
#pca2 <- pca[3:4]
pca2 <- pca[3:6]
head(pca2)
row.names(pca2) <- pca$HCPC
res.pca <- PCA(pca2, graph = FALSE)
res.hcpc <- HCPC(res.pca, graph = FALSE)


f<- fviz_dend(res.hcpc, 
              cex = 0.8,                     # Label size
              palette = "jco",               # Color palette see ?ggpubr::ggpar
              labels_track_height = 1.0,     # Augment the room for labels
              labels_cols = TRUE, pca$HCPC)
f + scale_color_manual(values=met.brewer("Hokusai1", 8))




#########################################
#SUPPLEMENTARY FIGURE 4
#########################################
#Colorado dataset

#SNPRelate
setwd("~/R/Cannabis_Genetics_Manuscript/Filtered_VCF_Files/")
vcf.fn <- "Colorado_n162_10chr_strictfilter.vcf"

#Imports variants & convert to Genlight object
snpgdsVCF2GDS(vcf.fn, "colorado_n162.gds", method="copy.num.of.ref")
snpgdsSummary("colorado_n162.gds")
genofile <- snpgdsOpen("colorado_n162.gds")
set.seed(1000)
snpset <- snpgdsLDpruning(genofile, ld.threshold=0.2,autosome.only = F) 
names(snpset)
snpset.id <- unlist(snpset)
pca <- snpgdsPCA(genofile, snp.id=snpset.id, num.thread=2, autosome.only = F)

#For PCA 
pc.percent <- pca$varprop*100
head(round(pc.percent, 2))

#make a data frame
tab <- data.frame(sample.id = pca$sample.id,
                  EV1 = pca$eigenvect[,1],
                  EV2 = pca$eigenvect[,2],
                  EV3 = pca$eigenvect[,3],
                  EV4 = pca$eigenvect[,4],
                  stringsAsFactors = FALSE)
head(tab)
write.csv(tab, "~/R/Cannabis_Genetics_Manuscript/csv_files/colorado_n162.csv")

#########################################
#Plot 2D PCA 
library("ggrepel")
pca <- read.csv("~/R/Cannabis_Genetics_Manuscript/csv_files/colorado_n162.csv")
head(pca)
p<-ggplot(pca,aes(x=EV1,y=EV2))
p<-p+geom_point()
p

pca <- read.csv("~/R/Cannabis_Genetics_Manuscript/csv_files/colorado_n162.csv")
p<-ggplot(data=pca, aes(x=EV1, y=EV2,group=Type,shape=Type,color=Type))+
  geom_point() +
  scale_shape_manual(values=seq(0,15)) +
  scale_color_manual(values = c("turquoise", "blue","red","orange", "green","grey")) +
  stat_ellipse() +
  theme_classic() +
  xlab("PC 1 (5.4%)") + 
  ylab("PC 2 (3.64%)")
p
p + theme(axis.text = element_text(size = 10))  + theme(text = element_text(size = 30)) 


#########################################
#Plot Hierarchical Cluster
pca <- read.csv("~/R/Cannabis_Genetics_Manuscript/csv_files/colorado_n162.csv")
head(pca)
pca2 <- pca[3:6]
head(pca2)
row.names(pca2) <- pca$HCPC
res.pca <- PCA(pca2, graph = FALSE)
res.hcpc <- HCPC(res.pca, graph = FALSE)


f<- fviz_dend(res.hcpc, 
              cex = 0.06,                     # Label size
              palette = "jco",               # Color palette see ?ggpubr::ggpar
              labels_track_height = 1.0,     # Augment the room for labels
              labels_cols = TRUE, pca$HCPC)
f + scale_color_manual(values=met.brewer("Hokusai1", 8))




#########################################
#SUPPLEMENTARY FIGURE 5
#########################################
#Courtagen

#SNPRelate
setwd("~/R/Cannabis_Genetics_Manuscript/Filtered_VCF_Files/")
vcf.fn <- "Courtagen_n58_10chr_strictfilter.vcf"

#Imports variants & convert to Genlight object
snpgdsVCF2GDS(vcf.fn, "courtagen_n58.gds", method="copy.num.of.ref")
snpgdsSummary("courtagen_n58.gds")
genofile <- snpgdsOpen("courtagen_n58.gds")
set.seed(1000)
snpset <- snpgdsLDpruning(genofile, ld.threshold=0.2,autosome.only = F) #ERROR: SNP on non-Autosomes***
names(snpset)
snpset.id <- unlist(snpset)
pca <- snpgdsPCA(genofile, snp.id=snpset.id, num.thread=2, autosome.only = F)

#For PCA 
pc.percent <- pca$varprop*100
head(round(pc.percent, 2))

#make a data frame
tab <- data.frame(sample.id = pca$sample.id,
                  EV1 = pca$eigenvect[,1],
                  EV2 = pca$eigenvect[,2],
                  EV3 = pca$eigenvect[,3],
                  EV4 = pca$eigenvect[,4],
                  stringsAsFactors = FALSE)
head(tab)
write.csv(tab, "~/R/Cannabis_Genetics_Manuscript/csv_files/courtagen_n58.csv")

#########################################
#Plot 2D PCA 
library("ggrepel")
pca <- read.csv("~/R/Cannabis_Genetics_Manuscript/csv_files/courtagen_n58.csv")
head(pca)
p<-ggplot(pca,aes(x=EV1,y=EV2))
p<-p+geom_point()
p

pca <- read.csv("~/R/Cannabis_Genetics_Manuscript/csv_files/courtagen_n58_typeknown.csv")
p<-ggplot(data=pca, aes(x=EV1, y=EV2,group=Type,shape=Type,color=Type))+
  geom_point() +
  scale_shape_manual(values=seq(0,15)) +
  scale_color_manual(values = c("turquoise", "red","green","grey")) +
  stat_ellipse() +
  theme_classic() +
  xlab("PC 1 (6.16%)") + 
  ylab("PC 2 (4.83%)")
p
p + theme(axis.text = element_text(size = 10))  + theme(text = element_text(size = 30)) 


#########################################
#Plot Hierarchical Cluster
pca <- read.csv("~/R/Cannabis_Genetics_Manuscript/csv_files/courtagen_n58_typeknown.csv")
head(pca)
pca2 <- pca[3:6]
head(pca2)
row.names(pca2) <- pca$HCPC
res.pca <- PCA(pca2, graph = FALSE)
res.hcpc <- HCPC(res.pca, graph = FALSE)


f<- fviz_dend(res.hcpc, 
              cex = 0.06,                     # Label size
              palette = "jco",               # Color palette see ?ggpubr::ggpar
              labels_track_height = 1.0,     # Augment the room for labels
              labels_cols = TRUE, pca$HCPC)
f + scale_color_manual(values=met.brewer("Hokusai1", 8))






#########################################
#SUPPLEMENTARY FIGURE 6
#########################################
#Kannapedia n61

#SNPRelate 
setwd("/Users/annamccormick/R/VCF_analysis/nuclear3/")

setwd("~/R/Cannabis_Genetics_Manuscript/Filtered_VCF_Files/")
vcf.fn <- "Kannapedia1_n61_10ch_strictfilter.vcf"

#Imports variants & convert to Genlight object
snpgdsVCF2GDS(vcf.fn, "kannapedia_n61.gds", method="copy.num.of.ref")
snpgdsSummary("kannapedia_n61.gds")
genofile <- snpgdsOpen("kannapedia_n61.gds")
set.seed(1000)
snpset <- snpgdsLDpruning(genofile, ld.threshold=0.2,autosome.only = F) #ERROR: SNP on non-Autosomes***
names(snpset)
snpset.id <- unlist(snpset)
pca <- snpgdsPCA(genofile, snp.id=snpset.id, num.thread=2, autosome.only = F)

#For PCA 
pc.percent <- pca$varprop*100
head(round(pc.percent, 2))

#make a data frame
tab <- data.frame(sample.id = pca$sample.id,
                  EV1 = pca$eigenvect[,1],
                  EV2 = pca$eigenvect[,2],
                  EV3 = pca$eigenvect[,3],
                  EV4 = pca$eigenvect[,4],
                  stringsAsFactors = FALSE)
head(tab)
write.csv(tab, "~/R/Cannabis_Genetics_Manuscript/csv_files/kannapedia_n61.csv")

#########################################
#Plot 2D PCA 
library("ggrepel")
pca <- read.csv("~/R/Cannabis_Genetics_Manuscript/csv_files/kannapedia_n61_typeknown.csv")
head(pca)
p<-ggplot(pca,aes(x=EV1,y=EV2))
p<-p+geom_point()
p

pca <- read.csv("~/R/Cannabis_Genetics_Manuscript/csv_files/kannapedia_n61_typeknown.csv")
p<-ggplot(data=pca, aes(x=EV1, y=EV2,group=Type,shape=Type,color=Type))+
  geom_point() +
  scale_shape_manual(values=seq(0,15)) +
  scale_color_manual(values = c("red","green", "grey")) +
  stat_ellipse() +
  theme_classic() +
  xlab("PC 1 (4.95%)") + 
  ylab("PC 2 (3.89%)")
p
p + theme(axis.text = element_text(size = 10))  + theme(text = element_text(size = 30)) 


#########################################
#Plot Hierarchical Cluster
pca <- read.csv("~/R/Cannabis_Genetics_Manuscript/csv_files/kannapedia_n61_typeknown.csv")
head(pca)
pca2 <- pca[3:6]
head(pca2)
row.names(pca2) <- pca$HCPC
res.pca <- PCA(pca2, graph = FALSE)
res.hcpc <- HCPC(res.pca, graph = FALSE)


f<- fviz_dend(res.hcpc, 
              cex = 0.06,                     # Label size
              palette = "jco",               # Color palette see ?ggpubr::ggpar
              labels_track_height = 1.0,     # Augment the room for labels
              labels_cols = TRUE, pca$HCPC)
f + scale_color_manual(values=met.brewer("Hokusai1", 8))











#########################
#Core collection identification from datasets
#########################
library(vcfR)
library(gdsfmt)
library(SNPRelate)
library(gwscaR)
library(gt)
library(snpsel)
library(ggplot2)
library(poppr)
library(reshape2)
library(ggplot2)
library(dartR)
library(adegenet)
library(poppr)
library(ape)
library(RColorBrewer)
library(adegenet)
library(rCNV)

#read in vcf
setwd("/Users/annamccormick/R/VCF_analysis/nuclear3/core_collection")

#Choose VCF
#vcf <- read.vcfR("Phylos_n845_10chr_strictfilter.vcf")
#vcf <- read.vcfR("Phylos_n1378_10chr_strictfilter_nometa.vcf")
#vcf <- read.vcfR("lw_n498_10chr_strictfilter.vcf")
#vcf <- read.vcfR("Soorni_n94_10chr_strictfilter.vcf")
#vcf <- read.vcfR("Sunrise_n25_10chr_strictfilter.vcf")
#vcf <- read.vcfR("Courtagen_n58_10chr_strictfilter.vcf")
#vcf <- read.vcfR("Colorado_n162_10chr_strictfilter.vcf")
#vcf <- read.vcfR("Kannapedia1_n61_10ch_strictfilter.vcf")
vcf <- read.vcfR("MG_V1_n289_10ch_strictfilter.vcf")

####convert genlight object
x <- vcfR2genlight(vcf)

#create distance matrix
x.dist <- dist(x)
ploidy(x) <- 2
#tree <- aboot(gl.rubi, tree = "upgma", distance = bitwise.dist, sample = 100, showtree = F, cutoff = 50, quiet = T)

#change to genid
x<- vcfR2genind(vcf)
#set ploidy
ploidy(x) <- 2
xmlg<-mlg(x, quiet = FALSE)

#genetic distance
x.dist <- poppr::bitwise.dist(x)
hist(x.dist)
heatmap(as.matrix(x.dist)) 

#write.csv(x.dist,"phylos_n845_distance_matrix.csv")

# Ward Hierarchical Clustering
fit <- hclust(x.dist, method="ward") 
plot(fit) # display dendogram
groups <- cutree(fit, k=11) # cut tree into 5 clusters
# draw dendogram with red borders around the 5 clusters 
rect.hclust(fit, k=11, border="red")


#write.csv(as.matrix(x.dist), "distance_sug.csv")
thresholds<-mlg.filter(x, threshold = 0.05, distance = "bitwise.dist", missing="mean",threads = 1L)
cutoff_predictor(thresholds, fraction = 0.05)

pthresh  <- filter_stats(x, distance = "bitwise.dist",
                         plot = TRUE, stats = "THRESHOLD", threads = 1L)
cutoff_predictor(pthresh$farthest)

# prediction for all algorithms
sapply(pthresh, cutoff_predictor) #take note of the farthest number for later use
tab<-mlg.table(x)
diversity_boot(tab, 10L)
diversity_ci(tab, 10L, rarefy = F, raw = FALSE)
diversity_stats(tab)

#The threshold used to define MLGs was predicted using cutoff_predictor with the “farthest” algorithm. 
cc_test <- clonecorrect(x)
cc_test
nInd(x) # 1903
# How many after we clone corrected for country, year, and month?
nInd(cc_test) # 1190

#write.csv(x$tab, "test_1488.csv")

sbs<-as.genclone(x)
mll(sbs) 
mlg.filter(sbs, threshold = 0.05, distance = "bitwise.dist", missing="mean",threads = 1L)
raw_dist <- function(sbs){
  dist(genind2df(sbs, usepop = FALSE))
}
(xdis <- raw_dist(sbs))


#mlg.filter(sbs)<- 0.04350649 #phylos_n845
#mlg.filter(sbs)<- 0.0461432507 #phylos_n1378
#mlg.filter(sbs)<- 0.093060498  #lw
#mlg.filter(sbs)<- 0.1162538  #soorni
#mlg.filter(sbs)<- 0.08500553  #sunrise
#mlg.filter(sbs)<- 0.040192926  #courtagen
#mlg.filter(sbs)<- 0.02050342  #colorado
#mlg.filter(sbs)<- 0.03404658  #kan_61
mlg.filter(sbs)<- 0.01073374  #MG_V1


dups_sug = mlg.id(sbs)

for (i in dups_sug){ # for each element in the list object
  if (length(dups_sug[i]) > 1){ # if the length is greater than 1
    print(i) # print individuals that are duplicates
  }
}


library(corehunter)
# precomputed distance matrix
my.data <- distances(as.matrix(x.dist))
core <- sampleCore(my.data, size = 25)
#core <- sampleCore(my.data, size = 10)
core

