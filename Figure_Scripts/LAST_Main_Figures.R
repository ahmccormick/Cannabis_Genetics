#########################################
#MAIN FIGURES
#########################################
library(dartR)
library(adegenet)
library(vcfR)
library(SNPRelate)
library(remotes)
library(gwscaR)
library("readxl")
library(ggplot2)
library(grid)
library(gridExtra)
library("FactoMineR")
library("factoextra")
library(MetBrewer)

#########################################
#FIGURE 1A
#########################################
#LeafWorks

#SNPRelate
setwd("~/R/Cannabis_Genetics_Manuscript/Filtered_VCF_Files/")
vcf.fn <- "lw_n498_10chr_strictfilter.vcf"

#Imports variants & convert to Genlight object
snpgdsVCF2GDS(vcf.fn, "lw498.gds", method="copy.num.of.ref")
snpgdsSummary("lw498.gds")
genofile <- snpgdsOpen("lw498.gds")
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
write.csv(tab, "~/R/Cannabis_Genetics_Manuscript/csv_files/lw_n498.csv")

#########################################
#Plot 2D PCA 
library("ggrepel")
pca <- read.csv("~/R/Cannabis_Genetics_Manuscript/csv_files/lw_n498_typeknown.csv")
head(pca)
p<-ggplot(pca,aes(x=EV1,y=EV2))
p<-p+geom_point()
p

#Plot by use-type 
p<-ggplot(data=pca, aes(x=EV1, y=EV2,group=Type,shape=Type,color=Type))+
  geom_point() +
  scale_shape_manual(values=seq(0,15)) +
  scale_color_manual(values = c("turquoise", "blue","red","orange", "green", "grey")) +
  stat_ellipse() +
  theme_classic() +
  xlab("PC 1 (5.52%)") + 
  ylab("PC 2 (3.57%)")
p
p + theme(axis.text = element_text(size = 10))  + theme(text = element_text(size = 30)) 

#########################################
#Plot Hierarchical Cluster
library("FactoMineR")
library("factoextra")
pca <- read.csv("~/R/Cannabis_Genetics_Manuscript/csv_files/lw_n498_typeknown.csv")
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
              labels_cols = TRUE, pca$Type)
f + scale_color_manual(values=met.brewer("Hokusai1", 8))

#########################################
#FIGURE 1C
#########################################
#Phylos 845

#SNPRelate 
setwd("~/R/Cannabis_Genetics_Manuscript/Filtered_VCF_Files/")
vcf.fn <- "Phylos_n845_10chr_strictfilter.vcf"

#Imports variants & convert to Genlight object
snpgdsVCF2GDS(vcf.fn, "phylos845.gds", method="copy.num.of.ref")
snpgdsSummary("phylos845.gds")
genofile <- snpgdsOpen("phylos845.gds")
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
write.csv(tab, "~/R/Cannabis_Genetics_Manuscript/csv_files/phylos_845_PRJNA347566.csv")

#########################################
#Plot 2D PCA
library("ggrepel")
pca <- read.csv("~/R/Cannabis_Genetics_Manuscript/csv_files/phylos_845_PRJNA347566.csv")
head(pca)
p<-ggplot(pca,aes(x=EV1,y=EV2))
p<-p+geom_point()
p


pca <- read.csv("~/R/Cannabis_Genetics_Manuscript/csv_files/phylos_845_PRJNA347566_typeknown_pca.csv")
p<-ggplot(data=pca, aes(x=EV1, y=EV2,group=Type_final,shape=Type_final,color=Type_final))+
  geom_point() +
  scale_shape_manual(values=seq(0,15)) +
  scale_color_manual(values = c("grey","turquoise", "blue" , "red", "green", "grey2")) +
  stat_ellipse() +
  theme_classic() +
  xlab("PC 1 (6.66%)") + 
  ylab("PC 2 (5.34%)")
p
p + theme(axis.text = element_text(size = 10))  + theme(text = element_text(size = 30)) 



#########################################
#Plot Hierarchical Cluster
library("FactoMineR")
library("factoextra")
pca <- read.csv("~/R/Cannabis_Genetics_Manuscript/csv_files/phylos_845_PRJNA347566_typeknown_hcpc.csv")
head(pca)
#pca2 <- pca[3:4]
pca2 <- pca[3:6]
head(pca2)
row.names(pca2) <- pca$Cultivar_new
res.pca <- PCA(pca2, graph = FALSE)
res.hcpc <- HCPC(res.pca, graph = FALSE)

library(MetBrewer)
f<- fviz_dend(res.hcpc, 
              cex = 0.06,                     # Label size
              palette = "jco",               # Color palette see ?ggpubr::ggpar
              labels_track_height = 1.0,     # Augment the room for labels
              labels_cols = TRUE, pca$Type)
f + scale_color_manual(values=met.brewer("Hokusai1", 8))




#########################################
#FIGURE 2A
#########################################
#Soorni et al 

#SNPRelate 
setwd("~/R/Cannabis_Genetics_Manuscript/Filtered_VCF_Files/")
vcf.fn <- "Soorni_n94_10chr_strictfilter.vcf"

#Imports variants & convert to Genlight object
snpgdsVCF2GDS(vcf.fn, "soorni_n94.gds", method="copy.num.of.ref")
snpgdsSummary("soorni_n94.gds")
genofile <- snpgdsOpen("soorni_n94.gds")
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
write.csv(tab, "~/R/Cannabis_Genetics_Manuscript/csv_files/soorni_n94.csv")

#########################################
#Plot 2D PCA
library("ggrepel")
pca <- read.csv("~/R/Cannabis_Genetics_Manuscript/csv_files/soorni_n94.csv")
head(pca)
p<-ggplot(pca,aes(x=EV1,y=EV2))
p<-p+geom_point()
p

pca <- read.csv("~/R/Cannabis_Genetics_Manuscript/csv_files/soorni_n94_typeknown_pca.csv")
p<-ggplot(data=pca, aes(x=EV1, y=EV2,group=Type,shape=Type,color=Type))+
  geom_point() +
  scale_shape_manual(values=seq(0,15)) +
  scale_color_manual(values = c("red", "orange","green")) +
  stat_ellipse() +
  theme_classic() +
  xlab("PC 1 (4.84%)") + 
  ylab("PC 2 (2.85%)")
p
p + theme(axis.text = element_text(size = 10))  + theme(text = element_text(size = 30)) 


#########################################
#Plot Hierarchical Cluster
library("FactoMineR")
library("factoextra")
library(MetBrewer)

pca <- read.csv("~/R/Cannabis_Genetics_Manuscript/csv_files/soorni_n94_typeknown_hcpc.csv")
head(pca)
#pca2 <- pca[3:4]
pca2 <- pca[3:6]
head(pca2)
row.names(pca2) <- pca$HCPC
res.pca <- PCA(pca2, graph = FALSE)
res.hcpc <- HCPC(res.pca, graph = FALSE)


f<- fviz_dend(res.hcpc, 
              cex = 0.50,                     # Label size
              palette = "jco",               # Color palette see ?ggpubr::ggpar
              labels_track_height = 1.0,     # Augment the room for labels
              labels_cols = TRUE, pca$Type)
f + scale_color_manual(values=met.brewer("Hokusai1", 8))



#########################################
#FIGURE 2C
#########################################
###MEDICINAL GENOMICS STRAINSEEK V1

#SNPRelate 
setwd("~/R/Cannabis_Genetics_Manuscript/Filtered_VCF_Files/")
vcf.fn <- "MG_V1_n289_10ch_strictfilter.vcf"

#Imports variants & convert to Genlight object
snpgdsVCF2GDS(vcf.fn, "mgv1_n289_2.gds", method="copy.num.of.ref")
snpgdsSummary("mgv1_n289_2.gds")
genofile <- snpgdsOpen("mgv1_n289_2.gds")
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
write.csv(tab, "~/R/Cannabis_Genetics_Manuscript/csv_files/MG_V1_n289.csv")

#########################################
#Plot 2D PCA
library("ggrepel")
pca <- read.csv("~/R/Cannabis_Genetics_Manuscript/csv_files/MG_V1_n289.csv")
head(pca)
p<-ggplot(pca,aes(x=EV1,y=EV2))
p<-p+geom_point()
p

pca <- read.csv("~/R/Cannabis_Genetics_Manuscript/csv_files/MG_V1_n289_typeknown_pca.csv")
p<-ggplot(data=pca, aes(x=EV1, y=EV2,group=Type,shape=Type,color=Type))+
  geom_point() +
  scale_shape_manual(values=seq(0,15)) +
  scale_color_manual(values = c("red", "orange","green")) +
  stat_ellipse() +
  theme_classic() +
  xlab("PC 1 (4.84%)") + 
  ylab("PC 2 (3.44%)")
p
p + theme(axis.text = element_text(size = 10))  + theme(text = element_text(size = 30)) 


#########################################
#Plot Hierarchical Cluster
library("FactoMineR")
library("factoextra")
pca <- read.csv("~/R/Cannabis_Genetics_Manuscript/csv_files/MG_V1_n289_typeknown_hcpc.csv")
head(pca)
pca2 <- pca[3:6]
head(pca2)
row.names(pca2) <- pca$Cultivar_Type
res.pca <- PCA(pca2, graph = FALSE)
res.hcpc <- HCPC(res.pca, graph = FALSE)


f<- fviz_dend(res.hcpc, 
              cex = 0.20,                     # Label size
              palette = "jco",               # Color palette see ?ggpubr::ggpar
              labels_track_height = 1.0,     # Augment the room for labels
              labels_cols = TRUE, pca$Type)
f + scale_color_manual(values=met.brewer("Hokusai1", 8))



#########################################
#FIGURE 3
#########################################
#Plot by use-type 

#LW 498
library("ggrepel")
pca <- read.csv("~/R/Cannabis_Genetics_Manuscript/csv_files/lw_n498_typeknown.csv")
p<-ggplot(data=pca, aes(x=EV1, y=EV2,group=Type,shape=Type,color=Type))+
  geom_point() +
  scale_shape_manual(values=seq(0,15)) +
  scale_color_manual(values = c("turquoise", "blue","red","orange", "green", "grey")) +
  stat_ellipse() +
  theme_classic() +
  xlab("PC 1 (5.52%)") + 
  ylab("PC 2 (3.57%)")
p
p + theme(axis.text = element_text(size = 10))  + theme(text = element_text(size = 30)) 



#Phylos 845
pca <- read.csv("~/R/Cannabis_Genetics_Manuscript/csv_files/phylos_845_PRJNA347566_typeknown_pca.csv")
p<-ggplot(data=pca, aes(x=EV1, y=EV2,group=Type_final,shape=Type_final,color=Type_final))+
  geom_point() +
  scale_shape_manual(values=seq(0,15)) +
  scale_color_manual(values = c("grey","turquoise", "blue" , "red", "green", "grey2")) +
  stat_ellipse() +
  theme_classic() +
  xlab("PC 1 (6.66%)") + 
  ylab("PC 2 (5.34%)")
p
p + theme(axis.text = element_text(size = 10))  + theme(text = element_text(size = 30)) 

#Soorni 
pca <- read.csv("~/R/Cannabis_Genetics_Manuscript/csv_files/soorni_n94_typeknown_pca.csv")
p<-ggplot(data=pca, aes(x=EV1, y=EV2,group=Type,shape=Type,color=Type))+
  geom_point() +
  scale_shape_manual(values=seq(0,15)) +
  scale_color_manual(values = c("red", "orange","green")) +
  stat_ellipse() +
  theme_classic() +
  xlab("PC 1 (4.84%)") + 
  ylab("PC 2 (2.85%)")
p
p + theme(axis.text = element_text(size = 10))  + theme(text = element_text(size = 30)) 


#MGV1
pca <- read.csv("~/R/Cannabis_Genetics_Manuscript/csv_files/MG_V1_n289_typeknown_pca.csv")
p<-ggplot(data=pca, aes(x=EV1, y=EV2,group=Type,shape=Type,color=Type))+
  geom_point() +
  scale_shape_manual(values=seq(0,15)) +
  scale_color_manual(values = c("red", "orange","green")) +
  stat_ellipse() +
  theme_classic() +
  xlab("PC 1 (4.84%)") + 
  ylab("PC 2 (3.44%)")
p
p + theme(axis.text = element_text(size = 10))  + theme(text = element_text(size = 30)) 




#########################################
#FIGURE 4
#########################################
#LANRDACE AND DOMESTICATES LINES

#Nucleotide plots
library(MetBrewer)

#Merging 10k pi
wild1 <- read.table("~/R/Cannabis_Genetics_Manuscript/pi/Phylos_landrace_n127_10chr_strictfilter_10k.windowed.pi", header=T)
head(wild1)
wild1$Name<-as.factor("Landrace")

wild2 <- read.table("~/R/Cannabis_Genetics_Manuscript/pi/Phylos_domesticate_n718_10chr_strictfilter_10k.windowed.pi", header=T)
wild2$Name<-as.factor("Domesticates")

wild_all<-rbind(wild1,wild2)
head(wild_all)
str(wild_all)
write.csv(wild_all, "~/R/Cannabis_Genetics_Manuscript/pi/both_phylos_pi.csv")


#Phylos
wild_all<-read.csv("~/R/Cannabis_Genetics_Manuscript/pi/both_phylos_pi.csv")

#Plot all the datasets for all chroms side by side for diversity
p1 <- ggplot(wild_all, aes(x = Name, y = PI, fill = CHROM)) +
  geom_boxplot(outlier.size = 0) +
  theme_bw() + theme(legend.text=element_text(size=20)) +
  scale_fill_manual(values=met.brewer("Hokusai1", 10)) +
  theme(axis.text = element_text(size = 20))  + theme(text = element_text(size = 20))
p1

wild_all$CHROM<-as.factor(wild_all$CHROM)
wild_all$Name<-as.factor(wild_all$Name)
wild_all$CHROM<- factor(wild_all$CHROM,levels=unique(wild_all$CHROM))
wild_all$BIN_START<-as.numeric(wild_all$BIN_START)
str(wild_all)

###Plot for all chromosomes side by side only
p2<- ggplot(data= wild_all, aes(x=BIN_START, y = PI, colour = Name, group = Name)) +
  geom_point() +
  #geom_smooth(se = T, method = loess) +
  facet_wrap(~CHROM) +
  theme_bw() + theme(legend.text=element_text(size=20)) +
  scale_color_manual(values=met.brewer("Hokusai1", 2)) +
  theme(axis.text = element_text(size = 10))  + theme(text = element_text(size = 20))
p2


#LeafWorks 
wild_all<-read.csv("~/R/Cannabis_Genetics_Manuscript/pi/both_leafworks_pi.csv")
#Plot all the datasets for all chroms side by side for diversity
p1 <- ggplot(wild_all, aes(x = Name, y = PI, fill = CHROM)) +
  geom_boxplot(outlier.size = 0) +
  theme_bw() + theme(legend.text=element_text(size=20)) +
  scale_fill_manual(values=met.brewer("Hokusai1", 10)) +
  theme(axis.text = element_text(size = 20))  + theme(text = element_text(size = 20))
p1

wild_all$CHROM<-as.factor(wild_all$CHROM)
wild_all$Name<-as.factor(wild_all$Name)
wild_all$CHROM<- factor(wild_all$CHROM,levels=unique(wild_all$CHROM))
wild_all$BIN_START<-as.numeric(wild_all$BIN_START)
str(wild_all)

###Plot for all chromosomes side by side only
p2<- ggplot(data= wild_all, aes(x=BIN_START, y = PI, colour = Name, group = Name)) +
  geom_point() +
  #geom_smooth(se = T, method = loess) +
  facet_wrap(~CHROM) +
  theme_bw() + theme(legend.text=element_text(size=20)) +
  scale_color_manual(values=met.brewer("Hokusai1", 2)) +
  theme(axis.text = element_text(size = 10))  + theme(text = element_text(size = 20))
p2


#########################################
#FIGURE 5
#########################################
#MAP for Landrace Locations

library(dismo)
library(raster)
library(maptools)
library(rasterVis)
data(wrld_simpl)
library(readr)
library(MetBrewer)

can <-read.csv("~/R/Cannabis_Genetics_Manuscript/csv_files/landrace_origins.csv",header=T, sep=",")
head(can) 
library(maptools)
data(wrld_simpl)
plot(wrld_simpl,xlim = c(50, 140), ylim = c(0, 50), axes=TRUE, col="grey") 
box()
points(can$lon, can$lat, col='black', pch=20, cex=0.75)
plot(wrld_simpl,xlim=c(50, 140),ylim=c(0, 50),col='olivedrab3',bg='lightblue', axes=TRUE)
points(can$lon, can$lat, col='black', pch=20, cex=0.75)

library(marmap)
eurasia <- getNOAA.bathy(lon1 = 50, lon2 = 140,
                         lat1 = 0, lat2 = 50, resolution = 10)
plot(eurasia)
plot(eurasia, image = TRUE, land = TRUE, axes = TRUE, lwd=0.1, bpal = list(c(0, max(eurasia), grey(.7), grey(.9), grey(.95)), c(min(eurasia), 0, "darkblue", "lightblue"))) 
points(can$lon, can$lat, col='black', pch=20, cex=1.5)
plot(eurasia, n = 1, lwd = 0.5, add = TRUE, axes=TRUE)

#FIGURE 5B
#HCPC clustering
pca <- read.csv("~/R/Cannabis_Genetics_Manuscript/csv_files/lw_nuclear_lr_subset_only.csv")
head(pca)
pca2 <- pca[3:6]
head(pca2)
row.names(pca2) <- pca$Cultivar
res.pca <- PCA(pca2, graph = FALSE)
fviz_screeplot(res.pca, addlabels = TRUE, ylim = c(0, 50))
res.hcpc <- HCPC(res.pca, graph = FALSE)


f<- fviz_dend(res.hcpc, 
              cex = 1.0,                     # Label size
              palette = "jco",               # Color palette see ?ggpubr::ggpar
              labels_track_height = 0.8,     # Augment the room for labels
              labels_cols = TRUE, pca$Type)
f + scale_color_manual(values=met.brewer("Hokusai1", 4))


#FIGURE 5C
#PCA (Regional Association) HK/LV/MB
pca <- read.csv("~/R/Cannabis_Genetics_Manuscript/csv_files/lw_nuclear_lr_subset_only.csv")
head(pca)

p<-ggplot(data=pca, aes(x=EV1, y=EV2,group=Ecotype,shape=Ecotype,color=Ecotype))+
  geom_point() + 
  theme_bw() +
  scale_shape_manual(values=seq(0,15)) +
  stat_ellipse() +
  theme(legend.text=element_text(size=20)) +
  scale_color_manual(values=met.brewer("Hokusai1", 10)) +
  theme(axis.text = element_text(size = 10))  + theme(text = element_text(size = 20)) +
  xlab("PC 1 (9.94%)") + 
  ylab("PC 2 (8.55%)")
p


#FIGURE 5D
#Plotting multiple data sets side by side
wild1 <- read.table("~/R/Cannabis_Genetics_Manuscript/pi/lw_nuclear_m_lr_hindu_strict_revcf_10k.windowed.pi", header=T)
head(wild1)
wild1$Name<-as.factor("Hindu_Kush")

wild2 <- read.table("~/R/Cannabis_Genetics_Manuscript/pi/lw_nuclear_m_lr_lolab_strict_revcf_10k.windowed.pi", header=T)
wild2$Name<-as.factor("Lolab_Valley")

wild3 <- read.table("~/R/Cannabis_Genetics_Manuscript/pi/lw_nuclear_m_lr_burmese_strict_revcf_10k.windowed.pi", header=T)
wild3$Name<-as.factor("Burma")

wild_all<-rbind(wild1,wild2,wild3)
head(wild_all)
str(wild_all)
write.csv(wild_all, "~/R/Cannabis_Genetics_Manuscript/pi/lw_landrace3_all_pi.csv")

library(ggplot2)
wild_all<-read.csv("~/R/Cannabis_Genetics_Manuscript/pi/lw_landrace3_all_pi.csv")
#Plot all the datasets for all chroms side by side for diversity
ggplot(wild_all, aes(x = Name, y = PI, fill = CHROM)) +
  geom_boxplot(outlier.size = 0) +
  #geom_point(pch = 21) +
  theme_bw()

p2<- ggplot(data= wild_all, aes(x=BIN_START, y = PI, colour = Name, group = Name)) +
  geom_point() +
  #geom_smooth(se = T, method = loess) +
  facet_wrap(~CHROM) +
  theme_bw() + theme(legend.text=element_text(size=20)) +
  scale_color_manual(values=met.brewer("Hokusai1", 4)) +
  theme(axis.text = element_text(size = 10))  + theme(text = element_text(size = 20))
p2











