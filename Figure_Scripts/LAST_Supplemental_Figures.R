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
pfin<- p1 + theme(axis.text = element_text(size = 8))  + theme(text = element_text(size = 20))

ggsave("SupB_FIN.pdf", plot = pfin, width = 20, height = 9, device = "pdf")

met.brewer("Hokusai3", 7)
met.brewer("Hokusai1", 7)
met.brewer("Hokusai1", 10)

#SUPPLEMENTARY FIGURE 1A
#Box Plot of Chromosomes
p2 <- ggplot(wild_all, aes(x = Name, y = PI, fill = CHROM)) +
  geom_boxplot(outlier.size = 0) +
  theme_bw() +
  theme(
    legend.text = element_text(size = 20),
    axis.text.x = element_text(size = 8)  # Set the size for x-axis text
  ) +
  scale_fill_manual(values = met.brewer("Hokusai1", 10))

# Display the plot
print(p2)

# Save the plot to a PDF file
ggsave("Sup1A_FIN.pdf", plot = p2, width = 20, height = 9, device = "pdf")


#########################################
#SUPPLEMENTARY FIGURE 2A Phylos n1378 PCA
#########################################

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

#####################
#Supplemental Figure 2A Phylos 1378 PCA
#####################

pca <- read.csv("~/R/Cannabis_Genetics_Manuscript/csv_files/phylos_1378_PRJNA510566_typeknown_clarke.csv")

p <- ggplot(data = pca, aes(x=EV1, y=EV2,group=Clarke_type,shape=Clarke_type,color=Clarke_type)) +
  geom_point(size=4) +
  scale_shape_manual(values = seq(0, 15)) +
  scale_color_manual(values = c("turquoise","blue","red", "green", "darkgreen","grey")) +
  stat_ellipse(size = 1.5) +  # Adjust ellipse line thickness
  theme_classic() +
  xlab("PC 1 (8.14%)") + 
  ylab("PC 2 (4.61%)") +
  labs(
    shape = "Type",  # Custom legend title for shape
    color = "Type"   # Custom legend title for color
  ) +
  theme(
    axis.text = element_text(size = 25),  # Increase axis text size
    axis.title = element_text(size = 25),  # Increase axis title size
    legend.text = element_text(size = 25),  # Increase legend text size
    legend.title = element_text(size = 25)  # Increase legend title size
  )
p
ggsave("Phylos1378_PCA.pdf", plot = p, width = 20, height = 13, device = "pdf")

#####################
#Supplemental Figure 2B Phylos 1378 Hierarchical clustering
#####################
pca <- read.csv("~/R/Cannabis_Genetics_Manuscript/csv_files/phylos_1378_PRJNA510566_typeknown_clarke.csv")

# Process the data for PCA
pca2 <- pca[3:6]
row.names(pca2) <- pca$Cultivar_full

# Perform PCA
res.pca <- PCA(pca2, graph = FALSE)

# Perform hierarchical clustering on the principal components
res.hcpc <- HCPC(res.pca, graph = FALSE)

# Read the order data for the color bar
order_data <- read.csv("~/R/Cannabis_Genetics_Manuscript/csv_files/Phylos1378_HCPC_order_plots.csv")

# Define the color palette for the dendrogram
hokusai_palette <- met.brewer("Hokusai1", 8)

# Plot the dendrogram with the specified color palette
f <- fviz_dend(res.hcpc, 
               cex = 0.06, 
               palette = "jco", 
               labels_track_height = 1.0, 
               labels_cols = TRUE) +
  scale_color_manual(values = hokusai_palette) +
  theme_minimal() +
  theme(legend.position = "none")

# Define a custom color palette for the types in the color bar
my_custom_palette <- c("Type I" = "red",
                       "Type II" = "orange",
                       "Type III" = "green",
                       "Type III / Hemp" = "darkgreen",
                       "Landrace" = "blue",
                       "Unknown" = "grey",
                       "Hemp" = "turquoise")

# Ensure that the BAR column is a factor with levels corresponding to the colors defined
order_data$BAR <- factor(order_data$BAR, levels = names(my_custom_palette))

# Create the color bar dataframe
color_bar_df <- data.frame(Type = order_data$BAR, x = 1:nrow(order_data))

# Adjust the y position to move the bar up
bar_position_y <- -0.5  # Adjust this value as needed to move the bar closer or further away

# Add the color bar with adjusted height
bar_height <- 0.2
f <- f + geom_tile(data = color_bar_df, aes(x = x, y = -0.5, fill = Type), 
                   width = 1, height = bar_height) +
  scale_fill_manual(values = my_custom_palette, 
                    name = "Type") +  # Include the legend
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.margin = margin(t = 1, b = 1, unit = "cm"),
        legend.position = "right")  # Adjust legend position if needed

# Plot the dendrogram with the color bar
print(f)
# Plot the dendrogram with the color bar
print(f)

# Save the plot to a PDF file
ggsave("Phylos1378_HCPC_FIN.pdf", plot = f, width = 20, height = 9, device = "pdf")

#####################
#Supplemental Figure 2C Phylos 1378 fastSTRUCTURE 
#####################
library(MetBrewer)
library(ggplot2)
library(MetBrewer)
library(ggplot2)

# Load the Hokusai palette
hiroshige_palette <- met.brewer("Hiroshige", 8)

selected_colors <- hiroshige_palette[c(1, 6)]
###########
# K2 Plot
###########
test <- read.csv("~/R/Cannabis_Genetics_Manuscript/fastSTRUCTURE/plots/Phylos1378_K2.csv")
long_data <- pivot_longer(test, cols = c("K1", "K2"), names_to = "category", values_to = "value")

K2 <- ggplot(long_data, aes(x = Name, y = value, fill = category, group = Type)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = selected_colors) +  # Use selected colors from Hokusai palette
  labs(title = "", x = "", y = "Percent Identity") +
  facet_grid(~ Type, scales = "free_x") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 0.5),
    axis.text.y = element_text(size = rel(0.8)),
    plot.margin = margin(5.5, 5.5, 5.5, 5.5, "pt"),
    strip.background = element_blank(),
    legend.position = "right"
  )
K2

ggsave("Phylos1378_K2_stru.pdf", plot = K2, width = 20, height = 3, device = "pdf")

###########
# K3 Plot
###########
hiroshige_palette <- met.brewer("Hiroshige", 8)
selected_colors <- hiroshige_palette[c(1, 4, 6)]

test <- read.csv("~/R/Cannabis_Genetics_Manuscript/fastSTRUCTURE/plots/Phylos1378_K3.csv")
long_data <- pivot_longer(test, cols = c("K1", "K2", "K3"), names_to = "category", values_to = "value")

K3 <- ggplot(long_data, aes(x = Name, y = value, fill = category, group = Type)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = selected_colors) +  # Use selected colors from Hokusai palette
  labs(title = "", x = "", y = "Percent Identity") +
  facet_grid(~ Type, scales = "free_x") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 0.5),
    axis.text.y = element_text(size = rel(0.8)),
    plot.margin = margin(5.5, 5.5, 5.5, 5.5, "pt"),
    strip.background = element_blank(),
    legend.position = "right"
  )
K3

ggsave("Phylos1378_K3_stru.pdf", plot = K3, width = 20, height = 3, device = "pdf")

###########
# K4 Plot
###########
hiroshige_palette <- met.brewer("Hiroshige", 8)
selected_colors <- hiroshige_palette[c(1, 4,6, 8)]

test <- read.csv("~/R/Cannabis_Genetics_Manuscript/fastSTRUCTURE/plots/Phylos1378_K4.csv")
long_data <- pivot_longer(test, cols = c("K1", "K2", "K3", "K4"), names_to = "category", values_to = "value")

K4 <- ggplot(long_data, aes(x = Name, y = value, fill = category, group = Type)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = selected_colors) +  # Use selected colors from Hokusai palette
  labs(title = "", x = "", y = "Percent Identity") +
  facet_grid(~ Type, scales = "free_x") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 0.5),
    axis.text.y = element_text(size = rel(0.8)),
    plot.margin = margin(5.5, 5.5, 5.5, 5.5, "pt"),
    strip.background = element_blank(),
    legend.position = "right"
  )
K4

ggsave("Phylos1378_K4_stru.pdf", plot = K4, width = 20, height = 3, device = "pdf")

###########
# K5 Plot
###########
hiroshige_palette <- met.brewer("Hiroshige", 8)
selected_colors <- hiroshige_palette[c(1, 3,5,7, 8)]

test <- read.csv("~/R/Cannabis_Genetics_Manuscript/fastSTRUCTURE/plots/Phylos1378_K5.csv")
long_data <- pivot_longer(test, cols = c("K1", "K2", "K3", "K4", "K5"), names_to = "category", values_to = "value")

K5 <- ggplot(long_data, aes(x = Name, y = value, fill = category, group = Type)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = selected_colors) +  # Use selected colors from Hokusai palette
  labs(title = "", x = "", y = "Percent Identity") +
  facet_grid(~ Type, scales = "free_x") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 0.5),
    axis.text.y = element_text(size = rel(0.8)),
    plot.margin = margin(5.5, 5.5, 5.5, 5.5, "pt"),
    strip.background = element_blank(),
    legend.position = "right"
  )
K5

ggsave("Phylos1378_K5_stru.pdf", plot = K5, width = 20, height = 3, device = "pdf")


#########################################
#SUPPLEMENTARY FIGURE 3A - Sunrise PCA
#########################################
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

pca <- read.csv("~/R/Cannabis_Genetics_Manuscript/csv_files/sunrise_n25.csv")
p <- ggplot(data = pca, aes(x=EV1, y=EV2,group=Type,shape=Type,color=Type)) +
  geom_point(size=4) +
  scale_shape_manual(values = seq(0, 15)) +
  scale_color_manual(values = c("red", "grey")) +
  stat_ellipse(size = 1.5) +  # Adjust ellipse line thickness
  theme_classic() +
  xlab("PC 1 (15.49%)") + 
  ylab("PC 2 (8.49%)") +
  labs(
    shape = "Type",  # Custom legend title for shape
    color = "Type"   # Custom legend title for color
  ) +
  theme(
    axis.text = element_text(size = 25),  # Increase axis text size
    axis.title = element_text(size = 25),  # Increase axis title size
    legend.text = element_text(size = 25),  # Increase legend text size
    legend.title = element_text(size = 25)  # Increase legend title size
  )
p
ggsave("Sunrise_PCA.pdf", plot = p, width = 20, height = 13, device = "pdf")

#########################################
#SUPPLEMENTARY FIGURE 3B  Sunrise - Hierarchical clustering
#########################################
pca <- read.csv("~/R/Cannabis_Genetics_Manuscript/csv_files/sunrise_n25.csv")
# Process the data for PCA
pca2 <- pca[3:6]
row.names(pca2) <- pca$HCPC

# Perform PCA
res.pca <- PCA(pca2, graph = FALSE)

# Perform hierarchical clustering on the principal components
res.hcpc <- HCPC(res.pca, graph = FALSE)

# Read the order data for the color bar
order_data <- read.csv("~/R/Cannabis_Genetics_Manuscript/csv_files/sunrise_HCPC_order_plots.csv")

# Define the color palette for the dendrogram
hokusai_palette <- met.brewer("Hokusai1", 8)

# Plot the dendrogram with the specified color palette
f <- fviz_dend(res.hcpc, 
               cex = 0.4, 
               palette = "jco", 
               labels_track_height = 1.0, 
               labels_cols = TRUE) +
  scale_color_manual(values = hokusai_palette) +
  theme_minimal() +
  theme(legend.position = "none")


# Define a custom color palette for the types in the color bar
my_custom_palette <- c("Type I" = "red",
                       "Type II" = "orange",
                       "Type III" = "green",
                       "Type III / Hemp" = "darkgreen",
                       "Landrace" = "blue",
                       "unknown" = "grey",
                       "Hemp" = "turquoise")

# Ensure that the BAR column is a factor with levels corresponding to the colors defined
order_data$BAR <- factor(order_data$BAR, levels = names(my_custom_palette))

# Create the color bar dataframe
color_bar_df <- data.frame(Type = order_data$BAR, x = 1:nrow(order_data))

# Adjust the y position to move the bar up
bar_position_y <- -0.5  # Adjust this value as needed to move the bar closer or further away

# Add the color bar with adjusted height
bar_height <- 0.2
f <- f + geom_tile(data = color_bar_df, aes(x = x, y = -0.5, fill = Type), 
                   width = 1, height = bar_height) +
  scale_fill_manual(values = my_custom_palette, 
                    name = "Type") +  # Include the legend
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.margin = margin(t = 1, b = 1, unit = "cm"),
        legend.position = "right")  # Adjust legend position if needed

# Plot the dendrogram with the color bar
print(f)


# Save the plot to a PDF file
ggsave("sunrise_HCPC_FIN.pdf", plot = f, width = 20, height = 9, device = "pdf")


#########################################
#SUPPLEMENTARY FIGURE 3C  #Sunrise fastSTRUCTURE
#########################################
# Load the Hokusai palette
hiroshige_palette <- met.brewer("Hiroshige", 8)

selected_colors <- hiroshige_palette[c(1, 6)]
###########
# K2 Plot
###########
test <- read.csv("~/R/Cannabis_Genetics_Manuscript/fastSTRUCTURE/plots/sunrise_K2.csv")
long_data <- pivot_longer(test, cols = c("K1", "K2"), names_to = "category", values_to = "value")

K2 <- ggplot(long_data, aes(x = Name, y = value, fill = category, group = Type)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = selected_colors) +  # Use selected colors from Hokusai palette
  labs(title = "", x = "", y = "Percent Identity") +
  facet_grid(~ Type, scales = "free_x") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 2),
    axis.text.y = element_text(size = rel(0.8)),
    plot.margin = margin(5.5, 5.5, 5.5, 5.5, "pt"),
    strip.background = element_blank(),
    legend.position = "right"
  )
K2

ggsave("sunrise_K2_stru.pdf", plot = K2, width = 20, height = 3, device = "pdf")

###########
# K3 Plot
###########
hiroshige_palette <- met.brewer("Hiroshige", 8)
selected_colors <- hiroshige_palette[c(1, 4, 6)]


test <- read.csv("~/R/Cannabis_Genetics_Manuscript/fastSTRUCTURE/plots/sunrise_K3.csv")
long_data <- pivot_longer(test, cols = c("K1", "K2", "K3"), names_to = "category", values_to = "value")

K3 <- ggplot(long_data, aes(x = Name, y = value, fill = category, group = Type)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = selected_colors) +  # Use selected colors from Hokusai palette
  labs(title = "", x = "", y = "Percent Identity") +
  facet_grid(~ Type, scales = "free_x") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 0.5),
    axis.text.y = element_text(size = rel(0.8)),
    plot.margin = margin(5.5, 5.5, 5.5, 5.5, "pt"),
    strip.background = element_blank(),
    legend.position = "right"
  )
K3

ggsave("sunrise_K3_stru.pdf", plot = K3, width = 20, height = 3, device = "pdf")

###########
# K4 Plot
###########
hiroshige_palette <- met.brewer("Hiroshige", 8)
selected_colors <- hiroshige_palette[c(1, 4,6, 8)]

test <- read.csv("~/R/Cannabis_Genetics_Manuscript/fastSTRUCTURE/plots/sunrise_K4.csv")
long_data <- pivot_longer(test, cols = c("K1", "K2", "K3", "K4"), names_to = "category", values_to = "value")

K4 <- ggplot(long_data, aes(x = Name, y = value, fill = category, group = Type)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = selected_colors) +  # Use selected colors from Hokusai palette
  labs(title = "", x = "", y = "Percent Identity") +
  facet_grid(~ Type, scales = "free_x") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 0.5),
    axis.text.y = element_text(size = rel(0.8)),
    plot.margin = margin(5.5, 5.5, 5.5, 5.5, "pt"),
    strip.background = element_blank(),
    legend.position = "right"
  )
K4

ggsave("sunrise_K4_stru.pdf", plot = K4, width = 20, height = 3, device = "pdf")

###########
# K5 Plot
###########
hiroshige_palette <- met.brewer("Hiroshige", 8)
selected_colors <- hiroshige_palette[c(1, 3,5,7, 8)]

test <- read.csv("~/R/Cannabis_Genetics_Manuscript/fastSTRUCTURE/plots/sunrise_K5.csv")
long_data <- pivot_longer(test, cols = c("K1", "K2", "K3", "K4", "K5"), names_to = "category", values_to = "value")

K5 <- ggplot(long_data, aes(x = Name, y = value, fill = category, group = Type)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = selected_colors) +  # Use selected colors from Hokusai palette
  labs(title = "", x = "", y = "Percent Identity") +
  facet_grid(~ Type, scales = "free_x") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 0.5),
    axis.text.y = element_text(size = rel(0.8)),
    plot.margin = margin(5.5, 5.5, 5.5, 5.5, "pt"),
    strip.background = element_blank(),
    legend.position = "right"
  )
K5

ggsave("sunrise_K5_stru.pdf", plot = K5, width = 20, height = 3, device = "pdf")

#########################################
#SUPPLEMENTARY FIGURE 4A - Colorado PCA
#########################################

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
pca <- read.csv("~/R/Cannabis_Genetics_Manuscript/csv_files/colorado_n162.csv")

p <- ggplot(data = pca, aes(x=EV1, y=EV2,group=Type,shape=Type,color=Type)) +
  geom_point(size=4) +
  scale_shape_manual(values = seq(0, 15)) +
  scale_color_manual(values = c("turquoise", "blue","red","orange", "green","grey")) +
  stat_ellipse(size = 1.5) +  # Adjust ellipse line thickness
  theme_classic() +
  xlab("PC 1 (5.40%)") + 
  ylab("PC 2 (3.64%)") +
  labs(
    shape = "Type",  # Custom legend title for shape
    color = "Type"   # Custom legend title for color
  ) +
  theme(
    axis.text = element_text(size = 25),  # Increase axis text size
    axis.title = element_text(size = 25),  # Increase axis title size
    legend.text = element_text(size = 25),  # Increase legend text size
    legend.title = element_text(size = 25)  # Increase legend title size
  )
p
ggsave("colorado_PCA.pdf", plot = p, width = 20, height = 13, device = "pdf")


#########################################
#SUPPLEMENTARY FIGURE 4B - Colorado Hierarchical clustering
#########################################

pca <- read.csv("~/R/Cannabis_Genetics_Manuscript/csv_files/colorado_n162.csv")
# Process the data for PCA
pca2 <- pca[3:6]
row.names(pca2) <- pca$HCPC

# Perform PCA
res.pca <- PCA(pca2, graph = FALSE)

# Perform hierarchical clustering on the principal components
res.hcpc <- HCPC(res.pca, graph = FALSE)

# Read the order data for the color bar
order_data <- read.csv("~/R/Cannabis_Genetics_Manuscript/csv_files/colorado_HCPC_order_plots.csv")

# Define the color palette for the dendrogram
hokusai_palette <- met.brewer("Hokusai1", 8)

# Plot the dendrogram with the specified color palette
f <- fviz_dend(res.hcpc, 
               cex = 0.4, 
               palette = "jco", 
               labels_track_height = 1.0, 
               labels_cols = TRUE) +
  scale_color_manual(values = hokusai_palette) +
  theme_minimal() +
  theme(legend.position = "none")

# Define a custom color palette for the types in the color bar
my_custom_palette <- c("Type I" = "red",
                       "Type II" = "orange",
                       "Type III" = "green",
                       "Type III / Hemp" = "darkgreen",
                       "Landrace" = "blue",
                       "unknown" = "grey",
                       "Hemp" = "turquoise")

# Ensure that the BAR column is a factor with levels corresponding to the colors defined
order_data$BAR <- factor(order_data$BAR, levels = names(my_custom_palette))

# Create the color bar dataframe
color_bar_df <- data.frame(Type = order_data$BAR, x = 1:nrow(order_data))

# Adjust the y position to move the bar up
bar_position_y <- -0.5  # Adjust this value as needed to move the bar closer or further away

# Add the color bar with adjusted height
bar_height <- 0.2
f <- f + geom_tile(data = color_bar_df, aes(x = x, y = -0.5, fill = Type), 
                   width = 1, height = bar_height) +
  scale_fill_manual(values = my_custom_palette, 
                    name = "Type") +  # Include the legend
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.margin = margin(t = 1, b = 1, unit = "cm"),
        legend.position = "right")  # Adjust legend position if needed

# Plot the dendrogram with the color bar
print(f)


# Save the plot to a PDF file
ggsave("colorado_HCPC_FIN.pdf", plot = f, width = 20, height = 9, device = "pdf")

#########################################
#SUPPLEMENTARY FIGURE 4C - Colorado fastSTRUCTURE
#########################################
# Load the Hokusai palette
hiroshige_palette <- met.brewer("Hiroshige", 8)

selected_colors <- hiroshige_palette[c(1, 6)]
###########
# K2 Plot
###########
test <- read.csv("~/R/Cannabis_Genetics_Manuscript/fastSTRUCTURE/plots/colorado_K2.csv")
long_data <- pivot_longer(test, cols = c("K1", "K2"), names_to = "category", values_to = "value")

K2 <- ggplot(long_data, aes(x = Name, y = value, fill = category, group = Type)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = selected_colors) +  # Use selected colors from Hokusai palette
  labs(title = "", x = "", y = "Percent Identity") +
  facet_grid(~ Type, scales = "free_x") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 2),
    axis.text.y = element_text(size = rel(0.8)),
    plot.margin = margin(5.5, 5.5, 5.5, 5.5, "pt"),
    strip.background = element_blank(),
    legend.position = "right"
  )
K2

ggsave("colorado_K2_stru.pdf", plot = K2, width = 20, height = 3, device = "pdf")

###########
# K3 Plot
###########
hiroshige_palette <- met.brewer("Hiroshige", 8)
selected_colors <- hiroshige_palette[c(1, 4, 6)]

test <- read.csv("~/R/Cannabis_Genetics_Manuscript/fastSTRUCTURE/plots/colorado_K3.csv")
long_data <- pivot_longer(test, cols = c("K1", "K2", "K3"), names_to = "category", values_to = "value")

K3 <- ggplot(long_data, aes(x = Name, y = value, fill = category, group = Type)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = selected_colors) +  # Use selected colors from Hokusai palette
  labs(title = "", x = "", y = "Percent Identity") +
  facet_grid(~ Type, scales = "free_x") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 0.5),
    axis.text.y = element_text(size = rel(0.8)),
    plot.margin = margin(5.5, 5.5, 5.5, 5.5, "pt"),
    strip.background = element_blank(),
    legend.position = "right"
  )
K3

ggsave("colorado_K3_stru.pdf", plot = K3, width = 20, height = 3, device = "pdf")

###########
# K4 Plot
###########
hiroshige_palette <- met.brewer("Hiroshige", 8)
selected_colors <- hiroshige_palette[c(1, 4,6, 8)]

test <- read.csv("~/R/Cannabis_Genetics_Manuscript/fastSTRUCTURE/plots/colorado_K4.csv")
long_data <- pivot_longer(test, cols = c("K1", "K2", "K3", "K4"), names_to = "category", values_to = "value")

K4 <- ggplot(long_data, aes(x = Name, y = value, fill = category, group = Type)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = selected_colors) +  # Use selected colors from Hokusai palette
  labs(title = "", x = "", y = "Percent Identity") +
  facet_grid(~ Type, scales = "free_x") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 0.5),
    axis.text.y = element_text(size = rel(0.8)),
    plot.margin = margin(5.5, 5.5, 5.5, 5.5, "pt"),
    strip.background = element_blank(),
    legend.position = "right"
  )
K4

ggsave("colorado_K4_stru.pdf", plot = K4, width = 20, height = 3, device = "pdf")

###########
# K5 Plot
###########
hiroshige_palette <- met.brewer("Hiroshige", 8)
selected_colors <- hiroshige_palette[c(1, 3,5,7, 8)]

test <- read.csv("~/R/Cannabis_Genetics_Manuscript/fastSTRUCTURE/plots/colorado_K5.csv")
long_data <- pivot_longer(test, cols = c("K1", "K2", "K3", "K4", "K5"), names_to = "category", values_to = "value")

K5 <- ggplot(long_data, aes(x = Name, y = value, fill = category, group = Type)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = selected_colors) +  # Use selected colors from Hokusai palette
  labs(title = "", x = "", y = "Percent Identity") +
  facet_grid(~ Type, scales = "free_x") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 0.5),
    axis.text.y = element_text(size = rel(0.8)),
    plot.margin = margin(5.5, 5.5, 5.5, 5.5, "pt"),
    strip.background = element_blank(),
    legend.position = "right"
  )
K5

ggsave("colorado_K5_stru.pdf", plot = K5, width = 20, height = 3, device = "pdf")



#########################################
#SUPPLEMENTARY FIGURE 5A - Courtagen PCA
#########################################

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

#Plot 2D PCA 
pca <- read.csv("~/R/Cannabis_Genetics_Manuscript/csv_files/courtagen_n58_typeknown.csv")
p <- ggplot(data = pca, aes(x=EV1, y=EV2,group=Type,shape=Type,color=Type)) +
  geom_point(size=4) +
  scale_shape_manual(values = seq(0, 15)) +
  scale_color_manual(values = c("turquoise", "red","green","grey")) +
  stat_ellipse(size = 1.5) +  # Adjust ellipse line thickness
  theme_classic() +
  xlab("PC 1 (6.16%)") + 
  ylab("PC 2 (4.83%)") +
  labs(
    shape = "Type",  # Custom legend title for shape
    color = "Type"   # Custom legend title for color
  ) +
  theme(
    axis.text = element_text(size = 25),  # Increase axis text size
    axis.title = element_text(size = 25),  # Increase axis title size
    legend.text = element_text(size = 25),  # Increase legend text size
    legend.title = element_text(size = 25)  # Increase legend title size
  )
p
ggsave("courtagen_PCA.pdf", plot = p, width = 20, height = 13, device = "pdf")

#########################################
#SUPPLEMENTARY FIGURE 5B - Courtagen Hierarchical Cluster
#########################################

pca <- read.csv("~/R/Cannabis_Genetics_Manuscript/csv_files/courtagen_n58_typeknown.csv")
# Process the data for PCA
pca2 <- pca[3:6]
row.names(pca2) <- pca$HCPC

# Perform PCA
res.pca <- PCA(pca2, graph = FALSE)

# Perform hierarchical clustering on the principal components
res.hcpc <- HCPC(res.pca, graph = FALSE)

# Read the order data for the color bar
order_data <- read.csv("~/R/Cannabis_Genetics_Manuscript/csv_files/courtagen_HCPC_order_plots.csv")

# Define the color palette for the dendrogram
hokusai_palette <- met.brewer("Hokusai1", 8)

# Plot the dendrogram with the specified color palette
f <- fviz_dend(res.hcpc, 
               cex = 0.4, 
               palette = "jco", 
               labels_track_height = 1.0, 
               labels_cols = TRUE) +
  scale_color_manual(values = hokusai_palette) +
  theme_minimal() +
  theme(legend.position = "none")

# Define a custom color palette for the types in the color bar
my_custom_palette <- c("Type I" = "red",
                       "Type II" = "orange",
                       "Type III" = "green",
                       "Type III / Hemp" = "darkgreen",
                       "Landrace" = "blue",
                       "unknown" = "grey",
                       "Hemp" = "turquoise")

# Ensure that the BAR column is a factor with levels corresponding to the colors defined
order_data$BAR <- factor(order_data$BAR, levels = names(my_custom_palette))

# Create the color bar dataframe
color_bar_df <- data.frame(Type = order_data$BAR, x = 1:nrow(order_data))

# Adjust the y position to move the bar up
bar_position_y <- -0.5  # Adjust this value as needed to move the bar closer or further away

# Add the color bar with adjusted height
bar_height <- 0.2
f <- f + geom_tile(data = color_bar_df, aes(x = x, y = -0.5, fill = Type), 
                   width = 1, height = bar_height) +
  scale_fill_manual(values = my_custom_palette, 
                    name = "Type") +  # Include the legend
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.margin = margin(t = 1, b = 1, unit = "cm"),
        legend.position = "right")  # Adjust legend position if needed

# Plot the dendrogram with the color bar
print(f)

# Save the plot to a PDF file
ggsave("courtagen_HCPC_FIN.pdf", plot = f, width = 20, height = 9, device = "pdf")

#########################################
#SUPPLEMENTARY FIGURE 5C - COURTAGEN FASTSTRUCTURE
#########################################

# Load the Hokusai palette
hiroshige_palette <- met.brewer("Hiroshige", 8)

selected_colors <- hiroshige_palette[c(1, 6)]
###########
# K2 Plot
###########
test <- read.csv("~/R/Cannabis_Genetics_Manuscript/fastSTRUCTURE/plots/courtagen_K2.csv")
long_data <- pivot_longer(test, cols = c("K1", "K2"), names_to = "category", values_to = "value")

K2 <- ggplot(long_data, aes(x = Name, y = value, fill = category, group = Type)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = selected_colors) +  # Use selected colors from Hokusai palette
  labs(title = "", x = "", y = "Percent Identity") +
  facet_grid(~ Type, scales = "free_x") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 2),
    axis.text.y = element_text(size = rel(0.8)),
    plot.margin = margin(5.5, 5.5, 5.5, 5.5, "pt"),
    strip.background = element_blank(),
    legend.position = "right"
  )
K2

ggsave("courtagen_K2_stru.pdf", plot = K2, width = 20, height = 3, device = "pdf")

###########
# K3 Plot
###########
hiroshige_palette <- met.brewer("Hiroshige", 8)
selected_colors <- hiroshige_palette[c(1, 4, 6)]


test <- read.csv("~/R/Cannabis_Genetics_Manuscript/fastSTRUCTURE/plots/courtagen_K3.csv")
long_data <- pivot_longer(test, cols = c("K1", "K2", "K3"), names_to = "category", values_to = "value")

K3 <- ggplot(long_data, aes(x = Name, y = value, fill = category, group = Type)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = selected_colors) +  # Use selected colors from Hokusai palette
  labs(title = "", x = "", y = "Percent Identity") +
  facet_grid(~ Type, scales = "free_x") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 0.5),
    axis.text.y = element_text(size = rel(0.8)),
    plot.margin = margin(5.5, 5.5, 5.5, 5.5, "pt"),
    strip.background = element_blank(),
    legend.position = "right"
  )
K3

ggsave("courtagen_K3_stru.pdf", plot = K3, width = 20, height = 3, device = "pdf")

###########
# K4 Plot
###########
hiroshige_palette <- met.brewer("Hiroshige", 8)
selected_colors <- hiroshige_palette[c(1, 4,6, 8)]

test <- read.csv("~/R/Cannabis_Genetics_Manuscript/fastSTRUCTURE/plots/courtagen_K4.csv")
long_data <- pivot_longer(test, cols = c("K1", "K2", "K3", "K4"), names_to = "category", values_to = "value")

K4 <- ggplot(long_data, aes(x = Name, y = value, fill = category, group = Type)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = selected_colors) +  # Use selected colors from Hokusai palette
  labs(title = "", x = "", y = "Percent Identity") +
  facet_grid(~ Type, scales = "free_x") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 0.5),
    axis.text.y = element_text(size = rel(0.8)),
    plot.margin = margin(5.5, 5.5, 5.5, 5.5, "pt"),
    strip.background = element_blank(),
    legend.position = "right"
  )
K4

ggsave("courtagen_K4_stru.pdf", plot = K4, width = 20, height = 3, device = "pdf")

###########
# K5 Plot
###########
hiroshige_palette <- met.brewer("Hiroshige", 8)
selected_colors <- hiroshige_palette[c(1, 3,5,7, 8)]

test <- read.csv("~/R/Cannabis_Genetics_Manuscript/fastSTRUCTURE/plots/courtagen_K5.csv")
long_data <- pivot_longer(test, cols = c("K1", "K2", "K3", "K4", "K5"), names_to = "category", values_to = "value")

K5 <- ggplot(long_data, aes(x = Name, y = value, fill = category, group = Type)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = selected_colors) +  # Use selected colors from Hokusai palette
  labs(title = "", x = "", y = "Percent Identity") +
  facet_grid(~ Type, scales = "free_x") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 0.5),
    axis.text.y = element_text(size = rel(0.8)),
    plot.margin = margin(5.5, 5.5, 5.5, 5.5, "pt"),
    strip.background = element_blank(),
    legend.position = "right"
  )
K5

ggsave("courtagen_K5_stru.pdf", plot = K5, width = 20, height = 3, device = "pdf")

#########################################
#SUPPLEMENTARY FIGURE 6A - Kannapedia PCA
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
pca <- read.csv("~/R/Cannabis_Genetics_Manuscript/csv_files/kannapedia_n61_typeknown.csv")


p <- ggplot(data = pca, aes(x=EV1, y=EV2,group=Type,shape=Type,color=Type)) +
  geom_point(size=4) +
  scale_shape_manual(values = seq(0, 15)) +
  scale_color_manual(values = c("red","green", "grey")) +
  stat_ellipse(size = 1.5) +  # Adjust ellipse line thickness
  theme_classic() +
  xlab("PC 1 (4.95%)") + 
  ylab("PC 2 (3.89%)") +
  labs(
    shape = "Type",  # Custom legend title for shape
    color = "Type"   # Custom legend title for color
  ) +
  theme(
    axis.text = element_text(size = 25),  # Increase axis text size
    axis.title = element_text(size = 25),  # Increase axis title size
    legend.text = element_text(size = 25),  # Increase legend text size
    legend.title = element_text(size = 25)  # Increase legend title size
  )
p
ggsave("kannapedia_PCA.pdf", plot = p, width = 20, height = 13, device = "pdf")

#########################################
#SUPPLEMENTARY FIGURE 6B - Kannapedia Hierarchical Cluster
#########################################

pca <- read.csv("~/R/Cannabis_Genetics_Manuscript/csv_files/kannapedia_n61_typeknown.csv")
# Process the data for PCA
pca2 <- pca[3:6]
row.names(pca2) <- pca$HCPC

# Perform PCA
res.pca <- PCA(pca2, graph = FALSE)

# Perform hierarchical clustering on the principal components
res.hcpc <- HCPC(res.pca, graph = FALSE)

# Read the order data for the color bar
order_data <- read.csv("~/R/Cannabis_Genetics_Manuscript/csv_files/kannapedia_HCPC_order_plots.csv")

# Define the color palette for the dendrogram
hokusai_palette <- met.brewer("Hokusai1", 8)

# Plot the dendrogram with the specified color palette
f <- fviz_dend(res.hcpc, 
               cex = 0.4, 
               palette = "jco", 
               labels_track_height = 1.0, 
               labels_cols = TRUE) +
  scale_color_manual(values = hokusai_palette) +
  theme_minimal() +
  theme(legend.position = "none")

# Define a custom color palette for the types in the color bar
my_custom_palette <- c("Type I" = "red",
                       "Type II" = "orange",
                       "Type III" = "green",
                       "Type III / Hemp" = "darkgreen",
                       "Landrace" = "blue",
                       "unknown" = "grey",
                       "Hemp" = "turquoise")

# Ensure that the BAR column is a factor with levels corresponding to the colors defined
order_data$BAR <- factor(order_data$BAR, levels = names(my_custom_palette))

# Create the color bar dataframe
color_bar_df <- data.frame(Type = order_data$BAR, x = 1:nrow(order_data))

# Adjust the y position to move the bar up
bar_position_y <- -0.5  # Adjust this value as needed to move the bar closer or further away

# Add the color bar with adjusted height
bar_height <- 0.2
f <- f + geom_tile(data = color_bar_df, aes(x = x, y = -0.5, fill = Type), 
                   width = 1, height = bar_height) +
  scale_fill_manual(values = my_custom_palette, 
                    name = "Type") +  # Include the legend
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.margin = margin(t = 1, b = 1, unit = "cm"),
        legend.position = "right")  # Adjust legend position if needed

# Plot the dendrogram with the color bar
print(f)

# Save the plot to a PDF file
ggsave("kannapedia_HCPC_FIN.pdf", plot = f, width = 20, height = 9, device = "pdf")


#########################################
#SUPPLEMENTARY FIGURE 6C - Kannapedia fastSTRUCTURE
#########################################
# Load the Hokusai palette
hiroshige_palette <- met.brewer("Hiroshige", 8)

selected_colors <- hiroshige_palette[c(1, 6)]
###########
# K2 Plot
###########
test <- read.csv("~/R/Cannabis_Genetics_Manuscript/fastSTRUCTURE/plots/kannapedia_K2.csv")
long_data <- pivot_longer(test, cols = c("K1", "K2"), names_to = "category", values_to = "value")

K2 <- ggplot(long_data, aes(x = Name, y = value, fill = category, group = Type)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = selected_colors) +  # Use selected colors from Hokusai palette
  labs(title = "", x = "", y = "Percent Identity") +
  facet_grid(~ Type, scales = "free_x") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 2),
    axis.text.y = element_text(size = rel(0.8)),
    plot.margin = margin(5.5, 5.5, 5.5, 5.5, "pt"),
    strip.background = element_blank(),
    legend.position = "right"
  )
K2

ggsave("kannapedia_K2_stru.pdf", plot = K2, width = 20, height = 3, device = "pdf")

###########
# K3 Plot
###########
hiroshige_palette <- met.brewer("Hiroshige", 8)
selected_colors <- hiroshige_palette[c(1, 4, 6)]


test <- read.csv("~/R/Cannabis_Genetics_Manuscript/fastSTRUCTURE/plots/kannapedia_K3.csv")
long_data <- pivot_longer(test, cols = c("K1", "K2", "K3"), names_to = "category", values_to = "value")

K3 <- ggplot(long_data, aes(x = Name, y = value, fill = category, group = Type)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = selected_colors) +  # Use selected colors from Hokusai palette
  labs(title = "", x = "", y = "Percent Identity") +
  facet_grid(~ Type, scales = "free_x") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 0.5),
    axis.text.y = element_text(size = rel(0.8)),
    plot.margin = margin(5.5, 5.5, 5.5, 5.5, "pt"),
    strip.background = element_blank(),
    legend.position = "right"
  )
K3

ggsave("kannapedia_K3_stru.pdf", plot = K3, width = 20, height = 3, device = "pdf")

###########
# K4 Plot
###########
hiroshige_palette <- met.brewer("Hiroshige", 8)
selected_colors <- hiroshige_palette[c(1, 4,6, 8)]

test <- read.csv("~/R/Cannabis_Genetics_Manuscript/fastSTRUCTURE/plots/kannapedia_K4.csv")
long_data <- pivot_longer(test, cols = c("K1", "K2", "K3", "K4"), names_to = "category", values_to = "value")

K4 <- ggplot(long_data, aes(x = Name, y = value, fill = category, group = Type)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = selected_colors) +  # Use selected colors from Hokusai palette
  labs(title = "", x = "", y = "Percent Identity") +
  facet_grid(~ Type, scales = "free_x") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 0.5),
    axis.text.y = element_text(size = rel(0.8)),
    plot.margin = margin(5.5, 5.5, 5.5, 5.5, "pt"),
    strip.background = element_blank(),
    legend.position = "right"
  )
K4

ggsave("kannapedia_K4_stru.pdf", plot = K4, width = 20, height = 3, device = "pdf")

###########
# K5 Plot
###########
hiroshige_palette <- met.brewer("Hiroshige", 8)
selected_colors <- hiroshige_palette[c(1, 3,5,7, 8)]

test <- read.csv("~/R/Cannabis_Genetics_Manuscript/fastSTRUCTURE/plots/kannapedia_K5.csv")
long_data <- pivot_longer(test, cols = c("K1", "K2", "K3", "K4", "K5"), names_to = "category", values_to = "value")

K5 <- ggplot(long_data, aes(x = Name, y = value, fill = category, group = Type)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = selected_colors) +  # Use selected colors from Hokusai palette
  labs(title = "", x = "", y = "Percent Identity") +
  facet_grid(~ Type, scales = "free_x") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 0.5),
    axis.text.y = element_text(size = rel(0.8)),
    plot.margin = margin(5.5, 5.5, 5.5, 5.5, "pt"),
    strip.background = element_blank(),
    legend.position = "right"
  )
K5

ggsave("kannapedia_K5_stru.pdf", plot = K5, width = 20, height = 3, device = "pdf")

#########################################
#SUPPLEMENTARY FIGURE 9/10
#########################################
library("FactoMineR")
library("factoextra")

#pca <- read.csv("~/R/Cannabis_Genetics_Manuscript/csv_files/lw_n498_typeknown.csv")
#pca <- read.csv("~/R/Cannabis_Genetics_Manuscript/csv_files/phylos_845_PRJNA347566_typeknown_hcpc.csv")
#pca <- read.csv("~/R/Cannabis_Genetics_Manuscript/csv_files/soorni_n94.csv")
#pca <- read.csv("~/R/Cannabis_Genetics_Manuscript/csv_files/MG_V1_n289.csv")

#pca <- read.csv("~/R/Cannabis_Genetics_Manuscript/csv_files/phylos_1378_PRJNA510566_typeknown.csv")
#pca <- read.csv("~/R/Cannabis_Genetics_Manuscript/csv_files/sunrise_n25.csv")
#pca <- read.csv("~/R/Cannabis_Genetics_Manuscript/csv_files/colorado_n162.csv")
#pca <- read.csv("~/R/Cannabis_Genetics_Manuscript/csv_files/courtagen_n58_typeknown.csv")
#pca <- read.csv("~/R/Cannabis_Genetics_Manuscript/csv_files/kannapedia_n61_typeknown.csv")

pca <- read.csv("~/R/Cannabis_Genetics_Manuscript/csv_files/lw_nuclear_lr_subset_only.csv")
head(pca)

#Select 3rd and 4th principal components
pca2 <- pca[3:4]

#select 3rd to 6th principal components
#pca2 <- pca[3:6]

head(pca2)
# Standardizing the PCA data before clustering
my_data <- scale(pca2)

#elbow 
a <- fviz_nbclust(my_data, kmeans, method = "wss") + ggtitle("the Elbow Method")

#ggsave("lw_elbow.pdf", plot = a, width = 6, height = 3, device = "pdf")
#ggsave("phylos845_elbow.pdf", plot = a, width = 6, height = 3, device = "pdf")
#ggsave("soorni_elbow.pdf", plot = a, width = 6, height = 3, device = "pdf")
#ggsave("mgv1_elbow.pdf", plot = a, width = 6, height = 3, device = "pdf")
#ggsave("phylos1378_elbow.pdf", plot = a, width = 6, height = 3, device = "pdf")
#ggsave("sunrise_elbow.pdf", plot = a, width = 6, height = 3, device = "pdf")
#ggsave("colorado_elbow.pdf", plot = a, width = 6, height = 3, device = "pdf")
#ggsave("courtagen_elbow.pdf", plot = a, width = 6, height = 3, device = "pdf")
#ggsave("kannapedia_elbow.pdf", plot = a, width = 6, height = 3, device = "pdf")
ggsave("lr_subset_elbow.pdf", plot = a, width = 6, height = 3, device = "pdf")

#silhouette
b<- fviz_nbclust(my_data, kmeans, method = "silhouette") + ggtitle("The Silhouette Plot")

#ggsave("lw_silhouette.pdf", plot = b, width = 6, height = 3, device = "pdf")
#ggsave("phylos845_silhouette.pdf", plot = b, width = 6, height = 3, device = "pdf")
#ggsave("soorni_silhouette.pdf", plot = b, width = 6, height = 3, device = "pdf")
#ggsave("mgv1_silhouette.pdf", plot = b, width = 6, height = 3, device = "pdf")
#ggsave("phylos1378_silhouette.pdf", plot = b, width = 6, height = 3, device = "pdf")
#ggsave("sunrise_silhouette.pdf", plot = b, width = 6, height = 3, device = "pdf")
#ggsave("colorado_silhouette.pdf", plot = b, width = 6, height = 3, device = "pdf")
#ggsave("courtagen_silhouette.pdf", plot = b, width = 6, height = 3, device = "pdf")
#ggsave("kannapedia_silhouette.pdf", plot = b, width = 6, height = 3, device = "pdf")
ggsave("lr_subset_silhouette.pdf", plot = b, width = 6, height = 3, device = "pdf")

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
setwd("~/R/Cannabis_Genetics_Manuscript/Filtered_VCF_Files/")

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

######################## FIN
