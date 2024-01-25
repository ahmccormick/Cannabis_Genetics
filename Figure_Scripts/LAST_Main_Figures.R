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

#####################
#Figure 1A - HCPC LeafWorks
#####################
pca <- read.csv("~/R/Cannabis_Genetics_Manuscript/csv_files/lw_n498_typeknown.csv")

# Process the data for PCA
pca2 <- pca[3:6]
row.names(pca2) <- pca$HCPC

# Perform PCA
res.pca <- PCA(pca2, graph = FALSE)

# Perform hierarchical clustering on the principal components
res.hcpc <- HCPC(res.pca, graph = FALSE)

# Read the order data for the color bar
order_data <- read.csv("~/R/Cannabis_Genetics_Manuscript/csv_files/LW_HCPC_order_plots.csv")

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
                       "Hemp" = "turquoise",
                       "Landrace" = "blue",
                       "Unknown" = "grey")

# Ensure that the BAR column is a factor with levels corresponding to the colors defined
order_data$BAR <- factor(order_data$BAR, levels = names(my_custom_palette))

# Create the color bar dataframe
color_bar_df <- data.frame(Type = order_data$BAR, x = 1:nrow(order_data))

# Adjust the y position to move the bar up
bar_position_y <- -0.5  # Adjust this value as needed to move the bar closer or further away

# Add the color bar with adjusted height
bar_height <- 0.2
f <- f + geom_tile(data = color_bar_df, aes(x = x, y = bar_position_y, fill = Type), 
                   width = 1, height = bar_height) +
  scale_fill_manual(values = my_custom_palette) +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.margin = margin(t = 1, b = 1, unit = "cm"))

# Plot the dendrogram with the color bar
print(f)

# Save the plot to a PDF file
ggsave("LW_HCPC_FIN.pdf", plot = f, width = 20, height = 9, device = "pdf")

#####################
#Figure 1B - fastSTRUCTURE LeafWorks
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
test <- read.csv("~/R/Cannabis_Genetics_Manuscript/fastSTRUCTURE/plots/LeafK2.csv")
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

ggsave("LWK2_stru.pdf", plot = K2, width = 20, height = 3, device = "pdf")

###########
# K3 Plot
###########
hiroshige_palette <- met.brewer("Hiroshige", 8)
selected_colors <- hiroshige_palette[c(1, 4, 6)]


test <- read.csv("~/R/Cannabis_Genetics_Manuscript/fastSTRUCTURE/plots/LeafK3.csv")
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

ggsave("LWK3_stru.pdf", plot = K3, width = 20, height = 3, device = "pdf")

###########
# K4 Plot
###########
hiroshige_palette <- met.brewer("Hiroshige", 8)
selected_colors <- hiroshige_palette[c(1, 4,6, 8)]

test <- read.csv("~/R/Cannabis_Genetics_Manuscript/fastSTRUCTURE/plots/LeafK4.csv")
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

ggsave("LWK4_stru.pdf", plot = K4, width = 20, height = 3, device = "pdf")

###########
# K5 Plot
###########
hiroshige_palette <- met.brewer("Hiroshige", 8)
selected_colors <- hiroshige_palette[c(1, 3,5,7, 8)]

test <- read.csv("~/R/Cannabis_Genetics_Manuscript/fastSTRUCTURE/plots/LeafK5.csv")
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

ggsave("LWK5_stru.pdf", plot = K5, width = 20, height = 3, device = "pdf")

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

#####################
#Figure 1C - Phylos
#####################
# Read the PCA data
pca <- read.csv("~/R/Cannabis_Genetics_Manuscript/csv_files/phylos_845_PRJNA347566_typeknown_hcpc.csv")

# Process the data for PCA
pca2 <- pca[3:6]
row.names(pca2) <- pca$Cultivar_new

# Perform PCA
res.pca <- PCA(pca2, graph = FALSE)

# Perform hierarchical clustering on the principal components
res.hcpc <- HCPC(res.pca, graph = FALSE)

# Read the order data for the color bar
order_data <- read.csv("~/R/Cannabis_Genetics_Manuscript/csv_files/Phylos845_HCPC_order_plots.csv")

# Define the color palette for the dendrogram
hokusai_palette <- met.brewer("Hokusai1", 8)

# Plot the dendrogram with the specified color palette
f <- fviz_dend(res.hcpc, 
               cex = 0.06, 
               palette = "jco", 
               labels_track_height = 1.0, 
               labels_cols = TRUE) +
  scale_color_manual(values = hokusai_palette) +
  theme_minimal()

# Define a custom color palette for the types in the color bar
my_custom_palette <- c("Type I" = "red",
                       "Type II" = "orange",
                       "Type III" = "green",
                       "Hemp" = "turquoise",
                       "Landrace" = "blue",
                       "Unknown" = "grey")

# Ensure that the BAR column is a factor with levels corresponding to the colors defined
order_data$BAR <- factor(order_data$BAR, levels = names(my_custom_palette))

# Create the color bar dataframe
color_bar_df <- data.frame(Type = order_data$BAR, x = 1:nrow(order_data))

# Add the color bar with adjusted height
bar_height <- 0.2
f <- f + geom_tile(data = color_bar_df, aes(x = x, y = -0.5, fill = Type), 
                   width = 1, height = bar_height) +
  scale_fill_manual(values = my_custom_palette, 
                    name = "Use Type") +  # Include the legend
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
ggsave("Phylos845_HCPC_FIN.pdf", plot = f, width = 20, height = 9, device = "pdf")

#####################
#Figure 1D - fastSTRUCTURE Phylos 845
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
test<-read.csv("~/R/Cannabis_Genetics_Manuscript/fastSTRUCTURE/plots/Phylos_K2.csv")

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

ggsave("Phylos845K2_stru.pdf", plot = K2, width = 20, height = 3, device = "pdf")

###########
# K3 Plot
###########
hiroshige_palette <- met.brewer("Hiroshige", 8)
selected_colors <- hiroshige_palette[c(1, 4, 6)]

test<-read.csv("~/R/Cannabis_Genetics_Manuscript/fastSTRUCTURE/plots/Phylos_K3.csv")

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

ggsave("Phylos845K3_stru.pdf", plot = K3, width = 20, height = 3, device = "pdf")

###########
# K4 Plot
###########
hiroshige_palette <- met.brewer("Hiroshige", 8)
selected_colors <- hiroshige_palette[c(1, 4,6, 8)]

test<-read.csv("~/R/Cannabis_Genetics_Manuscript/fastSTRUCTURE/plots/Phylos_K4.csv")
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

ggsave("Phylos845K4_stru.pdf", plot = K4, width = 20, height = 3, device = "pdf")

###########
# K5 Plot
###########
hiroshige_palette <- met.brewer("Hiroshige", 8)
selected_colors <- hiroshige_palette[c(1, 3,5,7, 8)]

test<-read.csv("~/R/Cannabis_Genetics_Manuscript/fastSTRUCTURE/plots/Phylos_K5.csv")
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

ggsave("Phylos845K5_stru.pdf", plot = K5, width = 20, height = 3, device = "pdf")

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

#####################
#Figure 2A - HCPC Soorni
#####################
pca <- read.csv("~/R/Cannabis_Genetics_Manuscript/csv_files/soorni_n94_typeknown_hcpc.csv")
# Process the data for PCA
pca2 <- pca[3:6]
row.names(pca2) <- pca$HCPC

# Perform PCA
res.pca <- PCA(pca2, graph = FALSE)

# Perform hierarchical clustering on the principal components
res.hcpc <- HCPC(res.pca, graph = FALSE)

# Read the order data for the color bar
order_data <- read.csv("~/R/Cannabis_Genetics_Manuscript/csv_files/Soorni_HCPC_order_plots.csv")

# Define the color palette for the dendrogram
hokusai_palette <- met.brewer("Hokusai1", 8)

# Plot the dendrogram with the specified color palette
f <- fviz_dend(res.hcpc, 
               cex = 0.4, 
               palette = "jco", 
               labels_track_height = 1.0, 
               labels_cols = TRUE) +
  scale_color_manual(values = hokusai_palette) +
  theme_minimal()
#f


# Define a custom color palette for the types in the color bar
my_custom_palette <- c("Type I" = "red",
                       "Type II" = "orange",
                       "Type III" = "green",
                       "Hemp" = "turquoise",
                       "Landrace" = "blue",
                       "Unknown" = "grey")

# Ensure that the BAR column is a factor with levels corresponding to the colors defined
order_data$BAR <- factor(order_data$BAR, levels = names(my_custom_palette))

# Create the color bar dataframe
color_bar_df <- data.frame(Type = order_data$BAR, x = 1:nrow(order_data))

# Add the color bar with adjusted height
bar_height <- 0.2
f <- f + geom_tile(data = color_bar_df, aes(x = x, y = -0.9, fill = Type), 
                   width = 1, height = bar_height) +
  scale_fill_manual(values = my_custom_palette, 
                    name = "Use Type") +  # Include the legend
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
ggsave("Soorni_HCPC_FIN.pdf", plot = f, width = 20, height = 9, device = "pdf")

#####################
#Figure 2B - fastSTRUCTURE Soorni
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
test<-read.csv("~/R/Cannabis_Genetics_Manuscript/fastSTRUCTURE/plots/Soorni_K2.csv")
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

ggsave("SoorniK2_stru.pdf", plot = K2, width = 20, height = 3, device = "pdf")

###########
# K3 Plot
###########
hiroshige_palette <- met.brewer("Hiroshige", 8)
selected_colors <- hiroshige_palette[c(1, 4, 6)]

test<-read.csv("~/R/Cannabis_Genetics_Manuscript/fastSTRUCTURE/plots/Soorni_K3.csv")
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

ggsave("SoorniK3_stru.pdf", plot = K3, width = 20, height = 3, device = "pdf")

###########
# K4 Plot
###########
hiroshige_palette <- met.brewer("Hiroshige", 8)
selected_colors <- hiroshige_palette[c(1, 4,6, 8)]

test<-read.csv("~/R/Cannabis_Genetics_Manuscript/fastSTRUCTURE/plots/Soorni_K4.csv")
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

ggsave("SoorniK4_stru.pdf", plot = K4, width = 20, height = 3, device = "pdf")

###########
# K5 Plot
###########
hiroshige_palette <- met.brewer("Hiroshige", 8)
selected_colors <- hiroshige_palette[c(1, 3,5,7, 8)]

test<-read.csv("~/R/Cannabis_Genetics_Manuscript/fastSTRUCTURE/plots/Soorni_K5.csv")
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

ggsave("SoorniK5_stru.pdf", plot = K5, width = 20, height = 3, device = "pdf")

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

#####################
#Figure 2C - HCPC MG V1
#####################
# Read the PCA data
pca <- read.csv("~/R/Cannabis_Genetics_Manuscript/csv_files/MG_V1_n289_typeknown_hcpc.csv")
# Process the data for PCA
pca2 <- pca[3:6]
row.names(pca2) <- pca$Cultivar_Type

# Perform PCA
res.pca <- PCA(pca2, graph = FALSE)

# Perform hierarchical clustering on the principal components
res.hcpc <- HCPC(res.pca, graph = FALSE)

# Read the order data for the color bar
order_data <- read.csv("~/R/Cannabis_Genetics_Manuscript/csv_files/MGV1_HCPC_order_plots.csv")

# Define the color palette for the dendrogram
hokusai_palette <- met.brewer("Hokusai1", 8)

# Plot the dendrogram with the specified color palette
f <- fviz_dend(res.hcpc, 
               cex = 0.2, 
               palette = "jco", 
               labels_track_height = 1.0, 
               labels_cols = TRUE) +
  scale_color_manual(values = hokusai_palette) +
  theme_minimal()
f


# Define a custom color palette for the types in the color bar
my_custom_palette <- c("Type I" = "red",
                       "Type II" = "orange",
                       "Type III" = "green",
                       "Hemp" = "turquoise",
                       "Landrace" = "blue",
                       "Unknown" = "grey")

# Ensure that the BAR column is a factor with levels corresponding to the colors defined
order_data$BAR <- factor(order_data$BAR, levels = names(my_custom_palette))

# Create the color bar dataframe
color_bar_df <- data.frame(Type = order_data$BAR, x = 1:nrow(order_data))

# Add the color bar with adjusted height
bar_height <- 0.125
f <- f + geom_tile(data = color_bar_df, aes(x = x, y = -0.5, fill = Type), 
                   width = 1, height = bar_height) +
  scale_fill_manual(values = my_custom_palette, 
                    name = "Use Type") +  # Include the legend
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
ggsave("MGV1_HCPC_FIN.pdf", plot = f, width = 20, height = 13, device = "pdf")

#####################
#Figure 2D - fastSTRUCTURE MG v1
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
test<-read.csv("~/R/Cannabis_Genetics_Manuscript/fastSTRUCTURE/plots/MG_K2.csv")
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

ggsave("MGK2_stru.pdf", plot = K2, width = 20, height = 3, device = "pdf")

###########
# K3 Plot
###########
hiroshige_palette <- met.brewer("Hiroshige", 8)
selected_colors <- hiroshige_palette[c(1, 4, 6)]

test<-read.csv("~/R/Cannabis_Genetics_Manuscript/fastSTRUCTURE/plots/MG_K3.csv")
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

ggsave("MGK3_stru.pdf", plot = K3, width = 20, height = 3, device = "pdf")

###########
# K4 Plot
###########
hiroshige_palette <- met.brewer("Hiroshige", 8)
selected_colors <- hiroshige_palette[c(1, 4,6, 8)]

test<-read.csv("~/R/Cannabis_Genetics_Manuscript/fastSTRUCTURE/plots/MG_K4.csv")
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

ggsave("MGK4_stru.pdf", plot = K4, width = 20, height = 3, device = "pdf")

###########
# K5 Plot
###########
hiroshige_palette <- met.brewer("Hiroshige", 8)
selected_colors <- hiroshige_palette[c(1, 3,5,7, 8)]

test<-read.csv("~/R/Cannabis_Genetics_Manuscript/fastSTRUCTURE/plots/MG_K5.csv")
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

ggsave("MGK5_stru.pdf", plot = K5, width = 20, height = 3, device = "pdf")

#########################################
#FIGURE 3
#########################################
#Plot by use-type 

#####################
#Figure 3A - PCA LW
#####################
pca <- read.csv("~/R/Cannabis_Genetics_Manuscript/csv_files/lw_n498_typeknown_only.csv")

# Create the plot
p <- ggplot(data = pca, aes(x = EV1, y = EV2, group = Type, shape = Type, color = Type)) +
  geom_point(size = 4) +
  scale_shape_manual(values = seq(0, 15)) +
  scale_color_manual(values = c("turquoise", "blue", "red", "orange", "green", "grey")) +
  stat_ellipse(size = 1.5) +  # Adjust ellipse line thickness
  theme_classic() +
  xlab("PC 1 (5.52%)") + 
  ylab("PC 2 (3.57%)") +
  theme(
    axis.text = element_text(size = 25),  # Increase axis text size
    axis.title = element_text(size = 25),  # Increase axis title size
    legend.text = element_text(size = 25),  # Increase legend text size
    legend.title = element_text(size = 25)  # Increase legend title size
  )

# Display the plot
print(p)

# Save the plot to a PDF file
ggsave("LW_PCA.pdf", plot = p, width = 20, height = 13, device = "pdf")

#####################
#Figure 3B - PCA Phylos
#####################
pca <- read.csv("~/R/Cannabis_Genetics_Manuscript/csv_files/phylos_845_PRJNA347566_typeknown_pca_update.csv")

# Create the plot with custom legend titles
p <- ggplot(data = pca, aes(x = EV1, y = EV2, group = Type_final, shape = Type_final, color = Type_final)) +
  geom_point(size=4) +
  scale_shape_manual(values = seq(0, 15)) +
  scale_color_manual(values = c("turquoise", "blue", "red", "orange", "green")) +
  stat_ellipse(size = 1.5) +  # Adjust ellipse line thickness
  theme_classic() +
  xlab("PC 1 (6.66%)") + 
  ylab("PC 2 (5.34%)") +
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

# Display the plot
print(p)

# Save the plot to a PDF file
ggsave("Phylos845_PCA.pdf", plot = p, width = 20, height = 13, device = "pdf")

#####################
#Figure 3C - PCA Soorni
#####################

pca <- read.csv("~/R/Cannabis_Genetics_Manuscript/csv_files/soorni_n94_typeknown_pca.csv")
p <- ggplot(data = pca, aes(x = EV1, y = EV2, group = Type, shape = Type, color = Type)) +
  geom_point(size=4) +
  scale_shape_manual(values = seq(0, 15)) +
  scale_color_manual(values = c("red", "orange","green")) +
  stat_ellipse(size = 1.5) +  # Adjust ellipse line thickness
  theme_classic() +
  xlab("PC 1 (4.84%)") + 
  ylab("PC 2 (2.85%)") +
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

# Display the plot
print(p)

# Save the plot to a PDF file
ggsave("Soorni_PCA.pdf", plot = p, width = 20, height = 13, device = "pdf")


#####################
#Figure 3D - PCA MGV1
#####################

pca <- read.csv("~/R/Cannabis_Genetics_Manuscript/csv_files/MG_V1_n289_typeknown_pca.csv")
p <- ggplot(data = pca, aes(x = EV1, y = EV2, group = Type, shape = Type, color = Type)) +
  geom_point(size=4) +
  scale_shape_manual(values = seq(0, 15)) +
  scale_color_manual(values = c("red", "orange","green")) +
  stat_ellipse(size = 1.5) +  # Adjust ellipse line thickness
  theme_classic() +
  xlab("PC 1 (4.84%)") + 
  ylab("PC 2 (3.44%)") +
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

# Display the plot
print(p)

# Save the plot to a PDF file
ggsave("MGV1_PCA.pdf", plot = p, width = 20, height = 13, device = "pdf")


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

ggsave("phylos_bar_crom.pdf", plot = p1, width = 20, height = 13, device = "pdf")


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

ggsave("phylos_dot_chrom.pdf", plot = p2, width = 20, height = 13, device = "pdf")

#LeafWorks 
wild_all<-read.csv("~/R/Cannabis_Genetics_Manuscript/pi/both_leafworks_pi.csv")
#Plot all the datasets for all chroms side by side for diversity
p1 <- ggplot(wild_all, aes(x = Name, y = PI, fill = CHROM)) +
  geom_boxplot(outlier.size = 0) +
  theme_bw() + theme(legend.text=element_text(size=20)) +
  scale_fill_manual(values=met.brewer("Hokusai1", 10)) +
  theme(axis.text = element_text(size = 20))  + theme(text = element_text(size = 20))
p1

ggsave("LW_bar_crom.pdf", plot = p1, width = 20, height = 13, device = "pdf")

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
ggsave("LW_dot_chrom.pdf", plot = p2, width = 20, height = 13, device = "pdf")

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

#####################
#FIGURE 5B
#####################
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
f1<- f + scale_color_manual(values=met.brewer("Hokusai1", 4))

ggsave("LR_HCPC.pdf", plot = f1, width = 20, height = 13, device = "pdf")

#####################
#FIGURE 5C
#####################
#PCA (Regional Association) HK/LV/MB
pca <- read.csv("~/R/Cannabis_Genetics_Manuscript/csv_files/lw_nuclear_lr_subset_only.csv")
p <- ggplot(data = pca, aes(x=EV1, y=EV2,group=Ecotype,shape=Ecotype,color=Ecotype)) +
  geom_point(size=4) +
  scale_shape_manual(values = seq(0, 15)) +
  scale_color_manual(values=met.brewer("Hokusai1", 10)) +
  stat_ellipse(size = 1.5) +  # Adjust ellipse line thickness
  theme_classic() +
  xlab("PC 1 (9.94%)") + 
  ylab("PC 2 (8.55%)") +
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

# Display the plot
print(p)
ggsave("LR_PCA.pdf", plot = p, width = 20, height = 13, device = "pdf")

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

ggsave("LR_chrom_dot.pdf", plot = p2, width = 20, height = 13, device = "pdf")

#####################
#Figure 5E - fastSTRUCTURE LR LeafWorks
#####################
###########
# K3 Plot
###########
hiroshige_palette <- met.brewer("Hiroshige", 8)
selected_colors <- hiroshige_palette[c(1, 4, 6)]

test<-read.csv("~/R/Cannabis_Genetics_Manuscript/fastSTRUCTURE/plots/LR_K3.csv")
long_data <- pivot_longer(test, cols = c("K1", "K2", "K3"), names_to = "category", values_to = "value")

K3 <- ggplot(long_data, aes(x = Name, y = value, fill = category, group = Type)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = selected_colors) +  # Use selected colors from Hokusai palette
  labs(title = "", x = "", y = "Percent Identity") +
  facet_grid(~ Type, scales = "free_x") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 20),
    axis.text.y = element_text(size = rel(0.8)),
    plot.margin = margin(5.5, 5.5, 5.5, 5.5, "pt"),
    strip.background = element_blank(),
    legend.position = "right"
  )
K3

ggsave("LR_K3_stru.pdf", plot = K3, width = 20, height = 10, device = "pdf")

##################################################
# FIN
##################################################
