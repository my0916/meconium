# Load necessary libraries
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(dendextend)
library(RColorBrewer)
library(cowplot)

# Read clinical information and proteome data
clinical_info <- read_tsv("data/clinical_info.txt") 
proteome_data <- read_csv("data/proteome_data.csv", check.names = FALSE)

# Scale proteome data and prepare for merging with clinical info
proteome_data.t <- data.frame(scale(t(proteome_data[, 5:ncol(proteome_data)]))) 
colnames(proteome_data.t) <- proteome_data$`Master accession Number`
proteome_df <- proteome_data.t %>% rownames_to_column("ID")

# Merge proteome data with clinical information
proteome_clinical <- left_join(clinical_info[, c("ID", "GA", "Gender")], proteome_df, by = "ID") %>% arrange(GA)
proteome_clinical$Gender <- factor(proteome_clinical$Gender)

# Transform data into a long format for plotting
proteome_clinical.tf <- gather(proteome_clinical, 4:ncol(proteome_clinical), key = "Protein", value = "Expression")

# Plot GA trajectories for 3,433 proteins
ggplot(proteome_clinical.tf, aes(x = GA, y = Expression, group = Protein)) +
  theme_classic(base_size = 16) +
  theme(aspect.ratio = 1, plot.title = element_text(face = "bold", hjust = 0.5)) +
  stat_smooth(method = "loess", se = FALSE, linewidth = .01, color = "grey40") +
  ylab("Protein levels (z-score)") + xlab("Gestational age (week)") + 
  coord_cartesian(ylim = c(-2.51, 2.51)) +
  ggtitle("GA trajectories (3,433 proteins)") +
  scale_y_continuous(breaks = c(-3, -2, -1, 0, 1, 2, 3))

# Initialize res_summary to store loess fitted values
res_summary <- NULL

# Fit loess model for each protein and extract fitted values
for (i in 4:ncol(proteome_clinical)) {
  loess_model <- loess(proteome_clinical[, i, drop = TRUE] ~ GA, data = proteome_clinical)
  fitted_values <- data.frame(loess_model$fitted)
  colnames(fitted_values) <- colnames(proteome_clinical)[i]
  if (is.null(res_summary)) {
    res_summary <- fitted_values
  } else {
    res_summary <- cbind(res_summary, fitted_values)
  }
}

# Transpose the matrix and set column names
loess.mat <- t(res_summary)
colnames(loess.mat) <- proteome_clinical$ID

# Perform hierarchical clustering
dist <- dist(loess.mat, method = "euclidean")
hclust <- hclust(dist, method = "complete")

# Color branches based on clusters
col6 <- brewer.pal(6, "Set2")
col.cl6 <- col6[dendextend::cutree(tree = hclust, k = 6, order_clusters_as_data = FALSE)]

# Plot dendrogram
as.dendrogram(hclust) %>%
  dendextend::set("labels_colors", value = col.cl6) %>%
  dendextend::set("branches_k_color", k = 6, col6) %>%
  plot(leaflab = "none")
legend("topright",
       legend = paste0("Cluster ", 1:6),
       box.lty = 0, bg = "transparent", pch = 15, col = col6)

# Join clustering results with proteome data
clustering.res.df <- data.frame(Protein = rownames(loess.mat), Cluster = cutree(hclust, k = 6))
proteome_df.cluster <- left_join(proteome_clinical.tf, clustering.res.df, by = c("Protein" = "Protein"))

# Count proteins in each cluster
cluster.count <- proteome_df.cluster %>% dplyr::count(Cluster)

# Convert Cluster to factor for plotting
proteome_df.cluster$Cluster <- as.factor(proteome_df.cluster$Cluster)

# Prepare plots for each cluster
cluster_gg <- list()
for (cluster in 1:6) {
  cluster_data <- filter(proteome_df.cluster, Cluster == cluster)
  cluster_gg[[cluster]] <- ggplot(cluster_data, aes(x = GA, y = Expression)) +
    theme_classic(base_size = 16) +
    theme(aspect.ratio = 1, plot.title = element_text(face = "bold", hjust = 0.5), legend.position = "none") +
    geom_smooth(method = "loess", aes(group = Protein), se = FALSE, color = col6[cluster], linewidth = 0.01) + 
    geom_smooth(method = "loess", se = FALSE, color = "White", linewidth = 2) + 
    geom_smooth(method = "loess", se = FALSE, color = col6[cluster], linewidth = 1) +
    ylab("Protein levels (z-score)") + xlab("Gestational age (week)") + 
    coord_cartesian(ylim = c(-2.51, 2.51)) +
    scale_y_continuous(breaks = c(-3, -2, -1, 0, 1, 2, 3)) +
    ggtitle(paste0("Cluster ", cluster, " (", cluster.count$n[cluster], " proteins)"))
}

# Combine and plot all cluster plots in a grid
plot_grid(cluster_gg[[1]], cluster_gg[[2]], cluster_gg[[3]], cluster_gg[[4]], cluster_gg[[5]], cluster_gg[[6]], nrow = 2, align = "hv", axis = "tblr")