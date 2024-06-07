# Load necessary libraries
library(readr)
library(tibble)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(stats)

# Read clinical information and proteome data
clinical_info <- read_tsv("data/clinical_info.txt") 
proteome_data <- read_csv("data/proteome_data.csv", check.names = FALSE)

# Initialize an empty list to store results
res_summary <- list()

# Loop through each protein in the proteome data
for (i in 1:nrow(proteome_data)) {
  # Prepare temporary data for linear model
  tmp_data <- tibble(
    clinical_info[, c("ID", "GA", "Gender")], 
    value = t(proteome_data[i, clinical_info$ID])
  )
  tmp_data$value <- as.numeric(tmp_data$value)
  
  # Fit a linear model to assess the impact of GA and Gender
  res <- lm(value ~ GA + Gender, data = tmp_data)
  
  # Extract relevant coefficients and p-values
  res_summary[[i]] <- c(
    proteome_data[i, c(2, 3, 4)], 
    summary(res)$coefficients[2, c(1, 4)],  # GA coefficient and p-value
    summary(res)$coefficients[3, c(1, 4)]   # Gender coefficient and p-value
  )
}

# Convert the list of results to a data frame
res_summary <- do.call(rbind, res_summary)
res_summary <- data.frame(res_summary)
colnames(res_summary) <- c("Master_accession_Number", "Master_gene_symbol", "Protein_name", "GA_coef", "GA_pval", "Gender_coef", "Gender_pval")

# Convert columns to numeric
res_summary$GA_coef <- as.numeric(res_summary$GA_coef)
res_summary$GA_pval <- as.numeric(res_summary$GA_pval)
res_summary$Gender_coef <- as.numeric(res_summary$Gender_coef)
res_summary$Gender_pval <- as.numeric(res_summary$Gender_pval)

# Adjust p-values using the Benjamini-Hochberg method
res_summary$GA_qval <- p.adjust(res_summary$GA_pval, method = "BH")
res_summary$Gender_qval <- p.adjust(res_summary$Gender_pval, method = "BH")

# Prepare data for Gender effect plot
sex_genes <- res_summary$Master_gene_symbol
sex_selected_genes <- res_summary[order(res_summary$Gender_pval), "Master_gene_symbol"][1:10]
sex_genes[!sex_genes %in% sex_selected_genes] <- ""

# Plot Gender effect
ggplot(data = res_summary) +
  geom_point(aes(x = Gender_coef, y = -log10(Gender_pval), 
                 fill = -log10(Gender_pval) * sign(Gender_coef), 
                 size = (-log10(Gender_pval))^2), 
             shape = 21, color = "black", stroke = 0.3) +
  geom_text_repel(aes(x = Gender_coef, y = -log10(Gender_pval), label = sex_genes), size = 3, max.overlaps = 99) +
  scale_size(range = c(0.1, 4)) +
  scale_fill_gradientn(colors = c("blue", "white", "red"), breaks = c(-20, 0, 20), limits = c(-20, 20), oob = scales::squish) +
  theme_classic(base_size = 16) +
  theme(aspect.ratio = 1, legend.position = "none") +
  ylab(bquote(~-log[10] ~ italic(P) ~ "-value")) +
  xlab("Coefficient") +
  ggtitle("Gender Effect on Protein Expression")

# Prepare data for GA effect plot
GA_genes <- res_summary$Master_gene_symbol
GA_selected_genes <- res_summary[order(res_summary$GA_pval), "Master_gene_symbol"][1:10]
GA_genes[!GA_genes %in% GA_selected_genes] <- ""

# Plot GA effect
ggplot(data = res_summary) +
  geom_point(aes(x = GA_coef, y = -log10(GA_pval), 
                 fill = -log10(GA_pval) * sign(GA_coef), 
                 size = (-log10(GA_pval))^2), 
             shape = 21, color = "black", stroke = 0.3) +
  geom_text_repel(aes(x = GA_coef, y = -log10(GA_pval), label = GA_genes), size = 3, max.overlaps = 99) +
  scale_size(range = c(0.1, 4)) +
  scale_fill_gradientn(colors = c("blue", "white", "red"), breaks = c(-20, 0, 20), limits = c(-20, 20), oob = scales::squish) +
  theme_classic(base_size = 16) +
  theme(aspect.ratio = 1, legend.position = "none") +
  ylab(bquote(~-log[10] ~ italic(P) ~ "-value")) +
  xlab("Coefficient") +
  ggtitle("Gestational Age Effect on Protein Expression")