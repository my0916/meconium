# Load necessary libraries
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(glmnet) 
library(cowplot) 
library(ggpubr)  
library(RColorBrewer)

# Read clinical information and proteome data
clinical_info <- read_tsv("data/clinical_info.txt")
proteome_data <- read_csv("data/proteome_data.csv", check.names = FALSE)  

# Add exclusion criteria based on clinical conditions
clinical_info <- clinical_info %>% mutate(
  Excl = case_when(
    CHD == 1 | chromosomal_abnormality == 1 | disease_of_pediatric_surgery == 1 | Congenital_infection == 1 ~ 1,
    TRUE ~ 0
  )
)

# Combine clinical info with proteome data
GA_proteome <- bind_cols(clinical_info[, c("ID", "Excl", "GA", "Gender")], as_tibble(t(proteome_data[, clinical_info$ID])))
colnames(GA_proteome) <- c("ID", "Excl", "GA", "Gender", proteome_data$`Master accession Number`)

# Define color palette for plots
col1 <- brewer.pal(5, "Set2")

# Split data into training and validation sets
set.seed(123456)
GA_proteome_train <- GA_proteome %>% filter(Excl == 0) %>% sample_n(149)
GA_proteome_valid <- GA_proteome[!(GA_proteome$ID %in% GA_proteome_train$ID),]

# Train LASSO model on training data
set.seed(123456)
lasso.model.cv <- cv.glmnet(x = as.matrix(GA_proteome_train[, 5:ncol(GA_proteome_train)]), y = GA_proteome_train$GA, family = "gaussian", alpha = 1)
lasso.model <- glmnet(x = as.matrix(GA_proteome_train[, 5:ncol(GA_proteome_train)]), y = GA_proteome_train$GA, family = "gaussian", lambda = lasso.model.cv$lambda.min, alpha = 1)

# Predict GA for training data
est.GA_train <- predict(lasso.model, newx = as.matrix(GA_proteome_train[, 5:ncol(GA_proteome_train)]), s = lasso.model.cv$lambda.min, type = 'response')
est_actual_GA_train <- data.frame(
  ID = GA_proteome_train$ID,
  actual_GA = GA_proteome_train$GA,
  predicted_GA = as.vector(est.GA_train),
  Gender = factor(GA_proteome_train$Gender, levels = c("0", "1"))
)
est_actual_GA_train <- est_actual_GA_train %>% mutate(GA_delta = predicted_GA - actual_GA)

# Predict GA for validation data
est.GA_valid <- predict(lasso.model, newx = as.matrix(GA_proteome_valid[, 5:ncol(GA_proteome_valid)]), s = lasso.model.cv$lambda.min, type = 'response')
est_actual_GA_valid <- data.frame(
  ID = GA_proteome_valid$ID,
  actual_GA = GA_proteome_valid$GA,
  predicted_GA = as.vector(est.GA_valid),
  Gender = factor(GA_proteome_valid$Gender, levels = c("0", "1")),
  Excl = factor(GA_proteome_valid$Excl, levels = c("0", "1"))
)
est_actual_GA_valid <- est_actual_GA_valid %>% mutate(GA_delta = predicted_GA - actual_GA)

# Prepare training data for plotting
est_actual_GA_train <- est_actual_GA_train %>% mutate(Sex = case_when(Gender == "1" ~ "Female", TRUE ~ "Male"))
est_actual_GA_train$Sex <- factor(est_actual_GA_train$Sex, levels = c("Male", "Female"))

# Plot for training data
ggplot(est_actual_GA_train) +
  geom_abline(colour = "grey", linetype = 2) +
  stat_smooth(aes(x = actual_GA, y = predicted_GA), method = "lm", se = TRUE, linewidth = .5, fill = "blue", alpha = .1) +
  geom_point(aes(x = actual_GA, y = predicted_GA, fill = Sex), shape = 21, size = 2, alpha = 0.7) +
  scale_fill_manual(values = col1[c(3, 4)]) +
  theme_classic() +
  theme(aspect.ratio = 1, plot.title = element_text(hjust = .5, face = "bold"), legend.position = "none") +
  coord_cartesian(xlim = c(22.5, 41.5), ylim = c(22.5, 41.5)) +
  xlab("Gestational age (week)") + ylab("Predicted gestational age (week)") +
  ggtitle("Training cohort\n(Non-disease 149 cases)") +
  annotate("text", x = 26.2, y = 38, label = paste0("RMSE = ", round(rmse(est_actual_GA_train$actual_GA, est_actual_GA_train$predicted_GA), digits = 2))) +
  stat_cor(aes(x = actual_GA, y = predicted_GA), label.x = 23, label.y = 39.5)

# Prepare validation data for plotting
est_actual_GA_valid <- est_actual_GA_valid %>% mutate(Sex = case_when(Gender == "1" ~ "Female", TRUE ~ "Male"))
est_actual_GA_valid$Sex <- factor(est_actual_GA_valid$Sex, levels = c("Male", "Female"))

# Plot for non-disease validation cohort
est_actual_GA_valid.no <- est_actual_GA_valid %>% filter(Excl == "0")
ggplot(est_actual_GA_valid.no) +
  geom_abline(colour = "grey", linetype = 2) +
  stat_smooth(aes(x = actual_GA, y = predicted_GA), method = "lm", se = TRUE, linewidth = .5, fill = "blue", alpha = .1) +
  geom_point(aes(x = actual_GA, y = predicted_GA, fill = Sex), shape = 21, size = 2, alpha = 0.7) +
  scale_fill_manual(values = col1[c(3, 4)]) +
  theme_classic() +
  theme(aspect.ratio = 1, plot.title = element_text(hjust = .5, face = "bold"), legend.position = "none") +
  coord_cartesian(xlim = c(22.5, 41.5), ylim = c(22.5, 41.5)) +
  xlab("Gestational age (week)") + ylab("") +
  ggtitle("Validation cohort\n(Non-disease 55 cases)") +
  annotate("text", x = 26.2, y = 38, label = paste0("RMSE = ", round(rmse(est_actual_GA_valid.no$actual_GA, est_actual_GA_valid.no$predicted_GA), digits = 2))) +
  stat_cor(aes(x = actual_GA, y = predicted_GA), label.x = 23, label.y = 39.5)

# Plot for specific-disease validation cohort
est_actual_GA_valid.disease <- est_actual_GA_valid %>% filter(Excl == "1")
ggplot(est_actual_GA_valid.disease) +
  geom_abline(colour = "grey", linetype = 2) +
  stat_smooth(aes(x = actual_GA, y = predicted_GA), method = "lm", se = TRUE, linewidth = .5, fill = "blue", alpha = .1) +
  geom_point(aes(x = actual_GA, y = predicted_GA, fill = Sex), shape = 21, size = 2, alpha = 0.7) +
  scale_fill_manual(values = col1[c(3, 4)]) +
  theme_classic() +
  theme(aspect.ratio = 1, plot.title = element_text(hjust = .5, face = "bold"), legend.title = element_blank()) +
  coord_cartesian(xlim = c(22.5, 41.5), ylim = c(22.5, 41.5)) +
  xlab("Gestational age (week)") + ylab("") +
  ggtitle("Validation cohort\n(Specific-disease 55 cases)") +
  annotate("text", x = 26.2, y = 38, label = paste0("RMSE = ", round(rmse(est_actual_GA_valid.disease$actual_GA, est_actual_GA_valid.disease$predicted_GA), digits = 2)))+
  stat_cor(aes(x=actual_GA,y=predicted_GA), label.x=23, label.y=39.5)
