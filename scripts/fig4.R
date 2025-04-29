library(dplyr)
library(randomForest)
library(themis)
library(caret)
library(pROC)
library(ggplot2)
library(ggsci)
library(MLmetrics)
library(reshape2)
library(ROSE)
library(doParallel)
library(metagenomeSeq)
source("scripts/methods.R")

#pre-process
{
train_genus <- readRDS("data/train_genus.rds")
train_family <- readRDS("data/train_family.rds")
train_response <- readRDS("data/train_response.rds")
colnames(train_genus) <- paste0("g_", colnames(train_genus), "")
colnames(train_family) <- paste0("f_", colnames(train_family), "")

test_genus <- readRDS("data/test_genus.rds")
test_family <- readRDS("data/test_family.rds")
test_response <- readRDS("data/test_response.rds")
colnames(test_genus) <- paste0("g_", colnames(test_genus), "")
colnames(test_family) <- paste0("f_", colnames(test_family), "")

train_feature <- cbind(train_genus)
test_feature <- cbind(test_genus)
all_feature <- rbind(train_feature)

ms_obj <- newMRexperiment(t(train_feature))
ms_obj <- cumNorm(ms_obj, p = 0.5)
train_feature <- t(MRcounts(ms_obj, norm = TRUE, log = FALSE))

ms_obj <- newMRexperiment(t(test_feature))
ms_obj <- cumNorm(ms_obj, p = 0.5)
test_feature <- t(MRcounts(ms_obj, norm = TRUE, log = FALSE))
}


{
    species_summed <- readRDS("data/species_summed.rds")
    genus_summed <- readRDS("data/genus_summed.rds")
    family_summed <- readRDS("data/family_summed.rds")

    species_summed <- species_summed[which(substr(rownames(species_summed), 1, 1) %in% c("N", "a", "T")), ]
    genus_summed <- genus_summed[which(substr(rownames(genus_summed), 1, 1) %in% c("N", "a", "T")), ]
    family_summed <- family_summed[which(substr(rownames(family_summed), 1, 1) %in% c("N", "a", "T")), ]

    species_summed <- species_summed[, which(colSums(species_summed > 0) > 0.05*nrow(species_summed))]
    genus_summed <- genus_summed[, which(colSums(genus_summed > 0) > 0.05*nrow(genus_summed))]
    family_summed <- family_summed[, which(colSums(family_summed > 0) > 0.05*nrow(family_summed))]
    response <- substr(rownames(family_summed), 1, 1)
    response = case_when(
        response == "N" ~ "Normal",
        response == "a" ~ "Polyp",
        response == "T" ~ "CRC"
    )
    response <- factor(response, levels = c("Normal", "Polyp", "CRC"))
    colnames(species_summed) <- paste0("s_", colnames(species_summed), "")
    colnames(genus_summed) <- paste0("g_", colnames(genus_summed), "")
    colnames(family_summed) <- paste0("f_", colnames(family_summed), "")
}


compares_list <- list(c("Normal", "Polyp"), c("Normal", "Polyp_CRC_1"),
 c("Polyp", "CRC"), c("Polyp", "CRC_1"),
 c("Normal_Polyp", "CRC"), c("Normal_Polyp", "CRC_1"), c("Normal", "CRC"), c("Normal", "CRC_1"))

compare_num <- 4
compares <- compares_list[[compare_num]]

tr_res <- deal_compares_tr(compares, train_feature, train_response)
table(tr_res$response)
te_res <- deal_compares_tr(compares, test_feature, test_response)
table(te_res$response)

compare_res <- deal_compares_tr(compares, genus_summed, response)

feature_summed <- compare_res$feature
response <- compare_res$response

set.seed(100)
model <- randomForest(x = compare_res$feature, y = compare_res$response, importance = TRUE)

importance <- importance(model) %>% 
    as.data.frame() %>%
    tibble::rownames_to_column("Feature") %>%
    arrange(desc(MeanDecreaseAccuracy)) %>%
    head(15)
  colnames(importance)[4] <- "Importance"

  p1 <- ggplot(importance, aes(x = reorder(Feature, Importance), y = Importance)) +
    geom_col(fill = "#4292C6", width = 0.7) +
    coord_flip() +
    labs(x = "Feature", y = "Mean Decrease in Accuracy",
         title = "Top 15 Feature Importance") +
    theme_minimal() +
    theme(
      panel.grid.major.y = element_blank(),
      axis.text.y = element_text(color = "black", size = 10),
      plot.title = element_text(hjust = 0.5, face = "bold")
    )
pdf("figs/fig4/com4_importance.pdf", width = 5, height = 4)
print(p1)
dev.off()

sig_feature <- importance$Feature

sig_data <- as.data.frame(feature_summed[, sig_feature])

sig_data$Class <- response

non_zero_proportions <- sapply(sig_feature, function(feature) {
  tapply(sig_data[[feature]] != 0, sig_data$Class, mean)
})


non_zero_proportions_melted <- melt(non_zero_proportions)
colnames(non_zero_proportions_melted) <- c("Class", "Feature", "NonZeroProportion")

non_zero_proportions_melted$Feature <- factor(non_zero_proportions_melted$Feature, levels = rev(sig_feature))


pdf("figs/fig4/com4_non_zeros.pdf", width = 5, height = 4)
ggplot(non_zero_proportions_melted, aes(x = Feature, y = NonZeroProportion, fill = Class)) +
  geom_bar(stat = "identity", position = "dodge") +
  coord_flip() +
  theme_minimal() +
  labs(title = "Non-zero Proportion for Top 15 Features by Class",
       x = "Feature",
       y = "Non-zero Proportion") +
  scale_fill_brewer(palette = "Set2")
dev.off()

my_heatmap <- function(train_feature, sig_feature, level, name) {

  sig_matrix <- train_feature[, sig_feature]
  cor_matrix <- cor(sig_matrix)

  melted_cormat <- melt(cor_matrix)
  melted_cormat$Var1 <- factor(melted_cormat$Var1, levels = level)
  melted_cormat$Var2 <- factor(melted_cormat$Var2, levels = level)

  pdf(paste0("figs/fig4/", name, ".pdf"), width = 5, height = 6)
  p <- ggplot(melted_cormat, aes(x = Var1, y = Var2, fill = value)) +
    geom_tile(color = "white") +
    scale_fill_gradient2(
      low = "blue",          
      high = "red",        
      mid = "white",     
      midpoint = 0,          
      limit = c(-1, 1)
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = -90),
      legend.position = "top" 
    ) +
    labs(
      x = "",
      y = "",
      fill = "Correlation"
    ) +
    coord_fixed()
  print(p)
  dev.off()
}


my_heatmap(train_feature, sig_feature, sig_feature, "com4_heatmap")

selected_feature <- sig_feature[c(2, 3, 4, 5, 7, 9)] #compare 4

selected_feature <- sig_feature[c(2, 4, 5, 6, 9, 12, 13)] #compare 1

selected_feature <- sig_feature[c(2, 5, 6, 7, 8, 10, 12, 13)] #compare 7

selected_feature <- sig_feature[c(2, 4, 5, 9, 13, 14)] #compare 3

auc_df <- data.frame(AUC = numeric(), feature = character(), type = character())
j <- 1
non_zeros <- c()
for (i in 1:length(selected_feature)) {
  if (i <= 1) {
    sig_feature2 <- selected_feature[1:i]
    rf_res <- glm_feaure(tr_res, te_res, sig_feature2)
    auc_df[j, ] <- c(rf_res$tr_auc, selected_feature[i], "Train")
    j <- j + 1
    auc_df[j, ] <- c(rf_res$te_auc, selected_feature[i], "Test")
    j <- j + 1
    if(i == 1) {
      non_zeros <- c(non_zeros, sum(all_feature[, sig_feature2] != 0))
    }else {
      non_zeros <- c(non_zeros, sum(rowSums(all_feature[, sig_feature2]) != 0))
    }
  }else {
    sig_feature2 <- selected_feature[1:i]
    rf_res <- my_rf_train(tr_res, te_res, sig_feature2)
    auc_df[j, ] <- c(rf_res$tr_auc, selected_feature[i], "Train")
    j <- j + 1
    auc_df[j, ] <- c(rf_res$te_auc, selected_feature[i], "Test")
    j <- j + 1
    non_zeros <- c(non_zeros, sum(rowSums(all_feature[, sig_feature2]) != 0))
  }
}
non_zeros <- non_zeros / nrow(all_feature)
auc_df$AUC <- as.numeric(auc_df$AUC)
auc_df$feature <- factor(auc_df$feature, levels = selected_feature)
auc_df$type <- factor(auc_df$type, levels = c("Train", "Test"))

pdf("figs/fig4/Com_4_auc.pdf", width = 4, height = 3)
ggplot(auc_df, aes(x = feature, y = AUC, color = type, group = type)) +
  geom_point(size = 2.5) +
  geom_line() +
  theme_minimal() +
  scale_color_d3() +
  theme(axis.text.x = element_text(angle = -45),
  axis.line = element_line(linewidth = 0.8))
dev.off()


####ROC

roc_problems <- data.frame()
name_ps <- c()
for (compare_num in c(1, 3, 7, 4)) {
    compares <- compares_list[[compare_num]]
    deal_res <- deal_compares_tr(compares, train_feature, train_response)
    response <- deal_res$response
    predictors <- deal_res$feature
    res_roc <- cv_roc_mean(response, predictors)
    roc_data <- res_roc$roc_data
    auc <- round(res_roc$auc, 3)
    problem <- paste(compares[1], "vs", compares[2])
    roc_data$problem <- rep(paste0(problem, "(AUC = ", auc, ")"), dim(roc_data)[1])
    name_ps <- c(name_ps, paste0(problem, "(AUC = ", auc, ")"))
    roc_problems <- rbind(roc_problems, roc_data)
}

roc_problems$problem <- factor(roc_problems$problem, levels = name_ps)

pdf("figs/fig4/Normal_Polyp_CRC.pdf", height = 2.4, width = 4.8)
ggplot(roc_problems, aes(x = 1 - fpr, y = tpr, color = problem)) +
  geom_line(size = 1) +
  geom_abline(linetype = "dashed") +
  labs(x = "1 - Specificity", y = "Sensitivity") +
  theme_classic()
dev.off()

roc_problems <- data.frame()
name_ps <- c()
for (compare_num in c(1, 3, 7, 4)) {
    compares <- compares_list[[compare_num]]
    deal_res <- deal_compares_tr(compares, train_feature, train_response)
    response <- deal_res$response
    predictors <- deal_res$feature
    test_res <- deal_compares_tr(compares, test_feature, test_response)
    res_roc <- test_roc(predictors, response, test_res$feature, test_res$response)
    roc_data <- res_roc$roc_data
    auc <- round(res_roc$auc, 3)
    problem <- paste(compares[1], "vs", compares[2])
    roc_data$problem <- rep(paste0(problem, "(AUC = ", auc, ")"), dim(roc_data)[1])
    name_ps <- c(name_ps, paste0(problem, "(AUC = ", auc, ")"))
    roc_problems <- rbind(roc_problems, roc_data)
}

roc_problems$problem <- factor(roc_problems$problem, levels = name_ps)

pdf("figs/fig4/Normal_Polyp_CRC_test.pdf", height = 2.4, width = 4.8)
ggplot(roc_problems, aes(x = 1 - fpr, y = tpr, color = problem)) +
  geom_line(size = 1) +
  geom_abline(linetype = "dashed") +
  labs(x = "1 - Specificity", y = "Sensitivity") +
  theme_classic()
dev.off()