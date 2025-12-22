library(dplyr)
library(randomForest)
library(themis)
library(caret)
library(pROC)
library(ggplot2)
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
names(train_response) <- rownames(train_genus)

test_genus <- readRDS("data/test_genus.rds")
test_family <- readRDS("data/test_family.rds")
test_response <- readRDS("data/test_response.rds")
colnames(test_genus) <- paste0("g_", colnames(test_genus), "")
colnames(test_family) <- paste0("f_", colnames(test_family), "")
names(test_response) <- rownames(test_genus)

train_feature <- cbind(train_genus)
test_feature <- cbind(test_genus)
all_feature <- rbind(train_feature)

ms_obj <- newMRexperiment(t(train_feature))
ms_obj <- cumNorm(ms_obj, p = 0.5)
train_feature <- t(MRcounts(ms_obj, norm = TRUE, log = TRUE))

ms_obj <- newMRexperiment(t(test_feature))
ms_obj <- cumNorm(ms_obj, p = 0.5)
test_feature <- t(MRcounts(ms_obj, norm = TRUE, log = TRUE))
}


{
    species_summed <- readRDS("adata2025/data/species_summed.rds")
    genus_summed <- readRDS("adata2025/data/genus_summed.rds")
    family_summed <- readRDS("adata2025/data/family_summed.rds")

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


my_rf_train <- function(tr_res, te_res, sig_feature) {
  train_feature <- tr_res$feature
  train_response <- tr_res$response
  test_feature <- te_res$feature
  test_response <- te_res$response
  sig_train_feature <- train_feature[, sig_feature]
  sig_test_feature <- test_feature[, sig_feature]
  sig_train_feature <- sig_train_feature / rowSums(sig_train_feature)
  sig_test_feature <- sig_test_feature / rowSums(sig_test_feature)
  sig_train_feature[is.na(sig_train_feature)] <- 0
  sig_test_feature[is.na(sig_test_feature)] <- 0
  train_control <- trainControl(method = "cv", number = 5, sampling = "up", classProbs = TRUE, summaryFunction = twoClassSummary)
  
  rf_model <- train(sig_train_feature, train_response, method = "rf", trControl = train_control, metric = "ROC")

  mean_auc_cv <- getTrainPerf(rf_model)$TrainROC

  rf_train_pred <- predict(rf_model, newdata = sig_train_feature, type = "prob")
  rf_train_pred_class <- predict(rf_model, newdata = sig_train_feature)
  confusion_matrix_train <- confusionMatrix(rf_train_pred_class, train_response)

  roc_train <- roc(train_response, rf_train_pred[, 2], levels = levels(train_response))
  roc_train_auc <- as.numeric(auc(roc_train))
  #roc_train_auc <- mean_auc_cv


  rf_test_pred <- predict(rf_model, newdata = sig_test_feature, type = "prob")
  rf_test_pred_class <- predict(rf_model, newdata = sig_test_feature)
  confusion_matrix_test <- confusionMatrix(rf_test_pred_class, test_response)
  # ROC AUC
  roc_test <- roc(test_response, rf_test_pred[, 2], levels = levels(test_response))
  roc_test_auc <- as.numeric(auc(roc_test))
  return(list(tr_auc = roc_train_auc, tr_mt = confusion_matrix_train,
                te_auc = roc_test_auc, te_mt = confusion_matrix_test))
}


glm_feaure <- function(tr_res, te_res, sig_feature) {
  train_feature <- tr_res$feature
  train_response <- tr_res$response
  test_feature <- te_res$feature
  test_response <- te_res$response
  up_data_full <- upSample(x = train_feature, y = train_response)
  X_train_full <- up_data_full[, -ncol(up_data_full)]
  y_train_full <- up_data_full$Class
  feature <- X_train_full[, sig_feature]
  train_df <- data.frame(feature)
  colnames(train_df) <- sig_feature
  train_df$response <- y_train_full
  model <- glm(response ~., data = train_df, family = binomial())
  train_predictions <- predict(model, type = "response")

  train_roc <- roc(y_train_full, train_predictions)

  train_auc <- auc(train_roc)
  test_df <- data.frame(test_feature[, sig_feature])
  colnames(test_df) <- sig_feature
  test_predictions <- predict(model,
   newdata = test_df, type = "response")

  test_roc <- roc(test_response, test_predictions)
  test_auc <- auc(test_roc)
  return(list(tr_auc = train_auc, te_auc = test_auc))
}


compares_list <- list(c("Normal", "Polyp"), c("Polyp", "CRC"), c("Normal", "CRC"))

compare_num <- 1
compares <- compares_list[[compare_num]]
sub_name <- "MSS"
{
    patient_info <- readRDS("adata2025/data/info_tb.rds")
    crc_id <- rownames(train_feature)[train_response == "CRC"]
    info_tb <- patient_info[, crc_id]
    locations <- get_group_label(info_tb, "MSI")
    table(locations)
    ids <- names(locations)[which(locations == "MSS")]
    control_ids <- rownames(train_feature)[train_response != "CRC"]
    train_feature <- train_feature[c(ids, control_ids), ]
    train_response <- train_response[c(ids, control_ids)]
    tr_res <- deal_compares_tr(compares, train_feature, train_response)
    table(tr_res$response)

    crc_id <- rownames(test_feature)[test_response == "CRC"]
    info_tb <- patient_info[, crc_id]
    locations <- get_group_label(info_tb, "MSI")
    table(locations)
    ids <- names(locations)[which(locations == "MSS")]
    control_ids <- rownames(test_feature)[test_response != "CRC"]
    test_feature <- test_feature[c(ids, control_ids), ]
    test_response <- test_response[c(ids, control_ids)]
    te_res <- deal_compares_tr(compares, test_feature, test_response)
    table(te_res$response)

    feature1 <- colnames(genus_summed)
    cleaned_string <- sub("_unclassified$", "", feature1)
    feature1 <- feature1[which(feature1 %in% cleaned_string)]

    feature_summed <- rbind(tr_res$feature)
    response <- c(tr_res$response)

    feature_summed <- feature_summed[, feature1]
}
set.seed(100)
model <- randomForest(x = feature_summed, y = response, importance = TRUE)

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

pdf(paste0("MicrobioInCRC/figs/figS10/", sub_name, "_", compare_num, "_importance.pdf"), width = 5, height = 4)
print(p1)
dev.off()

sig_feature <- importance$Feature

sig_data <- as.data.frame(feature_summed[, sig_feature])

sig_data$Class <- response

non_zero_proportions <- sapply(sig_feature, function(feature) {
  tapply(sig_data[[feature]] != 0, sig_data$Class, mean)
})
library(ggplot2)
library(reshape2)

non_zero_proportions_melted <- melt(non_zero_proportions)
colnames(non_zero_proportions_melted) <- c("Class", "Feature", "NonZeroProportion")

non_zero_proportions_melted$Feature <- factor(non_zero_proportions_melted$Feature, levels = rev(sig_feature))

pdf(paste0("MicrobioInCRC/figs/figS10/", sub_name, "_", compare_num, "_non_zeros.pdf"), width = 5, height = 4)
ggplot(non_zero_proportions_melted, aes(x = Feature, y = NonZeroProportion, fill = Class)) +
  geom_bar(stat = "identity", position = "dodge") +
  coord_flip() +
  theme_minimal() +
  labs(x = "Feature",
       y = "Non-zero Proportion") +
  scale_fill_brewer(palette = "Set2")
dev.off()

my_heatmap <- function(train_feature, sig_feature, level, name) {

  sig_matrix <- train_feature[, sig_feature]
  cor_matrix <- cor(sig_matrix)
  melted_cormat <- melt(cor_matrix)
  melted_cormat$Var1 <- factor(melted_cormat$Var1, levels = level)
  melted_cormat$Var2 <- factor(melted_cormat$Var2, levels = level)
  pdf(name, width = 5, height = 6)
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


name <- paste0("MicrobioInCRC/figs/figS10/", sub_name, "_", compare_num, "_heatmap.pdf")

my_heatmap(train_feature, sig_feature, sig_feature, name)

selected_feature <- sig_feature[c(1, 2, 3, 4, 5, 6, 8, 9, 10)] #Left compare2
selected_feature <- sig_feature[c(1, 2, 3, 5, 6, 7, 8, 9, 10)] #Left compare3
selected_feature <- sig_feature[c(1, 2, 3, 4, 5, 6, 7, 8, 10)] #Right compare2
selected_feature <- sig_feature[c(1, 2, 3, 5, 6, 7, 12)] #Right compare3
selected_feature <- sig_feature[c(1, 2, 3, 4, 5, 6, 7, 9)] #MSS compare2
selected_feature <- sig_feature[c(1, 2, 3, 4, 5, 6, 8, 9)] #MSS compare3

auc_df <- data.frame(AUC = numeric(), feature = character(), type = character())
j <- 1
non_zeros <- c()
set.seed(123)
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

library(ggsci)

pdf(paste0("MicrobioInCRC/figs/figS10/", sub_name, "_", compare_num, "_auc.pdf"), width = 4, height = 3)
ggplot(auc_df, aes(x = feature, y = AUC, color = type, group = type)) +
  geom_point(size = 2.5) +
  geom_line() +
  theme_minimal() +
  scale_color_d3() +
  theme(axis.text.x = element_text(angle = -45),
  axis.line = element_line(linewidth = 0.8))
dev.off()

left_2 <- selected_feature[1:7]
left_3 <- selected_feature[1:7]
right_2 <- selected_feature[1:5]
right_3 <- selected_feature[1:6]
mss_2 <- selected_feature[1:5]
mss_3 <- selected_feature[1:6]
all_2 <- c("Parvimonas", "Collinsella", "Pseudomonas", "Fusobacterium")
all_3 <- c("Fusobacterium", "Pseudomonas", "Collinsella", "Phascolarctobacterium", "Granulicatella", "Leptotrichia", "Streptococcus")
all_2 <- paste0("g_", all_2)
all_3 <- paste0("g_", all_3)
library(VennDiagram)
    x <- NULL
    x[["Left"]] <- left_3
    x[["Right"]] <- right_3
    x[["MSS"]] <- mss_3
    x[["All"]] <- all_3
venn.plot <- venn.diagram(
    x = x,
    filename = NULL,
    category.names = c("Left", "Right", "MSS", "All"),
    output = TRUE,
    height = 200, 
    width = 200, 
    resolution = 300,
    compression = "lzw",
    lwd = 1.5,
    col = "black",
    fill = c("#a6cee3", "#b2df8a", "#fb9a99", "#fdbf6f"),
    alpha = 0.50,
    cex = 1.2,
    fontface = "bold",
    fontfamily = "sans",
    cat.cex = 1.5,
    cat.fontface = "bold",
    cat.default.pos = "outer",
    cat.pos = c(-27, 27, 90, 90),
    cat.dist = c(0.06, 0.06, 0.06, 0.06),
    cat.just = list(c(0.5, 0.5), c(0.5, 0.5), c(0.5, 0.5), c(0.5, 0.5))
)
pdf(paste0("MicrobioInCRC/figs/figS10/veen_normal_crc.pdf"), height = 4, width = 4)
    grid.draw(venn.plot)
dev.off()