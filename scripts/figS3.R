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
library(gt)
library(patchwork)
library(tidyr)
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

train_feature <- cbind(train_genus, train_family)
test_feature <- cbind(test_genus, test_family)

ms_obj <- newMRexperiment(t(train_feature))
ms_obj <- cumNorm(ms_obj, p = 0.5)
train_feature <- t(MRcounts(ms_obj, norm = TRUE, log = FALSE))

ms_obj <- newMRexperiment(t(test_feature))
ms_obj <- cumNorm(ms_obj, p = 0.5)
test_feature <- t(MRcounts(ms_obj, norm = TRUE, log = FALSE))
}

cv_rf_three <- function(train_feature, train_response, test_feature, test_response) {

  class_levels <- levels(train_response)
  cv_metrics <- array(0, dim = c(5, 3, 3),
                      dimnames = list(
                        Fold = 1:5,
                        Class = class_levels,
                        Metric = c("Precision", "Recall", "F1")
                      ))
  

  set.seed(123)
  folds <- createFolds(train_response, k = 5)
  

  for (i in 1:5) {
    train_idx <- unlist(folds[-i])
    val_idx <- folds[[i]]
    
    up_data <- upSample(x = train_feature[train_idx, ], 
                        y = train_response[train_idx])
    X_train_up <- up_data[, -ncol(up_data)]
    y_train_up <- up_data$Class
    
    model <- randomForest(x = X_train_up, y = y_train_up)
    
    pred_val <- predict(model, train_feature[val_idx, ])
    
    cm <- confusionMatrix(pred_val, train_response[val_idx])
    
    for (cls in class_levels) {
      cls_idx <- which(class_levels == cls)
      precision <- cm$byClass[cls_idx, "Pos Pred Value"]
      recall <- cm$byClass[cls_idx, "Sensitivity"]
      f1 <- 2 * (precision * recall) / (precision + recall)
      
      cv_metrics[i, cls, "Precision"] <- precision
      cv_metrics[i, cls, "Recall"] <- recall
      cv_metrics[i, cls, "F1"] <- ifelse(is.nan(f1), 0, f1)
    }
  }
  

  cv_mean_metrics <- apply(cv_metrics, c(2, 3), mean, na.rm = TRUE)
  
  up_data_full <- upSample(x = train_feature, y = train_response)
  X_train_full <- up_data_full[, -ncol(up_data_full)]
  y_train_full <- up_data_full$Class
  
  final_model <- randomForest(x = X_train_full, y = y_train_full)
  

  pred_test <- predict(final_model, test_feature)
  cm_test <- confusionMatrix(pred_test, test_response)
  

  test_metrics <- matrix(0, nrow = 3, ncol = 3,
                        dimnames = list(
                          Class = class_levels,
                          Metric = c("Precision", "Recall", "F1")
                        ))
  
  for (cls in class_levels) {
    cls_idx <- which(class_levels == cls)
    precision <- cm_test$byClass[cls_idx, "Pos Pred Value"]
    recall <- cm_test$byClass[cls_idx, "Sensitivity"]
    f1 <- 2 * (precision * recall) / (precision + recall)
    
    test_metrics[cls, ] <- c(precision, recall, ifelse(is.nan(f1), 0, f1))
  }
  
  return(list(
    cross_validation_metrics = cv_mean_metrics,
    test_set_metrics = test_metrics,
    test_confusion_matrix = cm_test$table
  ))
}

results <- cv_rf_three(train_feature, train_response,
                         test_feature, test_response)

results$cross_validation_metrics
results$test_set_metrics

cm_data <- as.data.frame(results$test_confusion_matrix)
  colnames(cm_data) <- c("Predicted", "Actual", "Count")

pdf("figs/figS3/cm_1.pdf", height = 3, width = 3.2)
  heatmap_plot <- ggplot(cm_data, aes(x = Actual, y = Predicted, fill = Count)) +
    geom_tile(color = "white", linewidth = 0.8) +
    geom_text(aes(label = Count), color = "black", size = 5) +
    scale_fill_gradient(low = "#F7FBFF", high = "#2171B5") +
    labs(x = "Actual Class", y = "Predicted Class",
         title = "Confusion Matrix Heatmap") +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      panel.grid = element_blank(),
      axis.text = element_text(color = "black"),
      legend.position = "right"
    ) +
    coord_fixed()
print(heatmap_plot)
dev.off()

tr_res <- deal_compares_three(train_feature, train_response)
te_res <- deal_compares_three(test_feature, test_response)
results <- cv_rf_three(tr_res$feature, tr_res$response, te_res$feature, te_res$response)

results$test_confusion_matrix

cm_data <- as.data.frame(results$test_confusion_matrix)
  colnames(cm_data) <- c("Predicted", "Actual", "Count")
  
pdf("figs/figS3/cm_2.pdf", height = 3, width = 3.2)
  heatmap_plot <- ggplot(cm_data, aes(x = Actual, y = Predicted, fill = Count)) +
    geom_tile(color = "white", linewidth = 0.8) +
    geom_text(aes(label = Count), color = "black", size = 5) +
    scale_fill_gradient(low = "#F7FBFF", high = "#2171B5") +
    labs(x = "Actual Class", y = "Predicted Class",
         title = "Confusion Matrix Heatmap") +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      panel.grid = element_blank(),
      axis.text = element_text(color = "black"),
      legend.position = "right"
    ) +
    coord_fixed()
print(heatmap_plot)
dev.off()


visualize_features <- function(model, train_feature, train_response) {
  library(ggplot2)
  library(ggpubr)
  library(tidyr)
  library(dplyr)
  
  train_df <- as.data.frame(train_feature)
  if (is.null(colnames(train_df))) {
    colnames(train_df) <- paste0("V", 1:ncol(train_df))
  }
  
  importance <- importance(model) %>% 
    as.data.frame() %>%
    tibble::rownames_to_column("Feature") %>%
    arrange(desc(MeanDecreaseGini)) %>%
    head(15)
  colnames(importance)[2] <- "Importance"

  existing_features <- intersect(importance$Feature, colnames(train_df))
  if (length(existing_features) == 0) {
    stop("No matching features between importance matrix and training data")
  }
  
  p1 <- ggplot(importance, aes(x = reorder(Feature, Importance), y = Importance)) +
    geom_col(fill = "#4292C6", width = 0.7) +
    coord_flip() +
    labs(x = "Feature", y = "Mean Decrease in Gini",
         title = "Top 15 Feature Importance") +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid.major.y = element_blank(),
      axis.text.y = element_text(color = "black", size = 10),
      plot.title = element_text(hjust = 0.5, face = "bold")
    )
  
    feature_means <- train_df %>%
    select(all_of(existing_features)) %>%
    mutate(Class = train_response) %>%
    group_by(Class) %>%
    summarise(across(everything(), mean)) %>%
    pivot_longer(-Class, names_to = "Feature", values_to = "Mean") %>%
    pivot_wider(names_from = Class, values_from = Mean) %>%
    tibble::column_to_rownames("Feature")
  
  scaled_means <- as.data.frame(t(scale(t(feature_means))))

  df2 <- scaled_means %>%
    tibble::rownames_to_column("Feature") %>%
    pivot_longer(-Feature, names_to = "Class", values_to = "Value")
  df2$Class <- factor(df2$Class, levels = c("Normal", "Polyp", "CRC"))
  
  p2 <- df2 %>%
    ggplot(aes(x = Class, y = Feature, fill = Value)) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(aes(label = format(Value, digits = 2)), color = "black", size = 3) +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                        midpoint = 0, limits = c(-3, 3)) +
    labs(x = "Class", y = "Feature",
         title = "Standardized Feature Means by Class") +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
      axis.text.y = element_text(color = "black"),
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.position = "right"
    ) +
    coord_fixed(ratio = 0.6)
  
  ggarrange(p1, p2, ncol = 2, labels = c("A", "B"), 
           font.label = list(size = 14))
}

up_data_full <- upSample(x = train_feature, y = train_response)
  X_train_full <- up_data_full[, -ncol(up_data_full)]
  y_train_full <- up_data_full$Class
  
final_model <- randomForest(x = X_train_full, y = y_train_full)

model <- final_model


pdf("figs/figS3/feature_improtance.pdf", height = 6, width = 6)
print(p1)
dev.off()

pdf("figs/figS3/feature_heatmap.pdf", height = 6, width = 6)
print(p2)
dev.off()
