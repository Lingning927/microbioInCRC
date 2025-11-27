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


non_zeros <- c()
for (i in 1:length(selected_feature)) {
    sig_feature2 <- selected_feature[1:i]
    if(i == 1) {
      non_zeros <- c(non_zeros, sum(all_feature[, sig_feature2] != 0))
    }else {
      non_zeros <- c(non_zeros, sum(rowSums(all_feature[, sig_feature2]) != 0))
    }
}
non_zeros <- non_zeros / nrow(all_feature)


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

get_cv_auc <- function(genus_summed, response, k = 5) {
  set.seed(100)
  response <- as.factor(response)
  if (length(levels(response)) != 2) stop("response必须为二分类因子")
  
  train_control <- trainControl(
    method = "cv",
    number = 5,
    sampling = "up",
    classProbs = TRUE,
    summaryFunction = twoClassSummary,
    savePredictions = "final"
  )
  
  response <- relevel(response, ref = levels(response)[2])
  
  rf_model <- train(
    x = genus_summed,
    y = response,
    method = "rf",
    metric = "ROC",
    trControl = train_control
  )
  aucs <- rf_model$resample$ROC
  return(aucs)
}

train_cv_auc <- function(genus_summed, response) {
    set.seed(100) 
    folds <- createFolds(response, k = 5, returnTrain = TRUE)

    cv_results <- lapply(names(folds), function(fold_name) {

      train_index <- folds[[fold_name]]

      train_x <- genus_summed[train_index, ]
      train_y <- response[train_index]

      test_x <- genus_summed[-train_index, ]
      test_y <- response[-train_index]
      
      model_fit <- train(
          x = train_x,
          y = train_y,
          method = "rf",
          metric = "ROC",
          trControl = trainControl(method="none", classProbs = TRUE),
          nodesize = 5
      )

      pred_train <- predict(model_fit, train_x, type = "prob")
      roc_train <- roc(train_y, pred_train[, levels(response)[2]], 
                      levels = levels(response), direction = "<", quiet = TRUE)
      auc_train <- as.numeric(roc_train$auc)

      pred_test <- predict(model_fit, test_x, type = "prob")
      roc_test <- roc(test_y, pred_test[, levels(response)[2]], 
                      levels = levels(response), direction = "<", quiet = TRUE)
      auc_test <- as.numeric(roc_test$auc)
      
      return(c(Train_AUC = auc_train, Test_AUC = auc_test))
      })
    cv_results_df <- do.call(rbind, cv_results) %>% as.data.frame()
    cv_results_df$Train_AUC <- as.numeric(cv_results_df$Train_AUC)
    cv_results_df$Test_AUC <- as.numeric(cv_results_df$Test_AUC)
    return(cv_results_df$Train_AUC)
}

glm_cv <- function(genus_summed, response, sig_feature) {
    my_split <- function(response, t = 5, seed = 100) {
        values <- unique(response)
        i1 <- which(response == values[1])
        i2 <- which(response == values[2])
        k1 <- round(length(i1) / t)
        k2 <- round(length(i2) / t)
        set.seed(seed)
        res <- NULL
        for (i in 1 : (t-1)) {
          j1 <- sample(i1, k1)
          j2 <- sample(i2, k2)
          i1 <- setdiff(i1, j1)
          i2 <- setdiff(i2, j2)
          res[[i]] <- c(j1, j2)
        }
        res[[t]] <- c(i1, i2)
        return(res)
    }
    up_data_full <- upSample(x = genus_summed, y = response)
    X_train_full <- up_data_full[, -ncol(up_data_full)]
    y_train_full <- up_data_full$Class
    feature <- X_train_full[, sig_feature]
    train_df <- data.frame(feature)
    colnames(train_df) <- sig_feature
    train_df$response <- y_train_full
    train_df$response <- as.factor(train_df$response)
    index_list <- my_split(y_train_full, 5, 100)
    aucs <- c()
    for (i in 1 : length(index_list)) {
        sub_train_df <- train_df[-index_list[[i]], ]
        sub_val_df <- train_df[index_list[[i]], ]
        model <- glm(response ~., data = sub_train_df, family = binomial())
        train_predictions <- predict(model, newdata = sub_val_df, type = "response")
        val_roc <- roc(y_train_full[index_list[[i]]], train_predictions, quiet = TRUE)
        aucs <- c(aucs, auc(val_roc))
      }
    return(aucs)
}

train_glm <- function(genus_summed, response, sig_feature) {
    my_split <- function(response, t = 5, seed = 100) {
        values <- unique(response)
        i1 <- which(response == values[1])
        i2 <- which(response == values[2])
        k1 <- round(length(i1) / t)
        k2 <- round(length(i2) / t)
        set.seed(seed)
        res <- NULL
        for (i in 1 : (t-1)) {
          j1 <- sample(i1, k1)
          j2 <- sample(i2, k2)
          i1 <- setdiff(i1, j1)
          i2 <- setdiff(i2, j2)
          res[[i]] <- c(j1, j2)
        }
        res[[t]] <- c(i1, i2)
        return(res)
    }
    up_data_full <- upSample(x = genus_summed, y = response)
    X_train_full <- up_data_full[, -ncol(up_data_full)]
    y_train_full <- up_data_full$Class
    feature <- X_train_full[, sig_feature]
    train_df <- data.frame(feature)
    colnames(train_df) <- sig_feature
    train_df$response <- y_train_full
    train_df$response <- as.factor(train_df$response)
    index_list <- my_split(y_train_full, 5, 100)
    aucs <- c()
    for (i in 1 : length(index_list)) {
        sub_train_df <- train_df[-index_list[[i]], ]
        model <- glm(response ~., data = sub_train_df, family = binomial())
        train_predictions <- predict(model, newdata = sub_train_df, type = "response")
        train_roc <- roc(y_train_full[-index_list[[i]]], train_predictions, quiet = TRUE)
        aucs <- c(aucs, auc(train_roc))
      }
    return(aucs)
}

res_train <- data.frame()
for (i in 1:length(selected_feature)) {
  sig_feature2 <- selected_feature[1:i]
  tr_res <- deal_compares_tr(compares, feature_summed, response)
  if (i <= 1) {
    cv_res <- train_glm(tr_res$feature, tr_res$response, selected_feature[i])
  }else {
    cv_res <- get_cv_auc(tr_res$feature[, sig_feature2], tr_res$response, 6)
  }
  print(paste0("Feature: ", selected_feature[i]))
  tmp_df <- data.frame(Feature = selected_feature[i],
    AUCs = cv_res)
  res_train <- rbind(res_train, tmp_df)
  mean_auc <- mean(cv_res)
  sd <- sd(cv_res)
  ci_lower <- mean_auc - (sd / sqrt(length(cv_res)))
  ci_upper <- mean_auc + (sd / sqrt(length(cv_res)))
}



my_plot <- function(res_train, res_test, selected_feature) {
    summary_test <- res_test %>%
    group_by(Feature) %>%
    summarise(
        mean_auc = mean(AUCs),
        sd_auc = sd(AUCs),
        n = n(),
        se_auc = sd_auc / sqrt(n),
        ci_lower = mean_auc - se_auc,
        ci_upper = mean_auc + se_auc
    ) %>%
    mutate(type = "Validation")

    summary_train<- res_train %>%
    group_by(Feature) %>%
    summarise(
        mean_auc = mean(AUCs),
        sd_auc = sd(AUCs),
        n = n(),
        se_auc = sd_auc / sqrt(n),
        ci_lower = mean_auc - se_auc,
        ci_upper = mean_auc + se_auc
    ) %>%
    mutate(type = "Train")

    plot_df$type <- factor(plot_df$type, levels = c("Train", "Validation"))
    plot_df <- rbind(summary_train, summary_test)
    plot_df$Feature <- factor(plot_df$Feature, levels = selected_feature)
    
    return(plot_df)
}

res_test <- data.frame()
for (i in 1:length(selected_feature)) {
  sig_feature2 <- selected_feature[1:i]
  tr_res <- deal_compares_tr(compares, feature_summed, response)
  if (i <= 1) {
    cv_res <- glm_cv(tr_res$feature, tr_res$response, selected_feature[i])
  }else {
    cv_res <- get_cv_auc(tr_res$feature[, sig_feature2], tr_res$response, 6)
  }
  print(paste0("Feature: ", selected_feature[i]))
  tmp_df <- data.frame(Feature = selected_feature[i],
    AUCs = cv_res)
  res_test <- rbind(res_test, tmp_df)
  mean_auc <- mean(cv_res)
  sd <- sd(cv_res)
  ci_lower <- mean_auc - (sd / sqrt(length(cv_res)))
  ci_upper <- mean_auc + (sd / sqrt(length(cv_res)))
}

selected_feature <- sub("^g_", "", selected_feature)
plot_df <- my_plot(res_train, res_test, selected_feature)

pdf(paste0("figs/fig4/com", compare_num, ".pdf"), width = 4.7, height = 3)
ggplot(plot_df, aes(x = Feature, y = mean_auc, group = type)) +
  geom_point(aes(color = type), size = 1.5) +
  geom_line(aes(color = type), linewidth = 0.5) +
  geom_errorbar( data = plot_df[plot_df$type == "Train", ],
    aes(x = Feature, ymin = ci_lower, ymax = ci_upper),
    width = 0.3,
    color = "darkred",
    linewidth = 0.4
  ) +
  geom_text(data = plot_df[plot_df$type == "Train", ],
    aes(x = Feature, label = sprintf("%.3f", mean_auc)),
    vjust = 2,
    hjust = 0.25,
    size = 3
  ) +
  ggsci::scale_color_d3() +
  labs(
    y = "AUC"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = -45, vjust = -0.5, hjust = 0.5),
    axis.line = element_line(linewidth = 0.5)
  ) +
  scale_y_continuous(limits = c(0.8, 1))
dev.off()
