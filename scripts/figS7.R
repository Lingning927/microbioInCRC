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


{
  genus_summed <- readRDS("adata2025/data/genus_summed.rds")
  genus_summed <- genus_summed[which(substr(rownames(genus_summed), 1, 1) %in% c("N", "a", "T")), ]

  train_feature <- cbind(genus_summed)

  ms_obj <- newMRexperiment(t(train_feature))
  ms_obj <- cumNorm(ms_obj, p = 0.5)
  otu_table_norm <- t(MRcounts(ms_obj, norm = TRUE, log = TRUE))
  train_feature <- otu_table_norm
  feature_summed <- train_feature

  response <- substr(rownames(genus_summed), 1, 1)
  response = case_when(
        response == "N" ~ "Normal",
        response == "a" ~ "Polyp",
        response == "T" ~ "CRC"
    )
  response <- factor(response, levels = c("Normal", "Polyp", "CRC"))
}

normal_names <- rownames(feature_summed)[11:68]
crc_names <- sub("^N", "T", rownames(feature_summed)[11:68])

crc_names <- crc_names[crc_names  %in% rownames(feature_summed)]
normal_names <- sub("^T", "N", crc_names)


train_cv_auc <- function(genus_summed, response) {
    # 1. 手动创建固定的 Folds (5折)
    set.seed(100) # 设置种子以保证结果可复现
    folds <- createFolds(response, k = 5, returnTrain = TRUE)

    # 2. 定义一个循环函数，遍历每一折
    cv_results <- lapply(names(folds), function(fold_name) {
    
    # 获取当前折的训练集索引
    train_index <- folds[[fold_name]]
    
    # 分割数据：训练集 (Current Fold Training)
    train_x <- genus_summed[train_index, ]
    train_y <- response[train_index]
    
    # 分割数据：验证集 (Current Fold Validation)
    test_x <- genus_summed[-train_index, ]
    test_y <- response[-train_index]
    
    # 3. 在当前训练集上训练模型
    # 注意：为了加快速度并只看训练效果，这里我使用了 method="none" (即不进行内部参数调优)
    # 如果你需要在此过程中也进行参数调优，可以去掉 trControl 中的 method="none"，但速度会变慢
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


get_loo_auc <- function(genus_summed, response) {
  response <- as.factor(response)
  ref_level <- levels(response)[2]
  response <- relevel(response, ref = ref_level)
  train_control <- trainControl(
    method = "LOOCV",
    sampling = "up",           
    classProbs = TRUE,         
    savePredictions = "final", 
    allowParallel = TRUE
  )
  if (levels(response) == c("Normal", "CRC")) {
    pred_all <- rep(0, nrow(genus_summed))
    for (i in 1:nrow(genus_summed)) {
        if (rownames(genus_summed)[i] %in% c(crc_names, normal_names)) {
            if (response[i] == "Normal") {
                crc_tmp <- sub("^N", "T", rownames(genus_summed)[i])
                test_index <- c(i, which(rownames(genus_summed) == crc_tmp))
            }else {
                crc_tmp <- sub("^T", "N", rownames(genus_summed)[i])
                test_index <- c(i, which(rownames(genus_summed) == crc_tmp))
            }
        }else {
            test_index <- i
        }
        train_x <- genus_summed[-test_index, ]
        train_y <- response[-test_index]
        test_x <- genus_summed[test_index, ]
        test_y <- response[test_index]
        model_fit <- train(
            x = train_x,
            y = train_y,
            method = "rf",
            metric = "ROC",
            trControl = trainControl(method="none", classProbs = TRUE),
            nodesize = 5
        )
        pred_test <- predict(model_fit, test_x, type = "prob")
        pred_all[test_index] <- pred_test[, levels(response)[2]]
    }
    roc_val <- roc(response, pred_all, 
                        levels = levels(response), direction = "<", quiet = TRUE)
    val_auc <- as.numeric(roc_val$auc)
  }else {
    set.seed(100)
    rf_model <- train(
        x = genus_summed,
        y = response,
        method = "rf",
        metric = "Accuracy",
        trControl = train_control
    )
    selected_probs <- rf_model$pred[, ref_level]
    roc_obj <- roc(response = rf_model$pred$obs, 
                    predictor = selected_probs,
                    levels = rev(levels(response)))
    val_auc <- as.numeric(auc(roc_obj))
  }
  return(val_auc)
}

compares_list <- list(c("Normal", "Polyp"), c("Polyp", "CRC"),  c("Normal", "CRC"), c("Polyp", "CRC_1"))
polyp_crc1 <- c("Pseudomonas", "Megamonas", "Phascolarctobacterium", "Veillonella", "Parvimonas", "Leptotrichia")
sig_feature_list <- list(normal_polyp, polyp_crc, normal_crc, polyp_crc)

res_loo <- data.frame()
for (i in 1:4) {
  compares <- compares_list[[i]]
  sig_feature <- sig_feature_list[[i]]
  tr_res <- deal_compares_tr(compares, train_feature, response)
  loo_auc <- get_loo_auc(tr_res$feature[, sig_feature], tr_res$response)
  print(paste0("Comparison: ", paste(compares, collapse = " vs ")))
  tmp_df <- data.frame(Compare = paste(compares, collapse = " vs "),
    AUC = loo_auc)
  res_loo <- rbind(res_loo, tmp_df)
  print(paste0("LOO AUC: ", round(loo_auc, 4)))
}

saveRDS(res_loo, "data/res_loo.rds")


res_loo <- readRDS("data/res_loo.rds")
pdf("figs/figS7/loo_barplot.pdf", width = 3, height = 3)
ggplot(res_loo, aes(x = Compare, y = AUC)) +
    geom_col(fill = "steelblue", position="dodge",width =0.7) +
    geom_text(
        aes(label = sprintf("%.3f", AUC)),
        vjust = -1,
        size = 3
    ) +
    theme_classic() +
    theme(
    axis.text.x = element_text(angle = -45, vjust = -0.5, hjust = 0.5),
    axis.line = element_line(linewidth = 0.5)
  ) + 
  scale_y_continuous(limits = c(0, 1.1), expand = c(0, 0))
dev.off()
