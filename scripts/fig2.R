library(dplyr)
library(vegan)
library(ggplot2)
library(microeco)
library(magrittr)
library(ggtree)
library(randomForest)
library(pROC)
library(metagenomeSeq)
source("scripts/methods.R")
#load the data
{
  genus_summed <- readRDS("data/genus_summed.rds")
  genus_summed <- genus_summed[which(substr(rownames(genus_summed), 1, 1) %in% c("N", "a", "T")), ]
  genus_summed <- genus_summed[, which(colSums(genus_summed > 0) > 0.05*nrow(genus_summed))]
  feature_summed <- genus_summed
  patient_info <- readRDS("data/patient_info.rds")
}

#PCoA
for (problem in c("I.II.III.IV", "HD_MD", "Location", "age", "gender", "survival")) {
    deal_res <- deal_problem(feature_summed, patient_info, problem)
    response <- deal_res$response
    predictors <- deal_res$feature
    otu <- predictors
    otu.distance <- vegdist(otu)
    pcoa <- cmdscale(otu.distance,eig=TRUE)
    pc12 <- pcoa$points[,1:2]
    pc <- round(pcoa$eig/sum(pcoa$eig)*100,digits=2)
    pc12 <- as.data.frame(pc12)
    df <- cbind(pc12, data.frame(Group = response))

    p1 <- ggplot(data=df,aes(x=V1,y=V2,color=Group))+#指定数据、X轴、Y轴，颜色
        geom_point(size=2) +#绘制点图并设定大小
        theme_classic(base_size = 14) +
        labs(x=paste0("PCoA1 (",pc[1],"%)"),
            y=paste0("PCoA2 (",pc[2],"%)"))+#将x、y轴标题改为贡献度
        stat_ellipse(data=df,
                geom = "polygon",level=0.9,
                linetype = 1,size=0.8,
                aes(fill=Group),
                alpha=0.1,
                show.legend = T)  +
                theme(legend.position = "top")
    pdf(paste0("figs/fig2/pcoa", problem, "_genus.pdf"), width = 4, height = 4)
      print(p1)
    dev.off()
}


#ROC plot
{
    train_genus <- readRDS("data/train_genus.rds")
    colnames(train_genus) <- paste0("g_", colnames(train_genus), "")
    train_response <- readRDS("data/train_response.rds")

    test_genus <- readRDS("data/test_genus.rds")
    colnames(test_genus) <- paste0("g_", colnames(test_genus), "")
    test_response <- readRDS("data/test_response.rds")

    train_feature <- cbind(train_genus)
    test_feature <- cbind(test_genus)
    all_feature <- rbind(train_feature, test_feature)

    ms_obj <- newMRexperiment(t(train_feature))
    ms_obj <- cumNorm(ms_obj, p = 0.5)
    train_feature <- t(MRcounts(ms_obj, norm = TRUE, log = FALSE))

    ms_obj <- newMRexperiment(t(test_feature))
    ms_obj <- cumNorm(ms_obj, p = 0.5)
    test_feature <- t(MRcounts(ms_obj, norm = TRUE, log = FALSE))

    patient_info <- readRDS("data/patient_info.rds")
}


#ROC in train set
cv_roc_mean <- function(response, predictors) {
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
      index_list <- my_split(response, 5, 100)
      roc_predictors <- NULL
      roc_response <- NULL
      for (i in 1 : length(index_list)) {
        rf_model <- randomForest(x = predictors[-index_list[[i]], ],
            y = response[-index_list[[i]]], importance=TRUE, ntree = 500)
        pred_prob <- predict(rf_model, newdata=predictors[index_list[[i]], ], type="prob")
        positive_probs <- pred_prob[, 2]
        roc_predictors[[i]] <- positive_probs
        roc_response[[i]] <- response[index_list[[i]]]
      }
      roc_predictors <- do.call(c, roc_predictors)
      roc_response2 <- do.call(c, roc_response)
      roc_obj <- roc(roc_response2, roc_predictors)
      roc_data <- data.frame(
        tpr = roc_obj$sensitivities,
        fpr = roc_obj$specificities,
        thresholds = roc_obj$thresholds
      )
    roc_data <- roc_data[order(roc_data$tpr), ]
    return(list(roc_data = roc_data, auc = roc_obj$auc))
}

problem_list <- c("Location", "I.II_III.IV", "I.II.III_IV", "HD_MD", "survival")
roc_problems <- data.frame()
name_ps <- c()
for (problem in problem_list) {
    deal_res <- deal_problem(train_feature, patient_info, problem)
    response <- deal_res$response
    predictors <- deal_res$feature
    res_roc <- cv_roc_mean(response, predictors)
    roc_data <- res_roc$roc_data
    auc <- round(res_roc$auc, 3)
    roc_data$problem <- rep(paste0(problem, "(AUC = ", auc, ")"), dim(roc_data)[1])
    name_ps <- c(name_ps, paste0(problem, "(AUC = ", auc, ")"))
    roc_problems <- rbind(roc_problems, roc_data)
}
roc_problems$problem <- factor(roc_problems$problem, levels = name_ps)
pdf("figs/fig2/ROC_Problem_in_CRC_train.pdf", height = 3, width = 5)
ggplot(roc_problems, aes(x = 1 - fpr, y = tpr, color = problem)) +
  geom_line(size = 0.75) +
  geom_abline(linetype = "dashed") +
  labs(x = "1 - Specificity", y = "Sensitivity") +
  scale_color_npg() +
  theme_classic()
dev.off()

#ROC in test set
test_roc <- function(train_feature, train_response, test_feature, test_response) {
  #up_data_full <- upSample(x = train_feature, y = train_feature)
  #X_train_full <- up_data_full[, -ncol(up_data_full)]
  #y_train_full <- up_data_full$Class
  
  final_model <- randomForest(x = train_feature, y = train_response)
  pred_prob <- predict(final_model, newdata = test_feature, type="prob")
  positive_probs <- pred_prob[, 2]
  roc_obj <- roc(test_response, positive_probs)
  roc_data <- data.frame(
    tpr = roc_obj$sensitivities,
    fpr = roc_obj$specificities,
    thresholds = roc_obj$thresholds
  )
  roc_data <- roc_data[order(roc_data$tpr), ]
  return(list(roc_data = roc_data, auc = roc_obj$auc))
}

problem_list <- c("Location", "I.II_III.IV", "I.II.III_IV", "HD_MD", "survival")
roc_problems <- data.frame()
name_ps <- c()
for (problem in problem_list) {
    deal_res <- deal_problem(train_feature, patient_info, problem)
    response <- deal_res$response
    predictors <- deal_res$feature
    test_res <- deal_problem(test_feature, patient_info, problem)
    res_roc <- test_roc(predictors, response, test_res$feature, test_res$response)
    roc_data <- res_roc$roc_data
    auc <- round(res_roc$auc, 3)
    roc_data$problem <- rep(paste0(problem, "(AUC = ", auc, ")"), dim(roc_data)[1])
    name_ps <- c(name_ps, paste0(problem, "(AUC = ", auc, ")"))
    roc_problems <- rbind(roc_problems, roc_data)
}
roc_problems$problem <- factor(roc_problems$problem, levels = name_ps)
pdf("figs/fig2/ROC_Problem_in_CRC_test.pdf", height = 3, width = 5)
ggplot(roc_problems, aes(x = 1 - fpr, y = tpr, color = problem)) +
  geom_line(size = 0.75) +
  geom_abline(linetype = "dashed") +
  labs(x = "1 - Specificity", y = "Sensitivity") +
  scale_color_npg() +
  theme_classic()
dev.off()


#Heat map
library(pheatmap)
g_sum <- colSums(genus_summed)

patient_info <- readRDS("data/patient_info.rds")
top_genus <- names(g_sum)[order(g_sum, decreasing = TRUE)[1:30]]

g_sum[top_genus]

problem_list <- c("N_and_T", "gender", "age", "Location", "I.II.III.IV", "survival")

mean_df <- matrix(0, 30, 14)
rownames(mean_df) <- top_genus
colnames(mean_df) <- paste0("C", 1:14)
i <- 1
genus_summed <- log(genus_summed + 1)
for (problem in problem_list) {
  deal_res <- deal_problem(genus_summed, patient_info, problem)
  feature_summed <- deal_res$feature
  response <- deal_res$response
  for (class_name in levels(response)) {
    sub_genus <- feature_summed[response == class_name, top_genus]
    mean_df[, i] <- as.numeric(colMeans(sub_genus))
    colnames(mean_df)[i] <- class_name
    print(class_name)
    i <- i + 1
  }
}

rownames(mean_df) <- paste0("g_", top_genus)
colnames(mean_df)[3:4] <- c("Female", "Male")
pdf(paste0("figs/fig2/mean_top_genus_heatmap.pdf"), width = 7, height = 6)
ph <- pheatmap(mean_df, cluster_rows = TRUE, cluster_cols = FALSE, scale = "none", show_colnames = TRUE,
            show_rownames = TRUE, fontsize = 12)
dev.off()