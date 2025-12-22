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
  genus_summed <- readRDS("MicrobioInCRC/data/genus_summed.rds")
  genus_summed <- genus_summed[which(substr(rownames(genus_summed), 1, 1) %in% c("N", "a", "T")), ]
  genus_summed <- genus_summed[, which(colSums(genus_summed > 0) > 0.05*nrow(genus_summed))]
  feature_summed <- genus_summed
  patient_info <- readRDS("MicrobioInCRC/data/patient_info.rds")
}

#paired CRC samples
normal_names <- rownames(feature_summed)[11:68]
crc_names <- sub("^N", "T", rownames(feature_summed)[11:68])
feature_summed <- feature_summed[which(rownames(feature_summed) %in% crc_names), ]

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
    pdf(paste0("figs/figS3/Sub_pcoa", problem, "_genus.pdf"), width = 4, height = 4)
      print(p1)
    dev.off()
    colnames(df)[1:2] <- c("PCoA1", "PCoA2")
    write.csv(df, paste0("figdata/figS2_", problem, ".csv"))
}


#ROC plot
{
    all_feature <- feature_summed
    ms_obj <- newMRexperiment(t(all_feature))
    ms_obj <- cumNorm(ms_obj, p = 0.5)
    all_feature <- t(MRcounts(ms_obj, norm = TRUE, log = FALSE))

    patient_info <- patient_info[, rownames(all_feature)]
}
feature_summed <- feature_summed[crc_names, ]

#ROC in train set
problem_list <- c("Location", "I.II_III.IV", "I.II.III_IV", "HD_MD", "survival")
roc_problems <- data.frame()
name_ps <- c()
for (problem in problem_list) {
    deal_res <- deal_problem(all_feature, patient_info, problem)
    response <- deal_res$response
    predictors <- deal_res$feature
    res_roc <- cv_roc_mean(response, predictors)
    roc_data <- res_roc$roc_data
    auc <- round(res_roc$auc, 3)
    roc_data$problem <- rep(paste0(problem, "(AUC = ", auc, ")"), dim(roc_data)[1])
    print(paste0(problem, "(AUC = ", auc, ")"))
    name_ps <- c(name_ps, paste0(problem, "(AUC = ", auc, ")"))
    roc_problems <- rbind(roc_problems, roc_data)
}
roc_problems$problem <- factor(roc_problems$problem, levels = name_ps)
pdf("figs/figS2/ROC_Problem_in_CRC.pdf", height = 3, width = 5)
ggplot(roc_problems, aes(x = 1 - fpr, y = tpr, color = problem)) +
  geom_line(size = 0.75) +
  geom_abline(linetype = "dashed") +
  labs(x = "1 - Specificity", y = "Sensitivity") +
  scale_color_npg() +
  theme_classic()
dev.off()
