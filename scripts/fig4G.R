library(data.table)
library(dplyr)
library(metagenomeSeq)
library(randomForest)
library(themis)
library(caret)
library(pROC)
source("scripts/methods.R")

tb <- fread("data/changhai_otus.csv")
patient_id <- colnames(tb)[2:(ncol(tb) - 7)]
patient_info <- tb[1:7, 2:(ncol(tb) - 7)]
colnames(patient_info) <- patient_id
rownames(patient_info) <- as.character(t(tb[1:7, 1]))
patient_info <- as.data.frame(patient_info)

bio_id <- as.character(t(tb[9:nrow(tb), 1]))

bio_cnt <- tb[9:nrow(tb), 2:(ncol(tb) - 7)]
rownames(bio_cnt) <- bio_id
colnames(bio_cnt) <- patient_id
bio_cnt <- as.data.frame(bio_cnt)

bio_cat <- tb[9:nrow(tb), (ncol(tb) - 6):ncol(tb)]
colnames(bio_cat) <- as.character(tb[8, (ncol(tb) - 6):ncol(tb)])
rownames(bio_cat) <- bio_id
bio_cat <- as.data.frame(bio_cat)

bio_cnt_numeric <- bio_cnt %>%
  mutate(across(everything(), as.numeric))

genus_summed_new <- merge_to_species(bio_cnt_numeric, bio_cat, "Genus")

genus_summed <- readRDS("data/genus_summed.rds")
genus_summed <- genus_summed[which(substr(rownames(genus_summed), 1, 1) %in% c("N", "a", "T")), ]

train_feature <- cbind(genus_summed)
test_feature <- cbind(genus_summed_new)

ms_obj <- newMRexperiment(t(train_feature))
ms_obj <- cumNorm(ms_obj, p = 0.5)
otu_table_norm <- t(MRcounts(ms_obj, norm = TRUE, log = TRUE))
otu_table_pseudo <- otu_table_norm + 1

ms_obj <- newMRexperiment(t(test_feature))
ms_obj <- cumNorm(ms_obj, p = 0.5)
otu_table_norm <- t(MRcounts(ms_obj, norm = TRUE, log = TRUE))
otu_table_pseudo <- otu_table_norm + 1

response <- substr(rownames(genus_summed), 1, 1)
response = case_when(
        response == "N" ~ "Normal",
        response == "a" ~ "Polyp",
        response == "T" ~ "CRC"
)
response <- factor(response, levels = c("Normal", "Polyp", "CRC"))
test_response <- as.character(patient_info[7, ])
test_response <- factor(test_response, levels = c("Normal", "Polyp", "CRC"))


normal_polyp <- c("Fusobacterium", "Pseudomonas", "Pigmentiphaga", "Collinsella", "Parabacteroides", "Bacteroides")
normal_crc <- c("Fusobacterium", "Pseudomonas", "Collinsella", "Phascolarctobacterium", "Granulicatella", "Leptotrichia", "Streptococcus")
polyp_crc <- c("Parvimonas", "Collinsella", "Pseudomonas", "Fusobacterium")
compares_list <- list(c("Normal", "Polyp"), c("Polyp", "CRC"),  c("Normal", "CRC"))
sig_feature_list <- list(normal_polyp, polyp_crc, normal_crc)

roc_problems <- data.frame()
name_ps <- c()
for (i in 1:3) {
  compares <- compares_list[[i]]
  sig_feature <- sig_feature_list[[i]]
  tr_res <- deal_compares_tr(compares, train_feature, response)
  te_res <- deal_compares2(compares, test_feature, test_response)
  rf_res <- my_rf_train(tr_res, te_res, sig_feature)
  roc_test <- rf_res$te_data
  roc_test <- roc_test[order(roc_test$tpr), ]
  auc <- round(rf_res$te_auc, 3)
  problem <- paste(compares[1], "vs", compares[2])
  roc_test$problem <- rep(paste0(problem, "(AUC = ", auc, ")"), dim(roc_test)[1])
  name_ps <- c(name_ps, paste0(problem, "(AUC = ", auc, ")"))
  roc_problems <- rbind(roc_problems, roc_test)
  print(paste0("Comparison: ", paste(compares, collapse = " vs "), 
               ", Train AUC: ", round(rf_res$tr_auc, 3),
               ", Test AUC: ", round(rf_res$te_auc, 3)))
}

roc_problems$problem <- factor(roc_problems$problem, levels = name_ps)

pdf("figs/fig4/Normal_Polyp_CRC_test.pdf", height = 2.4, width = 4.8)
ggplot(roc_problems, aes(x = 1 - fpr, y = tpr, color = problem)) +
  geom_line(size = 1) +
  geom_abline(linetype = "dashed") +
  labs(x = "1 - Specificity", y = "Sensitivity") +
  theme_classic()
dev.off()
