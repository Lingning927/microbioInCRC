library(ggplot2)
library(vegan)
library(ggpubr)
library(dplyr)
library(rentrez)
source("scripts/methods.R")

# Load the data
{
genus_summed <- readRDS("data/genus_summed.rds")
genus_summed <- genus_summed[which(substr(rownames(genus_summed), 1, 1) %in% c("N", "a", "T")), ]

genus_summed <- log(genus_summed + 1)
genus_summed <- genus_summed / rowSums(genus_summed)
genus_summed <- genus_summed[, which(colSums(genus_summed > 0) > 0.05*nrow(genus_summed))]

feature_summed <- genus_summed
response <- substr(rownames(feature_summed), 1, 1)
response = case_when(
    response == "N" ~ "Normal",
    response == "a" ~ "Polyp",
    response == "T" ~ "CRC"
  )
response <- factor(response, levels = c("Normal", "Polyp", "CRC"))
}

load("data/CRC_Amplicon_genus.RData")
feat_list <- my_data0$feat_list
meta_list <- my_data0$meta_list

feat1 <- feat_list[[1]]
feat2 <- feat_list[[2]]

features_feat <- union(colnames(feat1), colnames(feat2))

ids <- sub("ncbi_", "", features_feat)


get_tax_name <- function(id) {
  tax_info <- entrez_summary(db = "taxonomy", id = id)
  tax_info$scientificname
}

names_feat <- sapply(ids, get_tax_name)
id_to_name <- unlist(names_feat)

id1 <- sub("ncbi_", "", colnames(feat1))
name1 <- as.character(id_to_name[id1])
colnames(feat1) <- name1
id2 <- sub("ncbi_", "", colnames(feat2))
name2 <- as.character(id_to_name[id2])
colnames(feat2) <- name2

feature1 <- colnames(genus_summed)
cleaned_string <- sub("_unclassified$", "", feature1)
feature1 <- feature1[which(feature1 %in% cleaned_string)]
feature_summed <- feature_summed[, feature1]
feature_summed <- feature_summed[which(substr(rownames(feature_summed), 1, 1) %in% c("N", "T")), ]

features2 <- as.character(unlist(names_feat))
all_features <- union(feature1, features2)

lack_features <- setdiff(all_features, colnames(feat1))
lack_feat <- matrix(0, nrow = nrow(feat1), ncol = length(lack_features))
colnames(lack_feat) <- lack_features
feat1 <- cbind(feat1, lack_feat)
feat1 <- feat1[, all_features]

lack_features <- setdiff(all_features, colnames(feat2))
lack_feat <- matrix(0, nrow = nrow(feat2), ncol = length(lack_features))
colnames(lack_feat) <- lack_features
feat2 <- cbind(feat2, lack_feat)
feat2 <- feat2[, all_features]

feat <- rbind(feat1, feat2)

lack_features <- setdiff(all_features, colnames(feature_summed))
lack_feat <- matrix(0, nrow = nrow(feature_summed), ncol = length(lack_features))
colnames(lack_feat) <- lack_features
feature_summed <- cbind(feature_summed, lack_feat)
feature_summed <- feature_summed[, all_features]

feat <- feat / rowSums(feat)
feature_summed <- feature_summed / rowSums(feature_summed)

meta_response <- c(meta_list[[1]]$Group, meta_list[[2]]$Group)
response <- as.character(response[which(response %in% c("Normal", "CRC"))])

tmp_data <- list(feat = feat, feature_summed = feature_summed, meta_response = meta_response, response = response)

saveRDS(tmp_data, file = "data/CRC_genus.rds")



tmp_data <- readRDS("data/CRC_genus.rds")
feat <- tmp_data$feat
feature_summed <- tmp_data$feature_summed
meta_response <- tmp_data$meta_response
response <- tmp_data$response

CRC_matrix1 <- feat[meta_response == "Case", ]
CRC_matrix2 <- feature_summed[response == "CRC", ]
CRC_matrix <- rbind(CRC_matrix1, CRC_matrix2)
CRC_response <- c(rep("Fecal", nrow(CRC_matrix1)), rep("Intratumoral", nrow(CRC_matrix2)))

Normal_matrix1 <- feat[meta_response == "Control", ]
Normal_matrix2 <- feature_summed[response == "Normal", ]
Normal_matrix <- rbind(Normal_matrix1, Normal_matrix2)
Normal_response <- c(rep("Fecal", nrow(Normal_matrix1)), rep("Intratumoral", nrow(Normal_matrix2)))

pcoa_df <- function(CRC_matrix, response) {
    otu <- CRC_matrix
     otu.distance <- vegdist(otu)
     #pcoa分析
     pcoa <- cmdscale(otu.distance,eig=TRUE)
     pc12 <- pcoa$points[,1:2]
     pc <- round(pcoa$eig/sum(pcoa$eig)*100,digits=2)#解释度
     #pc12原来是matrix,转化为data.frame
     pc12 <- as.data.frame(pc12)
     #给pc12添加samp1es变量
     df <- cbind(pc12, data.frame(Group = response))
     p1 <- ggplot(data=df,aes(x=V1,y=V2,color=Group))+#指定数据、X轴、Y轴，颜色
     geom_point(size=1.5) +#绘制点图并设定大小
     theme_classic(base_size = 14) +
     scale_color_manual(values = c("#1B9E77", "#D95F02")) +
     labs(x=paste0("PCoA1 (",pc[1],"%)"),
          y=paste0("PCoA2 (",pc[2],"%)"))+#将x、y轴标题改为贡献度
     stat_ellipse(data=df,
               geom = "polygon",level=0.9,
               linetype = 1,size=0.6,
               aes(fill=Group),
               alpha=0.1,
               show.legend = T)  +
               theme(legend.position = "top")
    return(p1)
}

p1 <- pcoa_df(CRC_matrix, CRC_response)
p2 <- pcoa_df(Normal_matrix, Normal_response)

pdf(paste0("figs/figS1/", "pcoa_CRC.pdf"), width = 4, height = 4)
    print(p1)
dev.off()
pdf(paste0("figs/figS1/", "pcoa_Normal.pdf"), width = 4, height = 4)
    print(p2)
dev.off()


calculate_shannon <- function(matrix_data) {
  data.frame(
    Sample = rownames(matrix_data),
    Shannon = diversity(matrix_data, index = "shannon"),
    Group = deparse(substitute(matrix_data))
  )
}

shannon_crc1 <- calculate_shannon(CRC_matrix1)
shannon_crc2 <- calculate_shannon(CRC_matrix2) 

combined_shannon <- bind_rows(shannon_crc1, shannon_crc2) %>%
  mutate(Group = CRC_response)

wilcox_test <- wilcox.test(Shannon ~ Group, data = combined_shannon)
p_value <- format.pval(wilcox_test$p.value, digits = 3, eps = 0.001)
p_label <- ifelse(wilcox_test$p.value < 0.001, "***",
                  ifelse(wilcox_test$p.value < 0.01, "**",
                         ifelse(wilcox_test$p.value < 0.05, "*", "ns")))

pdf("figs/figS1/shannon_compare_feat.pdf", width = 3, height = 4)
ggplot(combined_shannon, aes(x = Group, y = Shannon, fill = Group)) +

  geom_boxplot(width = 0.5, alpha = 0.8, outlier.shape = NA) +
  geom_signif(
    comparisons = list(c("Gut", "Intratumoral")),
    annotations = paste0("Wilcox, P ", p_value),
    y_position = max(combined_shannon$Shannon) * 1.05,
    tip_length = 0.01,                                 
    textsize = 5, 
    vjust = -0.5                                      
  ) +
  scale_fill_manual(values = c("#66C2A5", "#FC8D62")) +
  scale_color_manual(values = c("#1B9E77", "#D95F02")) +
  scale_y_continuous(limits = c(0, max(combined_shannon$Shannon) * 1.2)) +
  labs(
    x = "Group",
    y = "Shannon Index"
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "none", 
    axis.text.x = element_text(size = 12, color = "black")
  )
dev.off()


prepare_long_data <- function(matrix_data, group_name) {
  samples <- rownames(matrix_data)
  res <- data.frame(matrix_data) %>%
  mutate(Sample = samples)

  res <- res %>%
    tidyr::pivot_longer(
      cols = -"Sample",
      names_to = "Taxon",
      values_to = "Abundance"
    ) %>%
    mutate(Group = group_name)
    res
}

crc1_long <- prepare_long_data(CRC_matrix1, "CRC1")
crc2_long <- prepare_long_data(CRC_matrix2, "CRC2")
combined_long <- bind_rows(crc1_long, crc2_long)

length(unique(colnames(CRC_matrix1)))

################

library(viridis)

prepare_group_data <- function(matrix_data, group_name) {
    mean_abumdance <- colMeans(matrix_data)
    top_taxa <- names(sort(mean_abumdance, decreasing = TRUE)[1:15])
    abundance <- mean_abumdance[top_taxa]
    abundance <- c(abundance, sum(mean_abumdance[-which(names(mean_abumdance) %in% top_taxa)]))
    names(abundance)[length(abundance)] <- "Others"
    data.frame(
      Taxon = names(abundance),
      Abundance = abundance,
      Group = group_name
    )
}


crc1 <- prepare_group_data(CRC_matrix1, "Fecal")
crc2 <- prepare_group_data(CRC_matrix2, "Intratumoral")

crc1$Abundance <- crc1$Abundance / sum(crc1$Abundance)
crc2$Abundance <- crc2$Abundance / sum(crc2$Abundance)


all_top_taxa <- unique(c(
  unique(crc1$Taxon[crc1$Taxon != "Others"]),
  unique(crc2$Taxon[crc2$Taxon != "Others"])
))


combined_abundance <- bind_rows(crc1, crc2) %>%
  filter(Taxon != "Others") %>%
  group_by(Taxon) %>%
  summarise(Total = sum(Abundance)) %>%
  arrange(desc(Total))

col <- rev(c(
      "#9f9f9f", "#55838b", "#e0e096","#8dd3c7","#7ab3b2","#d1756d",
      "#a35e6f", "#c5af6f", "#7a906d","#FD8D3C","#3b8b94","#539faa",
      "#8c96c6","#ccb9d2", "#c9aea4","#9ECAE1", "#f7fcb9", "#43c1a2",
      "#A1D99B","#FED976","#ffb5c1", "#7fc7d4", "#7FC97F", "#FDC086"
    ))


taxa_colors <- setNames(
  col[1:length(all_top_taxa)],
  combined_abundance$Taxon
)
taxa_colors["Others"] <- "#9f9f9f" 


taxon_levels <- c(combined_abundance$Taxon, "Others")

crc1 <- crc1 %>%
  mutate(Taxon = factor(Taxon, levels = taxon_levels)) %>%
  arrange(Taxon)

crc2 <- crc2 %>%
  mutate(Taxon = factor(Taxon, levels = taxon_levels)) %>%
  arrange(Taxon)

create_plot <- function(data, group_name) {
  ggplot(data, aes(x = group_name, y = Abundance, fill = Taxon)) +
    geom_bar(stat = "identity", position = "fill") +
    scale_fill_manual(values = taxa_colors, drop = FALSE) +
    scale_y_continuous(labels = scales::percent_format()) +
    labs(
      x = "",
      y = "Relative Abundance",
      fill = "Taxon"
    ) +
    theme_classic(base_size = 12) +
    theme(
      axis.text.x = element_text(),
      legend.position = "right"
    )
}

plot_crc1 <- create_plot(crc1, "Gut") + 
  theme(legend.position = "none")

plot_crc2 <- create_plot(crc2, "Intratumoral")


library(patchwork)
combined_plot <- plot_crc1 | plot_crc2

pdf("figs/figS1/stacked_bar_plot.pdf", width = 7, height = 4)
print(combined_plot + plot_layout(widths = c(1, 1.2)))
dev.off()
