library(ggplot2)
library(vegan)
library(ggpubr)
library(dplyr)
library(rentrez)
library(MMUPHin)
source("scripts/methods.R")

{
genus_summed <- readRDS("adata2025/data/genus_summed.rds")
genus_summed <- genus_summed[which(substr(rownames(genus_summed), 1, 1) %in% c("N", "T")), ]

genus_summed <- log(genus_summed + 1)
genus_summed <- genus_summed / rowSums(genus_summed)
genus_summed <- genus_summed[, which(colSums(genus_summed > 0) > 0.05*nrow(genus_summed))]

feature_summed <- genus_summed
response <- substr(rownames(feature_summed), 1, 1)
response = case_when(
    response == "N" ~ "Normal",
    response == "T" ~ "CRC"
  )
response <- factor(response, levels = c("Normal", "CRC"))
}

load("/data/baining/ModelAssessment/micobio/feat_meta/CRC_Metagenomics_genus.RData")
feat_list <- my_data0$feat_list
meta_list <- my_data0$meta_list

features_feat <- c()
for (i in 1:7) {
  feat <- feat_list[[i]]
  features_feat <- union(features_feat, colnames(feat))
}

ids <- sub("ncbi_", "", features_feat)
names_feat <- readRDS("adata2025/data/names_feat_genus.rds")
id_to_name <- unlist(names_feat)
features_feat <- as.character(id_to_name[ids])

feature1 <- colnames(genus_summed)
cleaned_string <- sub("_unclassified$", "", feature1)
feature1 <- feature1[which(feature1 %in% cleaned_string)]
feature_summed <- feature_summed[, feature1]
feature1 <- colnames(feature_summed)

all_features <- intersect(feature1, features_feat)
all_features <- union(feature1, na.omit(features_feat)) #Fig S1


for (i in 1:7) {
    feati <- feat_list[[i]]
    metai <- meta_list[[i]]
    id <- sub("ncbi_", "", colnames(feati))
    name <- as.character(id_to_name[id])
    colnames(feati) <- name
    lack_features <- setdiff(all_features, colnames(feati))
    lack_feat <- matrix(0, nrow = nrow(feati), ncol = length(lack_features))
    colnames(lack_feat) <- lack_features
    feati <- cbind(feati, lack_feat)
    feati <- feati[, all_features]

    Y <- as.factor(ifelse(metai$Group == "Case", "CRC", "Normal"))
    meta_response <- factor(Y, levels = c("Normal", "CRC"))
    
    if (i == 1) {
        feat_all <- feati
        response_all <- meta_response
    }else {
        feat_all <- rbind(feat_all, feati)
        response_all <- c(response_all, meta_response)
    }
}
table(response_all)

feat_all <- feat_all / rowSums(feat_all)

lack_features <- setdiff(all_features, colnames(feature_summed))
lack_feat <- matrix(0, nrow = nrow(feature_summed), ncol = length(lack_features))
colnames(lack_feat) <- lack_features
feature_summed <- cbind(feature_summed, lack_feat)
feature_summed <- feature_summed[, all_features]

CRC_matrix1 <- feat_all[response_all == "CRC", ]
CRC_matrix2 <- feature_summed[response == "CRC", ]
CRC_matrix <- rbind(CRC_matrix1, CRC_matrix2)
CRC_response <- c(rep("Fecal", nrow(CRC_matrix1)), rep("Intratumoral", nrow(CRC_matrix2)))
meta2 <- data.frame(studyID = factor(CRC_response))
rownames(meta2) <- rownames(CRC_matrix)
CRC_abd <- t(CRC_matrix)
CRC_abd_adj <- adjust_batch(feature_abd = CRC_abd, batch = "studyID", data = meta2)$feature_abd_ad
CRC_matrix <- t(CRC_abd_adj[which(rowSums(CRC_abd_adj) != 0), ])


pcoa_df <- function(CRC_matrix, response) {
    otu <- CRC_matrix
     otu.distance <- vegdist(otu)
     pcoa <- cmdscale(otu.distance,eig=TRUE)
     pc12 <- pcoa$points[,1:2]
     pc <- round(pcoa$eig/sum(pcoa$eig)*100,digits=2)
     pc12 <- as.data.frame(pc12)
     df <- cbind(pc12, data.frame(Group = response))
     p1 <- ggplot(data=df,aes(x=V1,y=V2,color=Group))+
     geom_point(size=1.5) +
     theme_classic(base_size = 14) +
     scale_color_manual(values = c("#1B9E77", "#D95F02")) +
     labs(x=paste0("PCoA1 (",pc[1],"%)"),
          y=paste0("PCoA2 (",pc[2],"%)"))+
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

pdf(paste0("adata2025/figs/figS1/", "pcoa_CRC.pdf"), width = 4, height = 4)
    print(p1)
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


#FigS1 D
library(phyloseq)
library(ANCOMBC)

my_ancom <- function(feature_summed, response) {
    otu_mat <- as.matrix(feature_summed)
    OTU <- otu_table(otu_mat, taxa_are_rows = FALSE)
    samp_data <- data.frame(group = response, row.names = rownames(feature_summed))
    SAMP <- sample_data(samp_data)

    tax_mat <- matrix(colnames(feature_summed), ncol = 1)
    rownames(tax_mat) <- colnames(feature_summed)
    colnames(tax_mat) <- "Genus"
    TAX <- tax_table(tax_mat)

    # 创建phyloseq对象
    pseq <- phyloseq(OTU, TAX, SAMP)

    out <- ancombc2(
        data = pseq,
        tax_level = "Genus",
        p_adj_method = "fdr",    
        fix_formula = "group",
        prv_cut = 0.05,
        group = "group",
        struc_zero = TRUE,
        neg_lb = TRUE,
        alpha = 0.05,
        global = FALSE
    )

    res <- out$res
    
    diff_abn <- res$diff_groupCRC
    p_val <- res$p_groupCRC
    q_val <- res$q_groupCRC

    results_df <- data.frame(
        microbe = res$taxon,
        p_val = as.numeric(p_val),
        q_val = as.numeric(q_val),
        diff_abn = as.logical(diff_abn)
    )
    return(results_df)
}

response <- factor(response, levels = c("Normal", "CRC"))
ancom_res_tissue <- my_ancom(feature_summed, response)
ancom_res_feat <- my_ancom(feat_all, response_all)

combine_res <- data.frame(microbe = ancom_res_feat$microbe, feat_p = ancom_res_feat$p_val, feat_q = ancom_res_feat$q_val,
                          tissue_p = 0, tissue_q = 0)

for (i in 1:nrow(ancom_res_feat)) {
    microbe <- ancom_res_feat$microbe[i]
    idx <- which(ancom_res_tissue$microbe == microbe)
    if(length(idx) == 0) {
      next
    }
    combine_res$tissue_p[i] <- ancom_res_tissue$p_val[idx]
    combine_res$tissue_q[i] <- ancom_res_tissue$q_val[idx]
}

library(ggplot2)
library(ggrepel)


highlight_data <- combine_res[combine_res$microbe %in% c("Fusobacterium", "Pseudomonas", "Streptococcus"), ]
pdf("figs/figS1/ancom_feat_tissue.pdf", width = 4, height = 3)
ggplot(combine_res, aes(x = log10(feat_q), y = log10(tissue_q))) +
    geom_point(alpha = 0.5) +
    theme_classic(base_size = 14) +
    labs(x = "log10 P-value (Fecal)", y = "log10 P-value (Tissue)") +
    geom_hline(yintercept = log10(0.05), linetype = "dashed", color = "red") +
    geom_vline(xintercept = log10(0.05), linetype = "dashed", color = "red") +
    theme(legend.position = "top") +
    geom_text_repel(
        data = highlight_data,
        aes(label = microbe),
        color = "red",
        segment.color = "red",
        segment.size = 0.8,
        arrow = arrow(length = unit(0.02, "npc"), type = "closed"),
        box.padding = 0.8,
        point.padding = 0.5,
        min.segment.length = 0,
        force = 5,
        max.iter = 3000,
        bg.color = "white",
        bg.r = 0.15
    )
dev.off()
