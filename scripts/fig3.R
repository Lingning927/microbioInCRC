library(dplyr)
library(vegan)
library(ggplot2)
library(ggsci)
library(tidyr)

#load the data
{
  genus_summed <- readRDS("data/genus_summed.rds")
  genus_summed <- genus_summed[which(substr(rownames(genus_summed), 1, 1) %in% c("N", "a", "T")), ]
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

#PCoA
otu <- feature_summed
otu.distance <- vegdist(otu)
pcoa <- cmdscale(otu.distance,eig=TRUE)
pc12 <- pcoa$points[,1:2]
pc <- round(pcoa$eig/sum(pcoa$eig)*100,digits=2)
pc12 <- as.data.frame(pc12)
df <- cbind(pc12, data.frame(Group = response))
df_anosim <- anosim(otu.distance, df$Group, permutations = 999)

p1 <- ggplot(data=df,aes(x=V1,y=V2,color=Group))+
    geom_point(size=1.2) +#绘制点图并设定大小
    stat_ellipse(level = 0.95, linetype = 2)  +
               theme(legend.position = "top") +
    labs(x=paste0("PCoA1 (",pc[1],"%)"),
        y=paste0("PCoA2 (",pc[2],"%)"))+#将x、y轴标题改为贡献度
        theme_classic() +
     scale_color_d3()

pdf(paste0("figs/fig3/", "pcoa_fig3.pdf"), width = 5, height = 3.5)
    print(p1)
dev.off()

#DCA
dca_result <- decorana(feature_summed)
dca_points <- scores(dca_result, display = "sites")
dca_points <- data.frame(dca_points[, 1:2])
colnames(dca_points) <- c("DCA1", "DCA2")
dca_points$Group <- response

pdf(paste0("figs/fig3/", "dca_fig3.pdf"), width = 5, height = 3.5)
ggplot(dca_points, aes(x = DCA1, y = DCA2, color = Group)) +
  geom_point(size = 1.2) +
  stat_ellipse(level = 0.95, linetype = 2) +  # 添加95%置信椭圆
  labs(x = paste0("DCA1 (", round(dca_result$evals[1]/sum(dca_result$evals)*100, 1), "%)"),
       y = paste0("DCA2 (", round(dca_result$evals[2]/sum(dca_result$evals)*100, 1), "%)")) +
  theme_classic() +
  scale_color_d3()
dev.off()

#Stack bar plot
genus_summed <- readRDS("data/genus_summed.rds")
family_summed <- readRDS("data/family_summed.rds")

genus_summed <- log(genus_summed + 1)
family_summed <- log(family_summed + 1)

construct_stack_df <- function(genus_summed) {
    top_genus <- colnames(genus_summed)[order(colSums(genus_summed[substr(rownames(genus_summed), 1, 1) == "T", ]), decreasing = TRUE)[1:19]]
    top_genus_mat <- cbind(genus_summed[, top_genus], rowSums(genus_summed[, -which(colnames(genus_summed) %in% top_genus)]))
    colnames(top_genus_mat)[20] <- "Others"

    df <- as.data.frame(top_genus_mat)
    df$Genus <- rownames(df)

    group_means <- df %>%
    mutate(Group = substr(Genus, 1, 1)) %>%
    group_by(Group) %>%
    summarise(across(-Genus, mean, na.rm = TRUE))
    group_means_normalized <- group_means %>%
    mutate(Total = rowSums(select(., -Group))) %>%
    pivot_longer(-c(Group, Total), names_to = "Variable", values_to = "Mean") %>%
    mutate(Proportion = Mean / Total) %>%
    select(-Mean, -Total)

    level <- colnames(top_genus_mat)[order(group_means_normalized[group_means_normalized$Group == "T", 3])]
    level[1:19] <- rev(level[1:19])

    group_means_normalized$Variable <- factor(group_means_normalized$Variable,
        levels = colnames(top_genus_mat))

    group_means_normalized <- group_means_normalized %>%
    mutate(Group = case_when(
        Group == "C" ~ "Control",
        Group == "E" ~ "Control",
        Group == "F" ~ "Feacal",
        Group == "N" ~ "Normal",
        Group == "a" ~ "Polyp",
        Group == "T" ~ "CRC"
    ))
    group_means_normalized$Group <- factor(group_means_normalized$Group, levels = c("Control", "Feacal", "Normal", "Polyp", "CRC"))
    group_means_normalized <- group_means_normalized %>% filter(Group %in% c("Normal", "Polyp", "CRC"))
    group_means_normalized$Group <- factor(group_means_normalized$Group, levels = c("Normal", "Polyp", "CRC"))
    return(group_means_normalized)
}

col <- rev(c(
        "#9f9f9f", "#e0e096","#8dd3c7","#7ab3b2","#d1756d",
        "#a35e6f","#7a906d","#FD8D3C","#3b8b94","#539faa",
        "#8c96c6","#ccb9d2","#c9aea4","#9ECAE1",
        "#A1D99B","#FED976","#ffb5c1", "#7fc7d4", "#7FC97F", "#FDC086"
))

group_means_normalized <- construct_stack_df(genus_summed)
pdf("figs/fig3/stack_bar_genus.pdf", width = 7, height = 3.5)
ggplot(group_means_normalized, aes(x = Group, y = Proportion, fill = Variable)) +
  geom_bar(stat = "identity") +
  theme_classic(base_size = 12) +
  scale_fill_manual(values = col) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(y = "Proportion", x = "", fill = "Genus") +
  guides(fill = guide_legend(ncol = 2)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5))
dev.off()

group_means_normalized <- construct_stack_df(family_summed)
pdf("figs/fig3/stack_bar_family.pdf", width = 5.5, height = 3.5)
ggplot(group_means_normalized, aes(x = Group, y = Proportion, fill = Variable)) +
  geom_bar(stat = "identity") +
  theme_classic(base_size = 12) +
  scale_fill_manual(values = col) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(y = "Proportion", x = "", fill = "Genus") +
  guides(fill = guide_legend(ncol = 2)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5))
dev.off()


###Fig 3E
library(ggplot2)
library(reshape2)
library(ggplot2)
library(ggpubr)
library(tidyr)
library(dplyr)

genus_summed <- readRDS("data/genus_summed.rds")
genus_summed <- genus_summed[which(substr(rownames(genus_summed), 1, 1) %in% c("N", "a", "T")), ]
genus_summed <- genus_summed[, which(colSums(genus_summed > 0) > 0.05*nrow(genus_summed))]

feature_summed <- genus_summed

for (i in 1 : ncol(feature_summed)) {
  x <- feature_summed[, i]
  quantiles <- as.numeric(quantile(x, probs = c(0.25, 0.5, 0.75)))
  quantiles[1] <- quantiles[1] - 0.01
  quantiles[3] <- quantiles[3] + 0.01
  categorized_x <- cut(x, breaks = c(-Inf, quantiles, Inf), labels = c(1, 2, 3, 4), include.lowest = TRUE)
}

cochran_armitage_results <- apply(feature_summed, 2, function(x) {
  binary_x <- ifelse(x > median(x, na.rm = TRUE), 1, 0)
  result <- tryCatch({
    test_result <- DescTools::CochranArmitageTest(table(response, binary_x))
    test_result$p.value 
  }, error = function(e) {
    return(NA)
  })
  return(result)
})

significant_features_cochran <- cochran_armitage_results[order(cochran_armitage_results)[1:20]]

res_feature <- names(significant_features_cochran)
sig_data <- as.data.frame(feature_summed[, res_feature])
sig_data$Class <- response

non_zero_proportions <- sapply(res_feature, function(feature) {
  tapply(sig_data[[feature]] != 0, sig_data$Class, mean)
})


non_zero_proportions_melted <- melt(non_zero_proportions)
colnames(non_zero_proportions_melted) <- c("Class", "Feature", "NonZeroProportion")

non_zero_proportions_melted$Feature <- paste0("g_", non_zero_proportions_melted$Feature)
non_zero_proportions_melted$Feature <- factor(non_zero_proportions_melted$Feature, levels = paste0("g_", res_feature))

non_zero_proportions_melted$Class <- factor(non_zero_proportions_melted$Class, levels = rev(levels(non_zero_proportions_melted$Class)))
pdf("figs/fig3/non_zero_ca_genus.pdf", width = 5, height = 6)
ggplot(non_zero_proportions_melted, aes(x = Feature, y = NonZeroProportion, fill = Class)) +
  geom_bar(stat = "identity", position = "dodge") +
  coord_flip() +
  theme_minimal() +
  labs(title = "Non-zero Proportion",
       x = "",
       y = "Non-zero Proportion") +
  scale_fill_brewer(palette = "Set2") 
dev.off()

top_feature_data_melted <- melt(sig_data, id.vars = "Class", variable.name = "Feature")

top_feature_data_melted$Feature <- paste0("g_", top_feature_data_melted$Feature)
top_feature_data_melted$Feature <- factor(top_feature_data_melted$Feature, levels = paste0("g_", res_feature))

top_feature_data_melted$value <- log(top_feature_data_melted$value + 1)

top_feature_data_melted$Class <- factor(top_feature_data_melted$Class, levels = rev(levels(top_10_feature_data_melted$Class)))

pdf("figs/fig3/boxplot_ca_genus.pdf", width = 5, height = 6)
ggplot(top_feature_data_melted, aes(x = Feature, y = value, fill = Class)) +
  geom_boxplot(outlier.size = 0.8) +
  coord_flip() +
  theme_minimal() +
  labs(title = "Boxplot",
       x = "",
       y = "Log Abundance") +
  scale_fill_brewer(palette = "Set2")  # 可以选择不同的调色板
dev.off()

train_df <- as.data.frame(feature_summed)
feature_means <- train_df %>%
    select(all_of(res_feature)) %>%
    mutate(Class = response) %>%
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
  
df2$Feature <- paste0("g_", df2$Feature)
df2$Feature <- factor(df2$Feature, levels = paste0("g_", res_feature))
p2 <- df2 %>%
    ggplot(aes(x = Class, y = Feature, fill = Value)) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(aes(label = format(Value, digits = 1)), color = "black", size = 3) +
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


pdf("figs/fig3/feature_heatmap_ca.pdf", height = 6, width = 6)
print(p2)
dev.off()
