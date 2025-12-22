library(dplyr)
library(vegan)
library(ggplot2)
library(microeco)
library(magrittr)
library(ggtree)
source("scripts/methods.R")
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

  res <- deal_compares_tr(c("Normal", "CRC"), feature_summed, response)

  feature_summed <- res$feature
  response <- res$response
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
    geom_point(size=1.2) +
    stat_ellipse(level = 0.95, linetype = 2)  +
               theme(legend.position = "top") +
    labs(x=paste0("PCoA1 (",pc[1],"%)"),
        y=paste0("PCoA2 (",pc[2],"%)"))+
        theme_classic() +
     scale_color_manual(values = c("#159415", "#ff8000"))

pdf(paste0("figs/fig1/", "pcoa_fig1.pdf"), width = 5, height = 3.5)
    print(p1)
dev.off()

#DCA
dca_result <- decorana(feature_summed)
dca_points <- scores(dca_result, display = "sites")
dca_points <- data.frame(dca_points[, 1:2])
colnames(dca_points) <- c("DCA1", "DCA2")
dca_points$Group <- response


pdf(paste0("figs/fig1/", "dca_fig1.pdf"), width = 5, height = 3.5)
ggplot(dca_points, aes(x = DCA1, y = DCA2, color = Group)) +
  geom_point(size = 1.2) +
  stat_ellipse(level = 0.95, linetype = 2) +  # 添加95%置信椭圆
  labs(x = paste0("DCA1 (", round(dca_result$evals[1]/sum(dca_result$evals)*100, 1), "%)"),
       y = paste0("DCA2 (", round(dca_result$evals[2]/sum(dca_result$evals)*100, 1), "%)")) +
  theme_classic() +
  scale_color_manual(values = c("#159415", "#ff8000"))
dev.off()

#lefse
sample_table <- NULL

sample_table$Sample_ID <- rownames(feature_summed)
sample_table$group <- response

bio_cat <- readRDS("data/bio_cat.rds")
bio_cnt <- readRDS("data/bio_cnt.rds")

bio_cnt_numeric <- bio_cnt %>% 
  mutate(across(everything(), as.numeric))

bio_cnt_numeric <- as.data.frame(bio_cnt_numeric)
rownames(bio_cnt_numeric) <- rownames(bio_cnt)

bio_cnt_numeric <- bio_cnt_numeric[, rownames(feature_summed)]

head(colnames(bio_cnt_numeric))

sample_table$Sample_ID <- colnames(bio_cnt_numeric)
sample_table$group <- response

sample_table <- as.data.frame(sample_table)
sample_table$group <- as.factor(sample_table$group)
sample_table$Sample_ID <- as.factor(sample_table$Sample_ID)
rownames(sample_table) <- colnames(bio_cnt_numeric)

for (j in 1:ncol(bio_cat)) {
    first_letter <- substr(colnames(bio_cat)[j], 1, 1)
    first_letter <- tolower(first_letter)
    bio_cat[, j] <- paste0(first_letter, "_", bio_cat[, j])
}


mydataset <- microtable$new(sample_table = sample_table,
                          otu_table = bio_cnt_numeric,
                          tax_table = bio_cat)

lefse <- trans_diff$new(dataset = mydataset, 
                        method = "lefse", 
                        group = "group", 
                        filter_thres = 0.05,
                        alpha = 0.25,
                        lefse_subgroup = NULL)

head(lefse$res_diff)


pdf("figs/fig1/lefse_test.pdf", height = 5, width = 6)
lefse$plot_diff_bar(use_number = 1:20,
                      width = 0.8)
dev.off()


abund <- lefse$abund_table
res_diff <- lefse$res_diff
crc_abund <- abund[res_diff$Taxa[1:20], sample_table$group == "CRC"]
normal_abund <- abund[res_diff$Taxa[1:20], sample_table$group == "Normal"]

crc_abund <- crc_abund / lefse$lefse_norm
normal_abund <- normal_abund / lefse$lefse_norm

df_crc <- as.data.frame(t(crc_abund)) 
df_crc$Group <- "CRC"
df_crc$SampleID <- rownames(df_crc)

df_normal <- as.data.frame(t(normal_abund))
df_normal$Group <- "Normal"
df_normal$SampleID <- rownames(df_normal)

df_total <- rbind(df_crc, df_normal)

plot_data <- melt(df_total, 
                  id.vars = c("Group", "SampleID"),
                  variable.name = "Taxa_Long",
                  value.name = "Abundance")

plot_data$Taxa <- sub(".*\\|", "", plot_data$Taxa_Long)

raw_order <- rownames(crc_abund)
clean_order <- sub(".*\\|", "", raw_order)
plot_data$Taxa <- factor(plot_data$Taxa, levels = rev(clean_order))


p_final <- ggplot(plot_data, aes(x = Taxa, y = Abundance, fill = Group)) +
  geom_boxplot(
  outlier.size = 0.3,
               alpha = 0.6, 
               width = 0.8,
               lwd = 0.4) +
  coord_flip() + 
  scale_fill_manual(values = c("Normal" = "#00A087", "CRC" = "#E64B35")) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = "", y = "Relative Abundance") +
  theme_classic() +
  theme(
    panel.grid.major.x = element_line(color = "grey90", linewidth = 0.5), 
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.5),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(size = 10, color = "black"),
    axis.text.x = element_text(size = 10, color = "black"),
    axis.title.x = element_text(size = 12, color = "black"),

    legend.position = "right",
    legend.title = element_text(size = 11, face = "bold"),
    legend.text = element_text(size = 10)
  )

pdf("figs/fig1/lefse_abund_boxplot.pdf", height = 5, width = 6)
print(p_final)
dev.off()
