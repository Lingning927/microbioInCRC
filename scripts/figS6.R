library(ggplot2)
library(ggpubr)
library(rstatix)
library(dplyr)
library(ggtext)
source("scripts/methods.R")

# Load the data
{
    patient_info <- readRDS("data/patient_info.rds")
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

intersect_features <- c("Pseudomonas", "Fusobacterium", "Pigmentiphaga")

normal_crc_features <- c("Phascolarctobacterium", "Granulicatella", "Leptotrichia", "Streptococcus")

other_features <- c("Parvimonas", "Collinsella", "Parabacteroides", "Bacteroides")

color_list <- c("red", "red", "purple", rep("orange", 4), "#26a8e5", "#A40B5D", rep("#4fb94f", 2))

problem <- "survival"
problem_list <- c("HD_MD", "Location", "age", "gender", "survival", "I.II_III.IV")

for (problem in problem_list) {
    deal_res <- deal_problem(feature_summed, patient_info, problem)

    response <- deal_res$response
    predictors <- deal_res$feature

    cared_data <- data.frame(Genus = character(),
                            Abundance = numeric(),
                            Group = character())
for (genus in c(intersect_features, normal_crc_features, other_features)) {
    tmp_df <- data.frame(Genus = genus,
                            Abundance = log(predictors[, genus] + 1),
                            Group = response)
        cared_data <- rbind(cared_data, tmp_df)
    }
    cared_data$Genus <- factor(cared_data$Genus, levels = c(intersect_features, normal_crc_features, other_features))

    stat.test <- cared_data %>% 
    group_by(Genus) %>%
    wilcox_test(Abundance ~ Group) %>%
    adjust_pvalue(method = "BH") %>%
    add_significance("p.adj") %>%
    add_xy_position(x = "Genus",  step.increase = 0.15)  

    pdf(paste0("figs/figS6/boxplot_", problem, ".pdf"), width = 6.5, height = 3)
    p <- ggboxplot(cared_data, x = "Genus", y = "Abundance", fill = "Group",
         palette = "lancet", width = 0.8, outlier.shape = 20) +
    stat_pvalue_manual(
        stat.test, label = "p.adj.signif",
        step.increase = 0, hide.ns = FALSE, tip.length = 0.01,
        bracket.size = 1, size = 4, color = "black"
        ) +
    scale_fill_d3() +
    scale_y_continuous(limits = c(0, max(cared_data$Abundance) * 1.1)) +
    labs(x = "", y = "Abundance") +
    theme_classic() +
    theme(
        axis.text.x = element_markdown(
      angle = 45, 
      hjust = 1,
      color = color_list
    )
    )
    print(p)
    dev.off()
}
