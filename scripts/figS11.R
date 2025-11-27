library(metagenomeSeq)
library(dplyr)
library(ggplot2)
library(tidyr)
library(purrr)
fish_data <- read.csv("data/fish_data.csv")
genus_summed <- readRDS("data/genus_summed.rds")
genus_summed <- genus_summed[which(substr(rownames(genus_summed), 1, 1) %in% c("N", "a", "T")), ]
genus_summed <- genus_summed[, which(colSums(genus_summed > 0) > 0.05*nrow(genus_summed))]

ms_obj <- newMRexperiment(t(genus_summed))
ms_obj <- cumNorm(ms_obj, p = 0.5)
genus_summed <- t(MRcounts(ms_obj, norm = TRUE, log = FALSE))

response <- substr(rownames(genus_summed), 1, 1)
response = dplyr::case_when(
    response == "N" ~ "Normal",
    response == "a" ~ "Polyp",
    response == "T" ~ "CRC"
  )

cared_feature <- intersect(colnames(fish_data), colnames(genus_summed))
bio_data <- data.frame(genus_summed[, cared_feature])

bio_data$Tissue <- response

fish_data$Tissue[fish_data$Tissue == "Tumor"] <- "CRC"


fish_summary <- fish_data %>%
  pivot_longer(
    cols = -Tissue,
    names_to = "Bacterium",
    values_to = "Value"
  ) %>%
  mutate(Method = "Fluorescence")

bio_summary <- bio_data %>%
  pivot_longer(
    cols = -Tissue,
    names_to = "Bacterium",
    values_to = "Value"
  ) %>%
  mutate(Method = "Abundance")

combined_summary <- bind_rows(fish_summary, bio_summary) %>%
  group_by(Bacterium, Method) %>%
  mutate(ScaledValue = (Value - mean(Value)) / sd(Value)) %>%
  ungroup()


plot_data <- combined_summary %>%
  group_by(Tissue, Bacterium, Method) %>%
  summarise(
    Mean_Scaled = mean(ScaledValue),
    SE_Scaled = sd(ScaledValue) / sqrt(n()), # Standard Error
    .groups = 'drop'
  )

plot_data$Tissue <- factor(plot_data$Tissue, levels = c("Normal", "Polyp", "CRC"))

trend_plot <- ggplot(
  plot_data, 
  aes(x = Tissue, y = Mean_Scaled, group = Method, color = Method)
) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  geom_errorbar(
    aes(ymin = Mean_Scaled - SE_Scaled, ymax = Mean_Scaled + SE_Scaled),
    width = 0.1
  ) +
  facet_wrap(~ Bacterium, scales = "free_y") +
  labs(
    title = "Comparison of Trends Across Tissues",
    subtitle = "Data are Z-score scaled to compare relative patterns",
    x = "Tissue",
    y = "Mean Scaled Value",
    color = "Measurement Method"
  ) +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5),
    legend.position = "top"
  )
ggsave("figs/figS12/trend_plot.pdf", plot = trend_plot, width = 8, height = 4, dpi = 300)


point_plot_data <- data.frame(
  Fluorescence = rep(0, 9),
  Abundance = rep(0, 9),
  Bacterium = rep("a", 9),
  Tissue = rep(c("Normal", "Polyp", "CRC"), 3)
)

for (i in 1:9) {
  point_plot_data$Abundance[i] <- plot_data$Mean_Scaled[2*i - 1]
  point_plot_data$Fluorescence[i] <- plot_data$Mean_Scaled[2*i]
  point_plot_data$Bacterium[i] <- plot_data$Bacterium[2*i]
  point_plot_data$Tissue[i] <- plot_data$Tissue[2*i]
}
point_plot_data$Tissue <- rep(c("CRC", "Polyp", "Normal"), 3)
point_plot_data$Tissue <- factor(point_plot_data$Tissue, levels = c("Normal", "Polyp", "CRC"))
pdf("figs/figS12/point_plot.pdf", width = 4, height = 3)
point_plot <- ggplot(point_plot_data, aes(x = Abundance, y = Fluorescence, color = Tissue)) +
  geom_point(size = 3) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey") +
  labs(
    title = "Mean Scaled Values",
    x = "Abundance",
    y = "Fluorescence",
    color = "Tissue Type"
  ) +
  scale_x_continuous(limits = c(-1.2, 1.5)) +
  scale_y_continuous(limits = c(-1.2, 1.5)) +
  theme_classic()
print(point_plot)
dev.off()

cor.test(point_plot_data$Abundance, point_plot_data$Fluorescence)
