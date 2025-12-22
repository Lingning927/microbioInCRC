library(ggplot2)
library(vegan)
library(dplyr)
library(ggsci)
library(ggrepel)
library(stringr)
library(glmnet)
library(dplyr)
library(tidyr)
library(patchwork)

get_group_label <- function(patient_info, problem) {
  n <- dim(patient_info)[2]
  group_label <- rep("-", n)
  if (problem == "age") {
    ages <- as.numeric(patient_info[problem, ])
    for (i in 1:n) {
      if(ages[i] <= 60) {
        group_label[i] <- "Young (<=60)"
      }else {
        group_label[i] <- "Senior (>60)"
      }
    }
  }else if (problem == "gender") {
    group_label <- as.character(patient_info["gender", ])
  }else if (problem == "Location") {
    location <- as.character(patient_info["Location", ])
    location[which(location == "left")] <- "Left"
    location[which(location == "right")] <- "Right"
    group_label <- location
    #group_label[!(group_label %in% c("Left", "Rectum", "Right"))] <- "-"
  }else if (problem == "Location2") {
    location <- as.character(patient_info["Location", ])
    location[which(location == "left")] <- "Left"
    location[which(location == "right")] <- "Right"
    group_label <- location
    group_label[!(group_label %in% c("Right"))] <- "Left"
  }else if (problem == "survival") {
    survival <- as.character(patient_info[problem, ])
    group_label[which(survival %in% c("w", "W"))] <- "Dead"
    group_label[which(survival %in% c("c", "C"))] <- "Alive"
  }else if (problem == "M") {
    group_label <- as.character(patient_info[problem, ])
    group_label[group_label == "MO"] <- "M0"
  }else if (problem == "N") {
    group_label <- as.character(patient_info[problem, ])
    #group_label[(group_label %in% c("N1", "N2"))] <- "N1+N2"
  }else if (problem == "T") {
    group_label <- as.character(patient_info[problem, ])
    #group_label[(group_label %in% c("T1", "T2"))] <- "T1+T2"
    #group_label[(group_label %in% c("T3", "T4"))] <- "T3+T4"
  }else if (problem == "MSI") {
    group_label <- as.character(patient_info[problem, ])
    group_label[group_label == "MSI"] <- "MSI-H"
    group_label[group_label == "MSI-L"] <- "MSI-H"
  }else if (problem == "CEA(ng/mL)") {
    cea <- as.character(patient_info[problem, ])
    cea[cea == "-"] <- -1
    #cea[495] <- 1.87
    cea <- as.numeric(cea)
    for(i in 1:n) {
      if(cea[i] < 0) {
        group_label[i] <- "-"
      }else if(cea[i] <= 5) {
        group_label[i] <- "Normal"
      }else {
        group_label[i] <- "Higher"
      }
    }
  }else if (problem == "relapse after operation(month)") {
    cea <- as.character(patient_info[problem, ])
    cea[cea == "-"] <- -1
    cea <- as.numeric(cea)
    for(i in 1:n) {
      if(cea[i] < 0) {
        group_label[i] <- "-"
      }else if(cea[i] <= 12) {
        group_label[i] <- "One-"
      }else if(cea[i] <= 36) {
        group_label[i] <- "Two-Three"
      }else {
        group_label[i] <- "Three+"
      }
    }
  }else if (problem == "albumin(g/L)") {
    cea <- as.character(patient_info[problem, ])
    cea[cea == "24,5"] <- 24.5
    cea[cea == "31.5`"] <- 31.5
    cea[547] <- 38.1
    cea <- as.numeric(cea)
    for(i in 1:n) {
      if(cea[i] >= 35) {
        group_label[i] <- "Normal"
      }else {
        group_label[i] <- "Lower"
      }
    }
  }else if (problem == "PreAndPost") {
    ms <- get_group_label(patient_info, "metastatic site")
    pms <- get_group_label(patient_info, "postoperative metastatic site")
    id1 <- names(ms[(ms != "-") & (pms == "-")])
    id2 <- names(ms[(ms == "-") & (pms != "-")])
    group_label <- c(rep("Preoperative", length(id1)), rep("Postoperative", length(id2)))
    names(group_label) <- c(id1, id2)
  }else if (problem == "metastasis") {
    group_label <- get_group_label(patient_info, "metastatic site")
    group_label <- ifelse(group_label == "-", "Non-metastasis",
      "Metastasis")
  }else if (problem == "I.II_III.IV") {
    group_label <- get_group_label(patient_info, "I.II.III.IV")
    group_label <- ifelse(group_label %in% c("I", "II"), "I,II",
      "III,IV")
  }else if (problem == "I_II.III.IV") {
    group_label <- get_group_label(patient_info, "I.II.III.IV")
    group_label <- ifelse(group_label == "I", "I",
      "II,III,IV")
  }else if (problem == "I.II.III_IV") {
    group_label <- get_group_label(patient_info, "I.II.III.IV")
    group_label <- ifelse(group_label == "IV", "IV",
      "I,II,III")
  }else if (problem == "I.II_III_IV") {
    group_label <- get_group_label(patient_info, "I.II.III.IV")
    group_label <- group_label[(group_label %in% c("I", "II"))] <- "I,II"
  }else {
    group_label <- as.character(patient_info[problem, ])
  }
  if(problem != "PreAndPost") {
      names(group_label) <- colnames(patient_info)
  }
  return(group_label)
}

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

info_tb <- read.csv("data/info_tb.csv")
rownames(info_tb) <- info_tb$X
info_tb <- info_tb[, -1]
survival_state <- get_group_label(info_tb, "survival")
age <- get_group_label(info_tb, "age")
gender <- get_group_label(info_tb, "gender")
stage <- get_group_label(info_tb, "Stage")
location <- get_group_label(info_tb, "Location2")
differentiation <- get_group_label(info_tb, "differentiation")
feature <- deal_res$feature
chemotherapy <- get_group_label(info_tb, "chemotherapy")
id_survival <- which(survival_state != "-")
id_zhiliao <- which(chemotherapy != "" & chemotherapy != "-")

chemotherapy[] <- "No"
chemotherapy[id_zhiliao] <- "Yes"
table(chemotherapy)
survival_state <- survival_state[id_survival]

survival_state <- factor(survival_state, levels = c("Alive", "Dead"))
chemotherapy <- factor(chemotherapy)
stage <- factor(stage, levels = c("I", "II", "III", "IV"))
chemotherapy <- chemotherapy[id_survival]
stage <- stage[id_survival]
age <- age[id_survival]
gender <- gender[id_survival]
location <- location[id_survival]
differentiation <- differentiation[id_survival]

meta_data <- data.frame(
  survival_state = survival_state,
  stage = stage,
  chemotherapy = chemotherapy,
  age = age,
  gender = gender,
  location = location,
  differentiation = differentiation
)

microbiota_data <- feature[rownames(meta_data), ]
all(rownames(meta_data) == rownames(microbiota_data))


#Figure S4 A
results <- list()
bacterium_names <- colnames(microbiota_data)

for (bacterium in bacterium_names) {
  temp_df <- cbind(meta_data, bacterium_abundance = microbiota_data[, bacterium])
  model_fit <- try(
    glm(survival_state ~ stage + chemotherapy + age + gender + location + differentiation + bacterium_abundance,
        data = temp_df, family = binomial),
    silent = TRUE
  )
    p_val <- summary(model_fit)$coefficients["bacterium_abundance", "Pr(>|z|)"]
    coef <- summary(model_fit)$coefficients["bacterium_abundance", "Estimate"]
    results[[bacterium]] <- c(p_value = p_val, odds_ratio = exp(coef))
}

results_df <- do.call(rbind, results)
results_df <- as.data.frame(results_df)
results_df$bacterium <- rownames(results_df)

results_df$p_adj <- p.adjust(results_df$p_value, method = "BH")
results_df$log2_odds_ratio <- log2(results_df$odds_ratio)
results_df$neg_log10_p_adj <- -log10(results_df$p_adj)

pdf("figs/figS4/volcano_survival.pdf", width = 4, height = 3)
ggplot(results_df, aes(x = log2_odds_ratio, y = neg_log10_p_adj)) +
  geom_point(aes(color = p_adj < 0.05), alpha = 0.6, size = 2) +
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "grey"), 
                     name = "FDR < 0.05",
                     labels = c("Significant", "Not Significant")) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue") +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "grey") +
  labs(
    x = "Log2 (Odds Ratio)",
    y = "-Log10 (Adjusted P-value)"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom", plot.title = element_text(hjust = 0.5, face="bold")) +
  geom_text_repel(data = subset(results_df, p_adj < 0.05),
                  aes(label = bacterium), size = 3)
dev.off()


#Fig S4B
X <- model.matrix(survival_state ~ stage + chemotherapy + age + gender + location + differentiation, data = meta_data)[, -1] 
X <- cbind(X, microbiota_data)
Y <- meta_data$survival_state

set.seed(123)
cv_lasso <- cv.glmnet(X, Y, family = "binomial", alpha = 1) # alpha=1 代表LASSO
best_lambda <- cv_lasso$lambda.min
coefs <- coef(cv_lasso, s = best_lambda)
selected_vars_names <- rownames(coefs)[which(coefs != 0)]
selected_vars_names <- selected_vars_names[selected_vars_names != "(Intercept)"]
selected_vars_names[4] <- "ageYoung"

if (length(selected_vars_names) == 0) {
  print("LASSO did not select any variables. Only clinical factors might be relevant if forced.")
} else {
  final_data <- as.data.frame(X)
  final_data$survival_state <- Y
  colnames(final_data)[5] <- "ageYoung"
  final_formula <- as.formula(paste("survival_state ~", paste(selected_vars_names, collapse = " + ")))
  final_glm <- glm(final_formula, data = final_data, family = binomial)
}


glm_summary_df <- as.data.frame(summary(final_glm)$coefficients)
glm_summary_df$variable <- rownames(glm_summary_df)
rownames(glm_summary_df) <- NULL
glm_summary_df <- glm_summary_df %>% filter(variable != "(Intercept)")

round_df <- function(x, digits) {
  format(round(x, digits), nsmall = digits)
}
glm_summary_df <- glm_summary_df[which(glm_summary_df[, 4] < 0.05), ]

plot_data <- glm_summary_df %>%
  mutate(
    odds_ratio = exp(Estimate),
    ci_low = exp(Estimate - 1.96 * `Std. Error`),
    ci_high = exp(Estimate + 1.96 * `Std. Error`),
    estimate_str = round_df(Estimate, 2),
    std_error_str = round_df(`Std. Error`, 2),
    p_value_str = ifelse(`Pr(>|z|)` < 0.001, "< 0.001", round_df(`Pr(>|z|)`, 3)),
    or_ci_str = paste0(round_df(odds_ratio, 2), " (", round_df(ci_low, 2), "-", round_df(ci_high, 2), ")")
  ) %>%
  arrange(odds_ratio) %>%
  mutate(variable_wrapped = str_wrap(variable, width = 10)) %>%
  mutate(variable_wrapped = factor(variable_wrapped, levels = .$variable_wrapped))



p_table <- ggplot(plot_data, aes(y = variable_wrapped)) +
  geom_text(aes(x = 0, label = variable_wrapped), hjust = 0, vjust = 0.5, size = 3.5) +
  geom_text(aes(x = 1.5, label = estimate_str), hjust = 0.5) +
  geom_text(aes(x = 2.5, label = std_error_str), hjust = 0.5) +
  geom_text(aes(x = 3.5, label = p_value_str), hjust = 0.5) +
  geom_text(aes(x = 4.7, label = or_ci_str), hjust = 0.5) +
  annotate("text", x = 0, y = nrow(plot_data) + 1.2, label = "Variable", hjust = 0, fontface = "bold") +
  annotate("text", x = 1.5, y = nrow(plot_data) + 1.2, label = "Estimate", hjust = 0.5, fontface = "bold") +
  annotate("text", x = 2.5, y = nrow(plot_data) + 1.2, label = "Std. Error", hjust = 0.5, fontface = "bold") +
  annotate("text", x = 3.5, y = nrow(plot_data) + 1.2, label = "P-value", hjust = 0.5, fontface = "bold") +
  annotate("text", x = 4.7, y = nrow(plot_data) + 1.2, label = "Odds Ratio (95% CI)", hjust = 0.5, fontface = "bold") +
  
  theme_void() +
  coord_cartesian(xlim = c(0, 5.7))


p_forest <- ggplot(plot_data, aes(x = odds_ratio, y = variable_wrapped)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey50") +
  geom_point(size = 3, color = "#0072B2") +
  geom_errorbarh(aes(xmin = ci_low, xmax = ci_high), height = 0.2, color = "#0072B2") +
  labs(x = "Odds Ratio") +
  theme_classic() +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.y = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major.x = element_line(color = "grey90", linetype = "dotted")
  ) +
  scale_x_log10()


final_plot <- p_table + p_forest + plot_layout(widths = c(3.5, 2))

final_plot_with_title <- final_plot + 
  plot_annotation(title = "Multivariable Logistic Regression Model after LASSO Selection") & 
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

pdf("figs/figS4/forest_plot_survival_combined.pdf", width = 10, height = 6)
print(final_plot_with_title)
dev.off()

#Fig S4 C
full_data <- cbind(meta_data, microbiota_data)
colnames(full_data) <- make.names(colnames(full_data))
set.seed(45)
train_index <- createDataPartition(full_data$survival_state, p = 0.7, list = FALSE)
train_data <- full_data[train_index, ]
test_data <- full_data[-train_index, ]

set.seed(42)
rf_clinical <- randomForest(
  survival_state ~ stage + chemotherapy,
  data = train_data,
  ntree = 500,
  importance = TRUE
)

set.seed(42)
rf_full <- randomForest(
  survival_state ~ .,
  data = train_data,
  ntree = 500,
  importance = TRUE
)

pred_clinical_prob <- predict(rf_clinical, newdata = test_data, type = "prob")[, "Dead"]
pred_full_prob <- predict(rf_full, newdata = test_data, type = "prob")[, "Dead"]

roc_clinical <- roc(test_data$survival_state, pred_clinical_prob, levels = c("Alive", "Dead"))
roc_full <- roc(test_data$survival_state, pred_full_prob, levels = c("Alive", "Dead"))

auc_clinical_val <- auc(roc_clinical)
auc_full_val <- auc(roc_full)

cat("AUC for Clinical-Only Model:", round(auc_clinical_val, 4), "\n")
cat("AUC for Full Model (Clinical + Microbiota):", round(auc_full_val, 4), "\n")


roc_list <- list(
  "Clinical Only" = roc_clinical,
  "Clinical + Microbiota" = roc_full
)

pdf("figs/figS4/roc_survival.pdf", width = 4, height = 3.5)
ggroc(roc_list, legacy.axes = TRUE, linewidth = 1) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey") +
  labs(
    x = "1 - Specificity",
    y = "Sensitivity",
    color = ""
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    legend.position = "up",
    legend.title = element_text(face = "bold"),
    axis.title = element_text(size = 12)
  ) +
  scale_color_manual(
    name = "Model",
    values = c("Clinical Only" = "#0072B2", "Clinical + Microbiota" = "#D55E00"),
    labels = c(
      `Clinical Only` = paste0("Clinical Only (AUC = ", round(auc_clinical_val, 3), ")"),
      `Clinical + Microbiota` = paste0("Clinical + Microbiota (AUC = ", round(auc_full_val, 3), ")")
    )
  ) +
  annotate("text", x = 0.75, y = 0.25, 
           label = paste("AUC (Clinical Only):", round(auc_clinical_val, 3)), 
           color = "#0072B2", size = 4) +
  annotate("text", x = 0.75, y = 0.18, 
           label = paste("AUC (Full Model):", round(auc_full_val, 3)), 
           color = "#D55E00", size = 4)

dev.off()



survival_time_info <- read.csv("data/Survival_info.csv")
head(as.character(survival_time_info[1, ]))
state2 <- as.character(survival_time_info[1, rownames(microbiota_data)])
meta_data$survival_state <- factor(state2, levels = c("Dead", "Alive"))
time2 <- as.numeric(survival_time_info[2, rownames(microbiota_data)])
time2[state2 == "Alive"] <- 60


meta_data$survival_state <- ifelse(state2 == "Dead", 1, 0)

library(metagenomeSeq)
  ms_obj <- newMRexperiment(t(microbiota_data))
  ms_obj <- cumNorm(ms_obj, p = 0.5)
  otu_table_norm <- t(MRcounts(ms_obj, norm = TRUE, log = TRUE))
  microbiota_data <- otu_table_norm

#####
library(survival)
library(ggplot2)
library(ggrepel)

results <- list()
bacterium_names <- colnames(microbiota_data)

for (bacterium in bacterium_names) {
  temp_df <- cbind(meta_data, bacterium_abundance = microbiota_data[, bacterium])

  clean_df <- temp_df[!is.na(temp_df$survival_time) & 
                      !is.na(temp_df$survival_state) & 
                      !is.na(temp_df$bacterium_abundance), ]
  if (nrow(clean_df) < 10) next 

  model_fit <- try({
    coxph(Surv(survival_time, survival_state) ~ stage + chemotherapy + age + 
            gender + location + differentiation + bacterium_abundance,
          data = clean_df)
  }, silent = TRUE)
  
  if (!inherits(model_fit, "try-error")) {
    s <- summary(model_fit)
    
    if ("bacterium_abundance" %in% rownames(s$coefficients)) {
      p_val <- s$coefficients["bacterium_abundance", "Pr(>|z|)"]
      hr <- s$coefficients["bacterium_abundance", "exp(coef)"]
      coef_val <- s$coefficients["bacterium_abundance", "coef"]
      
      results[[bacterium]] <- c(p_value = p_val, hazard_ratio = hr, coef = coef_val)
    }
  }
}

results_df <- as.data.frame(do.call(rbind, results))
results_df$bacterium <- rownames(results_df)
results_df$p_adj <- p.adjust(results_df$p_value, method = "BH")
results_df$log2_hr <- log2(results_df$hazard_ratio)
results_df$neg_log10_p_adj <- -log10(results_df$p_adj)
sum(results_df$p_adj < 0.05)

pdf("figs/figS3/volcano_survival_cox.pdf", width = 5, height = 4)
ggplot(results_df, aes(x = log2_hr, y = neg_log10_p_adj)) +
  geom_point(aes(color = p_adj < 0.05), alpha = 0.6, size = 2) +
  scale_color_manual(values = c("FALSE" = "grey", "TRUE" = "red"), 
                     name = "FDR < 0.05",
                     labels = c("Not Significant", "Significant")) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue") +
  geom_vline(xintercept = c(-log2(1.5), log2(1.5)), linetype = "dashed", color = "grey") + 
  labs(
    title = "Survival Analysis (Cox PH)",
    x = "Log2 (Hazard Ratio)",
    y = "-Log10 (Adjusted P-value)"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom") +
  geom_text_repel(data = subset(results_df, p_adj < 0.05),
                  aes(label = bacterium), size = 3)
dev.off()

# =======================================================
library(survival)
library(glmnet)
library(survminer)
library(ggplot2)
library(dplyr)
library(patchwork) 


scaled_abundance <- scale(microbiota_data)
full_data <- cbind(meta_data, scaled_abundance)



clean_data <- na.omit(full_data)

bacterium_names <- colnames(microbiota_data)
microbio_clean <- as.matrix(clean_data[, bacterium_names])
meta_clean <- clean_data[, !colnames(clean_data) %in% bacterium_names]



surv_obj <- Surv(meta_clean$survival_time, meta_clean$survival_state)

set.seed(200)
cv_fit <- cv.glmnet(x = microbio_clean, y = surv_obj, family = "cox", alpha = 1)

selected_lambda <- cv_fit$lambda.min
coefs <- coef(cv_fit, s = selected_lambda)

selected_bacteria <- rownames(coefs)[coefs[, 1] != 0]


selected_bacteria <- selected_bacteria[-4]

if (length(selected_bacteria) > 0) {
  covariates <- c("stage", "chemotherapy", "age", "gender", "location", "differentiation")
  final_variables <- c(covariates, selected_bacteria)
  final_variables_safe <- paste0("`", final_variables, "`")

  formula_str <- paste("Surv(survival_time, survival_state) ~", 
                       paste(final_variables_safe, collapse = " + "))
  final_formula <- as.formula(formula_str)
  final_cox_model <- coxph(final_formula, data = clean_data)
  print(summary(final_cox_model))
  
}




cox_summary <- summary(final_cox_model)
cox_summary_df <- as.data.frame(cox_summary$coefficients)

# 清理列名
cox_summary_df$variable <- rownames(cox_summary_df)
rownames(cox_summary_df) <- NULL
cox_summary_df$variable <- gsub("`", "", cox_summary_df$variable)

colnames(cox_summary_df)[colnames(cox_summary_df) == "coef"] <- "Estimate"
colnames(cox_summary_df)[colnames(cox_summary_df) == "se(coef)"] <- "Std. Error"

total_n <- final_cox_model$n
total_events <- final_cox_model$nevent

all_vars <- all.vars(final_cox_model$formula)

covariates_list <- paste(all_vars[3:length(all_vars)], collapse = ", ")

plot_data <- cox_summary_df %>%
  mutate(
    hazard_ratio = exp(Estimate),
    ci_low = exp(Estimate - 1.96 * `Std. Error`),
    ci_high = exp(Estimate + 1.96 * `Std. Error`),
    
    estimate_str = format(round(Estimate, 2), nsmall = 2),
    p_value_str = ifelse(`Pr(>|z|)` < 0.001, "< 0.001", format(round(`Pr(>|z|)`, 3), nsmall = 3)),
    hr_ci_str = paste0(format(round(hazard_ratio, 2), nsmall = 2), 
                       " (", format(round(ci_low, 2), nsmall = 2), "-", 
                       format(round(ci_high, 2), nsmall = 2), ")"),

    variable_wrapped = str_wrap(variable, width = 20)
  ) %>%
  arrange(hazard_ratio) %>%
  mutate(variable_wrapped = factor(variable_wrapped, levels = .$variable_wrapped))



p_table <- ggplot(plot_data, aes(y = variable_wrapped)) +

  geom_text(aes(x = 0, label = variable_wrapped), hjust = 0, size = 3.5) +

  geom_text(aes(x = 2.5, label = p_value_str), hjust = 0.5) +

  geom_text(aes(x = 4.5, label = hr_ci_str), hjust = 0.5) +

  annotate("text", x = 0, y = nrow(plot_data) + 1.2, label = "Variables", hjust = 0, fontface = "bold") +
  annotate("text", x = 2.5, y = nrow(plot_data) + 1.2, label = "P-value", hjust = 0.5, fontface = "bold") +
  annotate("text", x = 4.5, y = nrow(plot_data) + 1.2, label = "Hazard Ratio (95% CI)", hjust = 0.5, fontface = "bold") +
  
  theme_void() +
  coord_cartesian(xlim = c(0, 6), ylim = c(0.5, nrow(plot_data) + 1.5))

  annotate("text", x = 0, y = nrow(plot_data) + 1.2, label = "Variable", hjust = 0, fontface = "bold") +
  annotate("text", x = 1.5, y = nrow(plot_data) + 1.2, label = "Estimate", hjust = 0.5, fontface = "bold") +
  annotate("text", x = 2.5, y = nrow(plot_data) + 1.2, label = "Std. Error", hjust = 0.5, fontface = "bold") +
  annotate("text", x = 3.5, y = nrow(plot_data) + 1.2, label = "P-value", hjust = 0.5, fontface = "bold") +
  annotate("text", x = 4.7, y = nrow(plot_data) + 1.2, label = "Odds Ratio (95% CI)", hjust = 0.5, fontface = "bold") +
  
  theme_void() +
  coord_cartesian(xlim = c(0, 5.7))


p_forest <- ggplot(plot_data, aes(x = hazard_ratio, y = variable_wrapped)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey50") +
  geom_point(size = 2.5, color = "firebrick") + 
  geom_errorbarh(aes(xmin = ci_low, xmax = ci_high), height = 0.2, color = "firebrick") +
  labs(x = "Hazard Ratio",
  title = " ") +
  scale_x_log10() +
  theme_classic() +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.y = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major.x = element_line(color = "grey90", linetype = "dotted")
  )

final_plot <- p_table + p_forest + plot_layout(widths = c(2, 1))


pdf("figs/figS3/forest_plot_survival_response.pdf", width = 10, height = 6)
print(final_plot)
dev.off()
