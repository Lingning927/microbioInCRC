

get_group_label <- function(patient_info, problem) {
  n <- dim(patient_info)[2]
  group_label <- rep(" ", n)
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
    group_label[!(group_label %in% c("Left", "Rectum", "Right"))] <- "-"
  }else if (problem == "Location2") {
    location <- as.character(patient_info["Location", ])
    location[which(location == "left")] <- "Left"
    location[which(location == "right")] <- "Right"
    group_label <- location
    group_label[!(group_label %in% c("Right"))] <- "Left"
  }else if (problem == "survival(C is alive, W is dead)") {
    group_label <- as.character(patient_info[problem, ])
    group_label <- ifelse(group_label == "W", "Dead", "Alive")
  }else if (problem == "M") {
    group_label <- as.character(patient_info[problem, ])
    group_label[group_label == "MO"] <- "M0"
  }else if (problem == "N") {
    group_label <- as.character(patient_info[problem, ])
    group_label[(group_label %in% c("N1", "N2"))] <- "N1+N2"
  }else if (problem == "T") {
    group_label <- as.character(patient_info[problem, ])
    group_label[(group_label %in% c("T1", "T2"))] <- "T1+T2"
    group_label[(group_label %in% c("T3", "T4"))] <- "T3+T4"
  }else if (problem == "MSI") {
    group_label <- as.character(patient_info[problem, ])
    group_label[group_label == "MSI-H"] <- "MSI"
    group_label[group_label == "MSI-L"] <- "MSI"
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

deal_problem <- function(feature_summed, patient_info, problem) {
  if (problem == "N_and_T") {
    feature_summed <- feature_summed[which(substr(rownames(feature_summed), 1, 1) %in% c("N", "T")), ]
    y <- ifelse(substr(rownames(feature_summed), 1, 1) == "N", "Noraml", "CRC")
    response <- as.factor(y)
    return(list(feature = feature_summed, response = response))
  }else {
    feature_summed <- feature_summed[which(substr(rownames(feature_summed), 1, 1) == "T"), ]
    patient_info <- patient_info[, rownames(feature_summed)]
  }
  if (problem == "I.II_III.IV") {
    y <- get_group_label(patient_info, "I.II.III.IV")
    response <- as.factor(ifelse((y == "III" | y == "IV"), "III.IV", "I.II"))
  }else if (problem == "I.II.III_IV") {
    y <- get_group_label(patient_info, "I.II.III.IV")
    response <- as.factor(ifelse(y == "IV", "IV", "I_II_III"))
  }else if (problem == "I_II.III.IV") {
    y <- get_group_label(patient_info, "I.II.III.IV")
    response <- as.factor(ifelse(y == "I", "I", "II_III_IV"))
  }else if (problem == "HD_MD") {
    y <- get_group_label(patient_info, "differentiation")
    response <- as.factor(ifelse((y == "HD" | y == "MD"), "HD_MD", "Other"))
  }else if (problem == "Location") {
    y <- get_group_label(patient_info, "Location2")
    response <- as.factor(y)
  }else if (problem == "age") {
    y <- get_group_label(patient_info, "age")
    response <- as.factor(y)
  }else if (problem == "survival") {
    y <- get_group_label(patient_info, "survival(C is alive, W is dead)")
    response <- as.factor(y)
  }else{
    y <- get_group_label(patient_info, problem)
    response <- as.factor(y)
  }
  return(list(feature = feature_summed, response = response))
}

deal_compares_tr <- function(compares, feature_summed, response) {
    response <- as.character(response)
    if(compares[1] == "Normal_Polyp") {
        response[response == "Normal"] <- "Normal_Polyp"
        response[response == "Polyp"] <- "Normal_Polyp"
    }
    if (compares[2] == "CRC_1") {
        patient_info <- readRDS("adata2025/data/patient_info.rds")
        patient_names <- rownames(feature_summed[response == "CRC", ])
        stage <- patient_info["I.II.III.IV", patient_names]
        CRC_1_names <- colnames(stage)[stage[1, ] == "I"]
        first_names <- rownames(feature_summed[response == compares[1], ])
        feature_new <- feature_summed[c(first_names, CRC_1_names), ]
        response_new <- substr(rownames(feature_new), 1, 1)
        response_new = case_when(
            response_new == "N" ~ "Normal",
            response_new == "a" ~ "Polyp",
            response_new == "T" ~ "CRC_1"
        )
        if (compares[1] == "Normal_Polyp") {
            response_new[response_new == "Normal"] <- compares[1]
            response_new[response_new == "Polyp"] <- compares[1]
        }
        response_new <- factor(response_new, levels = compares)
    }else if(compares[2] == "Polyp_CRC_1") {
        patient_info <- readRDS("adata2025/data/patient_info.rds")
        patient_names <- rownames(feature_summed[response == "CRC", ])
        stage <- patient_info["I.II.III.IV", patient_names]
        CRC_1_names <- colnames(stage)[stage[1, ] == "I"]
        first_names <- rownames(feature_summed[response == compares[1], ])
        polyp_names <- rownames(feature_summed[response == "Polyp", ])
        feature_new <- feature_summed[c(first_names, CRC_1_names, polyp_names), ]
        response_new <- substr(rownames(feature_new), 1, 1)
        response_new = case_when(
            response_new == "N" ~ "Normal",
            response_new == "a" ~ "Polyp_CRC_1",
            response_new == "T" ~ "Polyp_CRC_1"
        )
        response_new <- factor(response_new, levels = compares)
    }else {
        feature_new <- feature_summed[response %in% compares, ]
        response_new <- substr(rownames(feature_new), 1, 1)
        response_new = case_when(
            response_new == "N" ~ "Normal",
            response_new == "a" ~ "Polyp",
            response_new == "T" ~ "CRC"
        )
        if (compares[1] == "Normal_Polyp") {
            response_new[response_new == "Normal"] <- compares[1]
            response_new[response_new == "Polyp"] <- compares[1]
        }
        response_new <- factor(response_new, levels = compares)
    }
    return(list(response = response_new, feature = feature_new))
}

test_roc <- function(train_feature, train_response, test_feature, test_response) {
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