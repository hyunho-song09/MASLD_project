# =============================================================================
# MASLD (Metabolic Dysfunction-Associated Steatotic Liver Disease) 
# Metabolomics Analysis Pipeline
# =============================================================================
# Description: Comprehensive analysis pipeline for metabolomics data comparing
#              NC (Normal Control) vs MASLD groups including statistical testing,
#              correlation analysis, machine learning models, and validation
# Author: [Your Name]
# Date: [Current Date]
# =============================================================================

# Environment Setup
# =============================================================================
# Clear workspace and set working directory
rm(list = ls())
setwd("")

# Load required libraries
required_packages <- c(
  "xlsx", "dplyr", "ggplot2", "tidyverse", "pheatmap", "Hmisc", 
  "svglite", "psych", "LMSstat", "caret", "pROC", "ROCR", "readxl"
)

# Install missing packages and load all
for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}

# =============================================================================
# STEP 1: DATA IMPORT AND PREPROCESSING
# =============================================================================

#' Load and preprocess metabolomics data
#' @description Imports raw data, applies z-score normalization to metabolite columns
load_and_preprocess_data <- function() {
  cat("Loading and preprocessing data...\n")
  
  # Import raw data and mapping information
  raw_df <- read.xlsx("data/250721_MASLD.xlsx", sheetName = "df")
  mapping_df <- read.xlsx("data/250721_MASLD.xlsx", sheetName = "mapping")
  
  # Create normalized dataset (z-score normalization for metabolite columns)
  df <- raw_df
  df[, 3:ncol(df)] <- scale(df[, 3:ncol(df)])
  
  cat("Data loaded successfully!\n")
  cat("Samples:", nrow(df), "\n")
  cat("Metabolites:", ncol(df) - 2, "\n")
  cat("Groups:", unique(df$Group), "\n")
  
  return(list(raw = raw_df, normalized = df, mapping = mapping_df))
}

# =============================================================================
# STEP 2: STATISTICAL ANALYSIS AND VOLCANO PLOT
# =============================================================================

#' Perform comprehensive statistical analysis
#' @param df Normalized dataframe
#' @param raw_df Raw dataframe for fold change calculation
perform_statistical_analysis <- function(df, raw_df) {
  cat("Performing statistical analysis...\n")
  
  # LMSstat analysis with FDR correction
  df_stat <- All_stats(df, Adjust_p_value = TRUE, Adjust_method = "BH", parallel = FALSE)
  
  # Create boxplot visualization
  Boxplot(df_stat,
          color = c("darkgreen", "darkred"),
          width = 0.5,
          order = c("NC", "MASLD"),
          fig_width = 3,
          fig_height = 4,
          sig_int = c(0.15, 0.1),
          asterisk = "t_test")
  
  # Manual t-test and volcano plot preparation
  results <- perform_manual_ttest(df, raw_df)
  create_volcano_plot(results)
  
  return(list(lms_results = df_stat, manual_results = results))
}

#' Perform manual t-test for volcano plot
#' @param df Normalized dataframe
#' @param raw_df Raw dataframe
perform_manual_ttest <- function(df, raw_df) {
  # Separate groups (using normalized data for statistical testing)
  nc_data <- df[df$Group == "NC", -c(1:2)]
  masld_data <- df[df$Group == "MASLD", -c(1:2)]
  
  # Separate raw data for fold change calculation
  nc_raw <- raw_df[raw_df$Group == "NC", -c(1:2)]
  masld_raw <- raw_df[raw_df$Group == "MASLD", -c(1:2)]
  
  # Initialize results dataframe
  results <- data.frame(
    metabolite = colnames(nc_data),
    log2fc = numeric(ncol(nc_data)),
    pvalue = numeric(ncol(nc_data)),
    stringsAsFactors = FALSE
  )
  
  # Perform t-test for each metabolite
  for (i in 1:ncol(nc_data)) {
    # Extract values for current metabolite (normalized data)
    nc_values <- as.numeric(nc_data[, i])
    masld_values <- as.numeric(masld_data[, i])
    
    # Extract raw values for fold change calculation
    nc_raw_values <- as.numeric(nc_raw[, i])
    masld_raw_values <- as.numeric(masld_raw[, i])
    
    # Remove NA values
    nc_values <- nc_values[!is.na(nc_values)]
    masld_values <- masld_values[!is.na(masld_values)]
    nc_raw_values <- nc_raw_values[!is.na(nc_raw_values)]
    masld_raw_values <- masld_raw_values[!is.na(masld_raw_values)]
    
    # Calculate log2 fold change using raw data
    mean_masld_raw <- mean(masld_raw_values, na.rm = TRUE)
    mean_nc_raw <- mean(nc_raw_values, na.rm = TRUE)
    
    if (mean_nc_raw > 0) {
      results$log2fc[i] <- log2(mean_masld_raw / mean_nc_raw)
    } else {
      results$log2fc[i] <- NA
    }
    
    # Perform t-test using normalized data
    if (length(nc_values) > 1 && length(masld_values) > 1) {
      test_result <- t.test(masld_values, nc_values)
      results$pvalue[i] <- test_result$p.value
    } else {
      results$pvalue[i] <- NA
    }
  }
  
  # Remove first row (seems to be a header issue in original code)
  results <- results[-1, ]
  
  # FDR correction
  results$fdr <- p.adjust(results$pvalue, method = "BH")
  
  # Add significance classification
  results$significant <- ifelse(
    results$fdr < 0.15,
    ifelse(results$log2fc > 0, "Significant Positive", "Significant Negative"),
    "Not Significant"
  )
  
  # Calculate -log10(p-value) for volcano plot
  results$neg_log10_pvalue <- -log10(results$pvalue)
  
  return(results[complete.cases(results), ])
}

#' Create volcano plot
#' @param results Results from t-test analysis
create_volcano_plot <- function(results) {
  volcano_plot <- ggplot(results, aes(x = log2fc, y = neg_log10_pvalue)) +
    geom_point(aes(color = significant), alpha = 0.7, size = 2) +
    scale_color_manual(
      values = c("Not Significant" = "gray70", 
                 "Significant Positive" = "red", 
                 "Significant Negative" = "blue")
    ) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", alpha = 0.7) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50", alpha = 0.7) +
    labs(
      x = "Log2 Fold Change (MASLD vs NC)",
      y = "-log10(p-value)",
      title = "Volcano Plot: MASLD vs NC Metabolite Analysis",
      subtitle = paste(
        "Significant metabolites (FDR < 0.15): Positive =", 
        sum(results$significant == "Significant Positive"),
        ", Negative =", 
        sum(results$significant == "Significant Negative")
      ),
      color = "Significance"
    ) +
    theme_classic() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 12),
      legend.position = "top"
    )
  
  # Save volcano plot
  ggsave("outputs/volcano_plot.svg", plot = volcano_plot, 
         device = "svg", width = 5.2, height = 5)
  
  cat("Volcano plot saved to outputs/volcano_plot.svg\n")
  return(volcano_plot)
}

# =============================================================================
# STEP 3: CORRELATION ANALYSIS
# =============================================================================

#' Perform comprehensive correlation analysis
#' @param df Normalized dataframe
perform_correlation_analysis <- function(df) {
  cat("Performing correlation analysis...\n")
  
  # Create output directory
  if (!dir.exists("outputs")) dir.create("outputs")
  
  # Separate groups
  data_nc <- df %>% filter(Group == "NC") %>% select(-Sample, -Group)
  data_masld <- df %>% filter(Group == "MASLD") %>% select(-Sample, -Group)
  
  # Correlation analysis for each group
  cor_results <- calculate_group_correlations(data_nc, data_masld)
  
  # Create heatmaps
  create_correlation_heatmaps(cor_results$cor_nc, cor_results$cor_masld)
  
  # Find significantly different correlations between groups
  sig_pairs <- find_significant_correlation_differences(
    cor_results$cor_nc, cor_results$cor_masld, data_nc, data_masld, df
  )
  
  # Create scatter plots for significant pairs
  create_significant_scatter_plots(df, sig_pairs)
  
  return(sig_pairs)
}

#' Calculate correlations for each group
#' @param data_nc NC group data
#' @param data_masld MASLD group data
calculate_group_correlations <- function(data_nc, data_masld) {
  # Calculate correlations with p-values
  cor_nc <- rcorr(as.matrix(data_nc), type = "pearson")
  cor_masld <- rcorr(as.matrix(data_masld), type = "pearson")
  
  return(list(cor_nc = cor_nc, cor_masld = cor_masld))
}

#' Create correlation heatmaps for both groups
#' @param cor_nc NC correlation results
#' @param cor_masld MASLD correlation results
create_correlation_heatmaps <- function(cor_nc, cor_masld) {
  # NC group heatmap with clustering
  pheatmap_nc <- pheatmap(
    cor_nc$r,
    main = "Correlation Heatmap (NC Group)",
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    show_rownames = FALSE,
    show_colnames = FALSE,
    filename = "outputs/correlation_heatmap_NC.svg",
    width = 10,
    height = 10
  )
  
  # MASLD group heatmap using NC clustering order
  row_order <- pheatmap_nc$tree_row$order
  col_order <- pheatmap_nc$tree_col$order
  
  pheatmap(
    cor_masld$r[row_order, col_order],
    main = "Correlation Heatmap (MASLD Group - Ordered by NC)",
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    show_rownames = FALSE,
    show_colnames = FALSE,
    filename = "outputs/correlation_heatmap_MASLD.svg",
    width = 10,
    height = 10
  )
  
  cat("Correlation heatmaps saved to outputs/\n")
}

#' Find significantly different correlations between groups
#' @param cor_nc NC correlation results
#' @param cor_masld MASLD correlation results
#' @param data_nc NC group data
#' @param data_masld MASLD group data
#' @param df Original dataframe for additional statistics
find_significant_correlation_differences <- function(cor_nc, cor_masld, data_nc, data_masld, df) {
  metabolites <- colnames(data_nc)
  n_met <- length(metabolites)
  comparison_results <- list()
  
  n_nc <- nrow(data_nc)
  n_masld <- nrow(data_masld)
  
  # Calculate CV function
  cv_percent <- function(x) {
    m <- mean(x, na.rm = TRUE)
    if (isTRUE(all.equal(m, 0))) return(NA_real_)
    100 * sd(x, na.rm = TRUE) / abs(m)
  }
  
  # Calculate CVs for each group
  raw_nc <- df %>% filter(Group == "NC") %>% select(-Sample, -Group)
  raw_masld <- df %>% filter(Group == "MASLD") %>% select(-Sample, -Group)
  
  cv_nc_vec <- sapply(raw_nc, cv_percent)
  cv_masld_vec <- sapply(raw_masld, cv_percent)
  
  # Calculate metabolite statistics
  met_stats <- df %>%
    select(Group, all_of(metabolites)) %>%
    pivot_longer(-Group, names_to = "met", values_to = "value") %>%
    group_by(met) %>%
    summarise(
      mean_nc = mean(value[Group == "NC"], na.rm = TRUE),
      mean_masld = mean(value[Group == "MASLD"], na.rm = TRUE),
      fold_change = ifelse(isTRUE(all.equal(mean_nc, 0)), NA_real_, mean_masld / mean_nc),
      p_value = tryCatch(t.test(value ~ Group)$p.value, error = function(e) NA_real_),
      .groups = 'drop'
    )
  
  # Compare correlations between groups using Fisher's Z-test
  for (i in 1:(n_met - 1)) {
    for (j in (i + 1):n_met) {
      met1 <- metabolites[i]
      met2 <- metabolites[j]
      
      # Extract correlation coefficients
      r_nc_val <- cor_nc$r[met1, met2]
      r_masld_val <- cor_masld$r[met1, met2]
      
      # Fisher's Z-test for correlation difference
      z_nc <- psych::fisherz(r_nc_val)
      z_masld <- psych::fisherz(r_masld_val)
      se_diff <- sqrt(1/(n_nc - 3) + 1/(n_masld - 3))
      z_stat <- (z_nc - z_masld) / se_diff
      p_diff_val <- 2 * pnorm(abs(z_stat), lower.tail = FALSE)
      
      # Store results
      comparison_results[[paste(met1, met2, sep = "_vs_")]] <- data.frame(
        met1 = met1,
        met2 = met2,
        r_nc = r_nc_val,
        r_masld = r_masld_val,
        p_nc = cor_nc$P[met1, met2],
        p_masld = cor_masld$P[met1, met2],
        p_diff = p_diff_val
      )
    }
  }
  
  # Convert to dataframe and filter significant pairs
  comparison_df <- do.call(rbind, comparison_results)
  
  significant_pairs <- comparison_df %>%
    filter(p_diff < 0.05) %>%
    arrange(p_diff) %>%
    mutate(
      met1 = as.character(met1),
      met2 = as.character(met2)
    ) %>%
    rowwise() %>%
    mutate(
      met1_cv_nc = cv_nc_vec[[met1]],
      met2_cv_nc = cv_nc_vec[[met2]],
      met1_cv_masld = cv_masld_vec[[met1]],
      met2_cv_masld = cv_masld_vec[[met2]]
    ) %>%
    ungroup() %>%
    left_join(met_stats %>% select(met, fold_change, p_value), by = c("met1" = "met")) %>%
    rename(met1_FC = fold_change, met1_pval = p_value) %>%
    left_join(met_stats %>% select(met, fold_change, p_value), by = c("met2" = "met")) %>%
    rename(met2_FC = fold_change, met2_pval = p_value) %>%
    mutate(p_adj = p.adjust(p_diff, method = "BH"))
  
  # Save results
  write.xlsx(significant_pairs, "outputs/significant_correlation_pairs.xlsx")
  
  cat("Found", nrow(significant_pairs), "significantly different correlation pairs\n")
  return(significant_pairs)
}

#' Create scatter plots for significant correlation pairs
#' @param df Original dataframe
#' @param significant_pairs Significant correlation pairs
create_significant_scatter_plots <- function(df, significant_pairs) {
  if (!dir.exists("outputs/significant_plots")) {
    dir.create("outputs/significant_plots", recursive = TRUE)
  }
  
  cat("Creating scatter plots for significant pairs...\n")
  
  for (i in 1:min(10, nrow(significant_pairs))) {  # Limit to top 10 pairs
    pair_info <- significant_pairs[i, ]
    met1_name <- as.character(pair_info$met1)
    met2_name <- as.character(pair_info$met2)
    
    # Prepare plot data
    plot_data <- df %>%
      select(Group, all_of(c(met1_name, met2_name)))
    
    # Create annotation text
    annotation_text <- paste(
      sprintf("NC: r = %.3f, p = %.3f", pair_info$r_nc, pair_info$p_nc),
      sprintf("MASLD: r = %.3f, p = %.3f", pair_info$r_masld, pair_info$p_masld),
      sprintf("Difference p-value: %.4f", pair_info$p_diff),
      sep = "\n"
    )
    
    # Create scatter plot
    p <- ggplot(plot_data, aes_string(x = met1_name, y = met2_name, color = "Group")) +
      geom_point(alpha = 0.7, size = 2.5) +
      geom_smooth(method = "lm", se = TRUE, aes(group = Group), alpha = 0.3) +
      scale_color_manual(values = c("NC" = "darkgreen", "MASLD" = "darkred")) +
      labs(
        title = paste("Correlation Analysis:", met1_name, "vs", met2_name),
        x = met1_name,
        y = met2_name,
        caption = annotation_text
      ) +
      theme_classic() +
      theme(
        plot.caption = element_text(hjust = 0, size = 9),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.position = "top"
      )
    
    # Save plot
    file_name <- paste0("outputs/significant_plots/", met1_name, "_vs_", met2_name, ".svg")
    ggsave(file = file_name, plot = p, device = "svg", width = 7, height = 5)
  }
  
  cat("Scatter plots saved to outputs/significant_plots/\n")
}

# =============================================================================
# STEP 4: MACHINE LEARNING AND ROC ANALYSIS
# =============================================================================

#' Perform comprehensive machine learning analysis with overfitting detection
#' @param df Normalized dataframe
perform_ml_analysis <- function(df) {
  cat("Performing machine learning analysis...\n")
  
  # Define feature sets based on previous analysis
  mean_diff_features <- c("M_107", "M_109", "M_44", "M_185", "M_211", "M_111")
  interaction_features <- c("M_000", "M_999")  # These should be replaced with actual significant metabolites
  
  # Prepare datasets
  ml_data <- prepare_ml_datasets(df, mean_diff_features, interaction_features)
  
  # Train models
  models <- train_ml_models(ml_data)
  
  # Evaluate models and create ROC plots
  evaluation_results <- evaluate_models(models, ml_data)
  
  return(evaluation_results)
}

#' Prepare datasets for machine learning
#' @param df Normalized dataframe
#' @param mean_diff_features Features based on mean differences
#' @param interaction_features Features for interaction analysis
prepare_ml_datasets <- function(df, mean_diff_features, interaction_features) {
  # Create datasets
  m_df <- df[, c("Sample", "Group", mean_diff_features)]
  i_df <- df[, c("Sample", "Group", interaction_features)]
  
  # Convert Group to factor
  m_df$Group <- as.factor(m_df$Group)
  i_df$Group <- as.factor(i_df$Group)
  
  # Train/test split (70:30)
  set.seed(42)
  train_idx <- createDataPartition(i_df$Group, p = 0.7, list = FALSE)
  
  train_i_df <- i_df[train_idx, ]
  test_i_df <- i_df[-train_idx, ]
  train_m_df <- m_df[train_idx, ]
  test_m_df <- m_df[-train_idx, ]
  
  # Create derived features
  ml_data <- create_derived_features(train_i_df, test_i_df, train_m_df, test_m_df)
  
  return(ml_data)
}

#' Create derived features for enhanced model performance
#' @param train_i_df Training interaction dataset
#' @param test_i_df Test interaction dataset
#' @param train_m_df Training mean difference dataset
#' @param test_m_df Test mean difference dataset
create_derived_features <- function(train_i_df, test_i_df, train_m_df, test_m_df) {
  k <- 5  # Number of neighbors for local correlation
  
  # For interaction dataset (M_000, M_999)
  if (ncol(train_i_df) >= 4) {  # Ensure we have the expected columns
    # Calculate thresholds for local correlation
    thr_i_000 <- quantile(abs(diff(train_i_df[, 3])), 0.25, na.rm = TRUE)
    thr_i_999 <- quantile(abs(diff(train_i_df[, 4])), 0.25, na.rm = TRUE)
    
    # Local correlation for training set
    train_i_df$local_correlation <- calculate_local_correlation(
      train_i_df[, 3], train_i_df[, 4], thr_i_000, thr_i_999, k
    )
    
    # Group-centered interaction for training set
    train_i_df$group_mean_interaction <- calculate_group_interaction(
      train_i_df, colnames(train_i_df)[3], colnames(train_i_df)[4]
    )
    
    # Apply same transformations to test set
    test_i_df$local_correlation <- calculate_local_correlation_test(
      test_i_df[, 3], test_i_df[, 4], train_i_df[, 3], train_i_df[, 4], 
      thr_i_000, thr_i_999, k
    )
    
    grp_mu_i <- aggregate(cbind(train_i_df[, 3], train_i_df[, 4]) ~ Group, 
                          data = train_i_df, mean)
    colnames(grp_mu_i)[2:3] <- colnames(train_i_df)[3:4]
    
    test_i_df$group_mean_interaction <- mapply(
      function(g, x, y) {
        mu_x <- grp_mu_i[grp_mu_i$Group == g, 2]
        mu_y <- grp_mu_i[grp_mu_i$Group == g, 3]
        (x - mu_x) * (y - mu_y)
      },
      test_i_df$Group, test_i_df[, 3], test_i_df[, 4]
    )
  }
  
  # Similar process for mean difference dataset
  if (ncol(train_m_df) >= 4) {
    thr_m_1 <- quantile(abs(diff(train_m_df[, 3])), 0.25, na.rm = TRUE)
    thr_m_2 <- quantile(abs(diff(train_m_df[, 4])), 0.25, na.rm = TRUE)
    
    train_m_df$local_correlation <- calculate_local_correlation(
      train_m_df[, 3], train_m_df[, 4], thr_m_1, thr_m_2, k
    )
    
    train_m_df$group_mean_interaction <- calculate_group_interaction(
      train_m_df, colnames(train_m_df)[3], colnames(train_m_df)[4]
    )
    
    test_m_df$local_correlation <- calculate_local_correlation_test(
      test_m_df[, 3], test_m_df[, 4], train_m_df[, 3], train_m_df[, 4],
      thr_m_1, thr_m_2, k
    )
    
    grp_mu_m <- aggregate(cbind(train_m_df[, 3], train_m_df[, 4]) ~ Group,
                          data = train_m_df, mean)
    colnames(grp_mu_m)[2:3] <- colnames(train_m_df)[3:4]
    
    test_m_df$group_mean_interaction <- mapply(
      function(g, x, y) {
        mu_x <- grp_mu_m[grp_mu_m$Group == g, 2]
        mu_y <- grp_mu_m[grp_mu_m$Group == g, 3]
        (x - mu_x) * (y - mu_y)
      },
      test_m_df$Group, test_m_df[, 3], test_m_df[, 4]
    )
  }
  
  return(list(
    train_i = train_i_df, test_i = test_i_df,
    train_m = train_m_df, test_m = test_m_df
  ))
}

#' Calculate local correlation for a sample
#' @param x First metabolite values
#' @param y Second metabolite values
#' @param thr_x Threshold for first metabolite
#' @param thr_y Threshold for second metabolite
#' @param k Number of neighbors
calculate_local_correlation <- function(x, y, thr_x, thr_y, k) {
  sapply(seq_along(x), function(row) {
    nb <- which(abs(x - x[row]) < thr_x & abs(y - y[row]) < thr_y)
    nb <- nb[nb != row]
    if (length(nb) > k) nb <- nb[1:k]
    if (length(nb) > 1) cor(x[nb], y[nb], use = "complete.obs") else 0
  })
}

#' Calculate local correlation for test set using training set neighbors
calculate_local_correlation_test <- function(test_x, test_y, train_x, train_y, thr_x, thr_y, k) {
  sapply(seq_along(test_x), function(row) {
    nb <- which(abs(train_x - test_x[row]) < thr_x & abs(train_y - test_y[row]) < thr_y)
    if (length(nb) > k) nb <- nb[1:k]
    if (length(nb) > 1) cor(train_x[nb], train_y[nb], use = "complete.obs") else 0
  })
}

#' Calculate group-centered interaction terms
#' @param data Dataset
#' @param var1 First variable name
#' @param var2 Second variable name
calculate_group_interaction <- function(data, var1, var2) {
  grp_means <- aggregate(data[, c(var1, var2)], by = list(data$Group), mean, na.rm = TRUE)
  colnames(grp_means) <- c("Group", var1, var2)
  
  mapply(
    function(g, x, y) {
      mu_x <- grp_means[grp_means$Group == g, var1]
      mu_y <- grp_means[grp_means$Group == g, var2]
      (x - mu_x) * (y - mu_y)
    },
    data$Group, data[[var1]], data[[var2]]
  )
}

#' Train multiple machine learning models
#' @param ml_data Prepared datasets
train_ml_models <- function(ml_data) {
  # Cross-validation setup
  ctrl <- trainControl(
    method = "cv",
    number = 5,
    repeats = 3,
    summaryFunction = twoClassSummary,
    classProbs = TRUE,
    savePredictions = TRUE
  )
  
  cat("Training machine learning models...\n")
  
  # Train models with different approaches
  models <- list()
  
  set.seed(42)
  models$mean_diff <- train(
    Group ~ .,
    data = ml_data$train_m[, -1],  # Remove Sample column
    method = "glm",
    family = "binomial",
    trControl = ctrl,
    metric = "ROC"
  )
  
  set.seed(42)
  models$interaction_glm <- train(
    Group ~ . + I(local_correlation^2),  # Add polynomial term
    data = ml_data$train_i[, c("Group", "local_correlation", "group_mean_interaction")],
    method = "glm",
    family = "binomial",
    trControl = ctrl,
    metric = "ROC"
  )
  
  set.seed(42)
  models$rf_mean <- train(
    Group ~ .,
    data = ml_data$train_m[, -1],
    method = "rf",
    trControl = ctrl,
    metric = "ROC",
    ntree = 500,
    tuneLength = 3
  )
  
  set.seed(42)
  models$rf_interaction <- train(
    Group ~ .,
    data = ml_data$train_i[, c("Group", "local_correlation", "group_mean_interaction")],
    method = "rf",
    trControl = ctrl,
    metric = "ROC",
    ntree = 500,
    tuneLength = 3
  )
  
  cat("Models trained successfully!\n")
  return(models)
}

#' Evaluate models and create ROC analysis
#' @param models Trained models
#' @param ml_data ML datasets
evaluate_models <- function(models, ml_data) {
  cat("Evaluating model performance...\n")
  
  # Make predictions on test sets
  predictions <- make_predictions(models, ml_data)
  
  # Calculate ROC curves
  roc_curves <- calculate_roc_curves(predictions, ml_data)
  
  # Create performance comparison
  performance_table <- create_performance_table(roc_curves, predictions, ml_data)
  
  # Create ROC plot
  create_roc_plot(roc_curves)
  
  # Check for overfitting
  overfitting_analysis <- check_overfitting(models, ml_data)
  
  cat("Model evaluation completed!\n")
  
  return(list(
    predictions = predictions,
    roc_curves = roc_curves,
    performance = performance_table,
    overfitting = overfitting_analysis
  ))
}

#' Make predictions on test sets
#' @param models Trained models
#' @param ml_data ML datasets
make_predictions <- function(models, ml_data) {
  predictions <- list()
  
  # Mean difference model predictions
  predictions$mean_diff_train <- predict(models$mean_diff, ml_data$train_m, type = "prob")[, 2]
  predictions$mean_diff_test <- predict(models$mean_diff, ml_data$test_m, type = "prob")[, 2]
  
  # Interaction model predictions
  pred_data_i_train <- ml_data$train_i[, c("Group", "local_correlation", "group_mean_interaction")]
  pred_data_i_test <- ml_data$test_i[, c("Group", "local_correlation", "group_mean_interaction")]
  
  predictions$interaction_train <- predict(models$interaction_glm, pred_data_i_train, type = "prob")[, 2]
  predictions$interaction_test <- predict(models$interaction_glm, pred_data_i_test, type = "prob")[, 2]
  
  # Random Forest predictions
  predictions$rf_mean_train <- predict(models$rf_mean, ml_data$train_m, type = "prob")[, 2]
  predictions$rf_mean_test <- predict(models$rf_mean, ml_data$test_m, type = "prob")[, 2]
  
  predictions$rf_interaction_train <- predict(models$rf_interaction, pred_data_i_train, type = "prob")[, 2]
  predictions$rf_interaction_test <- predict(models$rf_interaction, pred_data_i_test, type = "prob")[, 2]
  
  return(predictions)
}

#' Calculate ROC curves for all models
#' @param predictions Model predictions
#' @param ml_data ML datasets
calculate_roc_curves <- function(predictions, ml_data) {
  roc_curves <- list()
  
  # Training ROC curves
  roc_curves$mean_diff_train <- roc(ml_data$train_m$Group, predictions$mean_diff_train)
  roc_curves$interaction_train <- roc(ml_data$train_i$Group, predictions$interaction_train)
  roc_curves$rf_mean_train <- roc(ml_data$train_m$Group, predictions$rf_mean_train)
  roc_curves$rf_interaction_train <- roc(ml_data$train_i$Group, predictions$rf_interaction_train)
  
  # Test ROC curves
  roc_curves$mean_diff_test <- roc(ml_data$test_m$Group, predictions$mean_diff_test)
  roc_curves$interaction_test <- roc(ml_data$test_i$Group, predictions$interaction_test)
  roc_curves$rf_mean_test <- roc(ml_data$test_m$Group, predictions$rf_mean_test)
  roc_curves$rf_interaction_test <- roc(ml_data$test_i$Group, predictions$rf_interaction_test)
  
  return(roc_curves)
}

#' Create comprehensive performance table
#' @param roc_curves ROC curve results
#' @param predictions Model predictions
#' @param ml_data ML datasets
create_performance_table <- function(roc_curves, predictions, ml_data) {
  # Helper function to calculate metrics
  calc_metrics <- function(actual, predicted_prob, threshold = 0.5) {
    predicted_class <- ifelse(predicted_prob > threshold, 
                              levels(actual)[2], levels(actual)[1])
    predicted_class <- factor(predicted_class, levels = levels(actual))
    
    cm <- confusionMatrix(predicted_class, actual)
    
    return(list(
      Sensitivity = cm$byClass["Sensitivity"],
      Specificity = cm$byClass["Specificity"],
      F1 = cm$byClass["F1"]
    ))
  }
  
  # Calculate metrics for all models
  model_names <- c("Mean_Diff_GLM", "Interaction_GLM", "Mean_Diff_RF", "Interaction_RF")
  
  performance_table <- data.frame(
    Model = rep(model_names, each = 2),
    Dataset = rep(c("Train", "Test"), times = 4),
    AUC = c(
      auc(roc_curves$mean_diff_train), auc(roc_curves$mean_diff_test),
      auc(roc_curves$interaction_train), auc(roc_curves$interaction_test),
      auc(roc_curves$rf_mean_train), auc(roc_curves$rf_mean_test),
      auc(roc_curves$rf_interaction_train), auc(roc_curves$rf_interaction_test)
    ),
    stringsAsFactors = FALSE
  )
  
  # Add detailed metrics
  train_metrics <- list(
    calc_metrics(ml_data$train_m$Group, predictions$mean_diff_train),
    calc_metrics(ml_data$train_i$Group, predictions$interaction_train),
    calc_metrics(ml_data$train_m$Group, predictions$rf_mean_train),
    calc_metrics(ml_data$train_i$Group, predictions$rf_interaction_train)
  )
  
  test_metrics <- list(
    calc_metrics(ml_data$test_m$Group, predictions$mean_diff_test),
    calc_metrics(ml_data$test_i$Group, predictions$interaction_test),
    calc_metrics(ml_data$test_m$Group, predictions$rf_mean_test),
    calc_metrics(ml_data$test_i$Group, predictions$rf_interaction_test)
  )
  
  all_metrics <- c(train_metrics, test_metrics)
  
  performance_table$Sensitivity <- sapply(all_metrics, function(x) x$Sensitivity)
  performance_table$Specificity <- sapply(all_metrics, function(x) x$Specificity)
  performance_table$F1_Score <- sapply(all_metrics, function(x) x$F1)
  
  # Save performance table
  write.xlsx(performance_table, "outputs/model_performance_comparison.xlsx")
  
  cat("\n=== Model Performance Comparison ===\n")
  print(performance_table)
  
  return(performance_table)
}

#' Create ROC curve visualization
#' @param roc_curves ROC curve results
create_roc_plot <- function(roc_curves) {
  # Base R ROC plot
  svg("outputs/ROC_curves_test.svg", width = 6, height = 6)
  
  plot(roc_curves$mean_diff_test, col = "blue", main = "ROC Curves (Test Data)", 
       lwd = 2, cex.main = 1.2)
  plot(roc_curves$interaction_test, col = "green", add = TRUE, lwd = 2)
  plot(roc_curves$rf_mean_test, col = "red", add = TRUE, lwd = 2)
  plot(roc_curves$rf_interaction_test, col = "purple", add = TRUE, lwd = 2)
  
  # Add legend with AUC values
  legend("bottomright", 
         legend = c(
           paste("Mean_Diff_GLM (AUC =", round(auc(roc_curves$mean_diff_test), 3), ")"),
           paste("Interaction_GLM (AUC =", round(auc(roc_curves$interaction_test), 3), ")"),
           paste("Mean_Diff_RF (AUC =", round(auc(roc_curves$rf_mean_test), 3), ")"),
           paste("Interaction_RF (AUC =", round(auc(roc_curves$rf_interaction_test), 3), ")")
         ),
         col = c("blue", "green", "red", "purple"), 
         lwd = 2,
         cex = 0.8)
  
  dev.off()
  
  # ggplot2 version for better aesthetics
  create_ggplot_roc(roc_curves)
  
  cat("ROC curves saved to outputs/\n")
}

#' Create ggplot2 ROC curve visualization
#' @param roc_curves ROC curve results
create_ggplot_roc <- function(roc_curves) {
  # Prepare data for ggplot
  roc_data <- data.frame(
    FPR = c(1 - roc_curves$mean_diff_test$specificities,
            1 - roc_curves$interaction_test$specificities,
            1 - roc_curves$rf_mean_test$specificities,
            1 - roc_curves$rf_interaction_test$specificities),
    TPR = c(roc_curves$mean_diff_test$sensitivities,
            roc_curves$interaction_test$sensitivities,
            roc_curves$rf_mean_test$sensitivities,
            roc_curves$rf_interaction_test$sensitivities),
    Model = c(rep("Mean_Diff_GLM", length(roc_curves$mean_diff_test$specificities)),
              rep("Interaction_GLM", length(roc_curves$interaction_test$specificities)),
              rep("Mean_Diff_RF", length(roc_curves$rf_mean_test$specificities)),
              rep("Interaction_RF", length(roc_curves$rf_interaction_test$specificities)))
  )
  
  # Create ggplot ROC curve
  p <- ggplot(roc_data, aes(x = FPR, y = TPR, color = Model)) +
    geom_line(size = 1.2) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray50") +
    scale_color_manual(
      values = c("Mean_Diff_GLM" = "blue", "Interaction_GLM" = "green",
                 "Mean_Diff_RF" = "red", "Interaction_RF" = "purple"),
      labels = c(
        paste("Mean_Diff_GLM (AUC =", round(auc(roc_curves$mean_diff_test), 3), ")"),
        paste("Interaction_GLM (AUC =", round(auc(roc_curves$interaction_test), 3), ")"),
        paste("Mean_Diff_RF (AUC =", round(auc(roc_curves$rf_mean_test), 3), ")"),
        paste("Interaction_RF (AUC =", round(auc(roc_curves$rf_interaction_test), 3), ")")
      )
    ) +
    labs(
      x = "False Positive Rate (1 - Specificity)",
      y = "True Positive Rate (Sensitivity)",
      title = "ROC Curves for Test Data",
      subtitle = "Comparison of Different Model Approaches"
    ) +
    theme_classic() +
    theme(
      legend.position = "bottom",
      legend.title = element_blank(),
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 12, hjust = 0.5)
    ) +
    coord_equal()
  
  ggsave("outputs/ROC_curves_ggplot.svg", plot = p, device = "svg", width = 8, height = 6)
}

#' Check for overfitting using various methods
#' @param models Trained models
#' @param ml_data ML datasets
check_overfitting <- function(models, ml_data) {
  cat("Checking for overfitting...\n")
  
  overfitting_results <- list()
  
  # Calculate AUC differences between train and test
  train_aucs <- c(
    auc(roc(ml_data$train_m$Group, predict(models$mean_diff, ml_data$train_m, type = "prob")[, 2])),
    auc(roc(ml_data$train_m$Group, predict(models$rf_mean, ml_data$train_m, type = "prob")[, 2]))
  )
  
  test_aucs <- c(
    auc(roc(ml_data$test_m$Group, predict(models$mean_diff, ml_data$test_m, type = "prob")[, 2])),
    auc(roc(ml_data$test_m$Group, predict(models$rf_mean, ml_data$test_m, type = "prob")[, 2]))
  )
  
  auc_differences <- train_aucs - test_aucs
  
  overfitting_results$auc_diff <- data.frame(
    Model = c("Mean_Diff_GLM", "Mean_Diff_RF"),
    Train_AUC = train_aucs,
    Test_AUC = test_aucs,
    Difference = auc_differences,
    Overfitting = ifelse(auc_differences > 0.1, "Likely", "Minimal")
  )
  
  cat("\n=== Overfitting Analysis ===\n")
  print(overfitting_results$auc_diff)
  
  return(overfitting_results)
}

# =============================================================================
# STEP 5: VALIDATION ANALYSIS
# =============================================================================

#' Perform validation analysis using external datasets
perform_validation_analysis <- function() {
  cat("Performing validation analysis...\n")
  
  # Load validation datasets
  val_01 <- read_excel("data/250726_val_01_plasma.xlsx", sheet = "test_all")[, -c(1:2)]
  val_03 <- read_excel("data/250726_val_03_plasma.xlsx", sheet = "test_all")[, -c(1:2)]
  
  # Clean validation dataset 03
  val_03 <- val_03[, -3]
  val_03 <- val_03[, c(2, 1)]
  
  # Standardize column names
  colnames(val_01) <- c("M_263", "M_180")
  colnames(val_03) <- c("M_263", "M_180")
  
  # Perform correlation analysis for each validation set
  validation_results <- list()
  
  validation_results$val_01 <- analyze_validation_dataset(val_01, "Validation Set 01")
  validation_results$val_03 <- analyze_validation_dataset(val_03, "Validation Set 03")
  
  # Create comparison plot
  create_validation_comparison_plot(validation_results)
  
  return(validation_results)
}

#' Analyze individual validation dataset
#' @param val_data Validation dataset
#' @param dataset_name Name of the dataset
analyze_validation_dataset <- function(val_data, dataset_name) {
  cat(paste("Analyzing", dataset_name, "...\n"))
  
  # Correlation analysis
  cor_test <- cor.test(val_data$M_180, val_data$M_263, method = "pearson")
  cor_matrix <- rcorr(as.matrix(val_data), type = "pearson")
  
  # Create scatter plot
  p <- ggplot(val_data, aes(x = M_180, y = M_263)) +
    geom_point(size = 3, alpha = 0.7, color = "steelblue") +
    geom_smooth(method = "lm", se = TRUE, color = "red", size = 1.2, alpha = 0.3) +
    labs(
      x = "M_180",
      y = "M_263", 
      title = paste("Scatter Plot:", dataset_name),
      subtitle = "M_263 vs M_180 with Linear Trend"
    ) +
    theme_classic() +
    theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 12, hjust = 0.5),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10)
    ) +
    annotate("text",
             x = min(val_data$M_180) + 0.05 * (max(val_data$M_180) - min(val_data$M_180)),
             y = max(val_data$M_263) - 0.05 * (max(val_data$M_263) - min(val_data$M_263)),
             label = paste0("Pearson r = ", round(cor_test$estimate, 4), "\n",
                            "p-value = ", format(cor_test$p.value, scientific = TRUE, digits = 3)),
             hjust = 0, vjust = 1, size = 4, fontface = "bold",
             color = "darkblue")
  
  # Save plot
  filename <- paste0("outputs/", gsub(" ", "_", tolower(dataset_name)), "_scatter.svg")
  ggsave(filename, plot = p, device = "svg", width = 6, height = 6)
  
  return(list(
    correlation = cor_test$estimate,
    p_value = cor_test$p.value,
    cor_matrix = cor_matrix,
    plot = p
  ))
}

#' Create comparison plot for validation results
#' @param validation_results Results from validation analysis
create_validation_comparison_plot <- function(validation_results) {
  # Combine validation data for comparison
  comparison_data <- data.frame(
    Dataset = c("Validation_01", "Validation_03"),
    Correlation = c(validation_results$val_01$correlation, 
                    validation_results$val_03$correlation),
    P_Value = c(validation_results$val_01$p_value,
                validation_results$val_03$p_value)
  )
  
  # Create comparison bar plot
  p <- ggplot(comparison_data, aes(x = Dataset, y = Correlation, fill = Dataset)) +
    geom_col(alpha = 0.7, width = 0.6) +
    geom_text(aes(label = paste("r =", round(Correlation, 3), "\np =", 
                                format(P_Value, scientific = TRUE, digits = 2))),
              vjust = -0.5, size = 4, fontface = "bold") +
    scale_fill_manual(values = c("steelblue", "darkorange")) +
    labs(
      x = "Validation Dataset",
      y = "Pearson Correlation Coefficient",
      title = "Validation Results Comparison",
      subtitle = "M_180 vs M_263 Correlation in External Datasets"
    ) +
    theme_classic() +
    theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 12, hjust = 0.5),
      legend.position = "none"
    ) +
    ylim(0, max(comparison_data$Correlation) * 1.2)
  
  ggsave("outputs/validation_comparison.svg", plot = p, device = "svg", width = 7, height = 5)
  
  cat("Validation comparison plot saved to outputs/validation_comparison.svg\n")
}

# =============================================================================
# MAIN EXECUTION PIPELINE
# =============================================================================

#' Main function to execute the entire analysis pipeline
#' @description Orchestrates the complete MASLD metabolomics analysis
main_analysis_pipeline <- function() {
  cat("=============================================================================\n")
  cat("STARTING MASLD METABOLOMICS ANALYSIS PIPELINE\n")
  cat("=============================================================================\n")
  
  # Create output directories
  if (!dir.exists("outputs")) dir.create("outputs")
  if (!dir.exists("outputs/significant_plots")) dir.create("outputs/significant_plots", recursive = TRUE)
  
  # Step 1: Data Loading and Preprocessing
  cat("\n--- STEP 1: DATA LOADING AND PREPROCESSING ---\n")
  data_list <- load_and_preprocess_data()
  
  # Step 2: Statistical Analysis
  cat("\n--- STEP 2: STATISTICAL ANALYSIS ---\n")
  stat_results <- perform_statistical_analysis(data_list$normalized, data_list$raw)
  
  # Step 3: Correlation Analysis
  cat("\n--- STEP 3: CORRELATION ANALYSIS ---\n")
  correlation_results <- perform_correlation_analysis(data_list$normalized)
  
  # Step 4: Machine Learning Analysis
  cat("\n--- STEP 4: MACHINE LEARNING ANALYSIS ---\n")
  ml_results <- perform_ml_analysis(data_list$normalized)
  
  # Step 5: Validation Analysis
  cat("\n--- STEP 5: VALIDATION ANALYSIS ---\n")
  validation_results <- perform_validation_analysis()
  
  # Summary Report
  cat("\n=============================================================================\n")
  cat("ANALYSIS PIPELINE COMPLETED SUCCESSFULLY!\n")
  cat("=============================================================================\n")
  cat("Results saved in 'outputs/' directory:\n")
  cat("- volcano_plot.svg: Differential metabolite analysis\n")
  cat("- correlation_heatmap_*.svg: Correlation heatmaps\n")
  cat("- ROC_curves_*.svg: Model performance visualization\n")
  cat("- significant_plots/: Scatter plots for significant correlations\n")
  cat("- *.xlsx: Detailed numerical results\n")
  cat("=============================================================================\n")
  
  return(list(
    data = data_list,
    statistics = stat_results,
    correlations = correlation_results,
    machine_learning = ml_results,
    validation = validation_results
  ))
}

# =============================================================================
# UTILITY FUNCTIONS
# =============================================================================

#' Create learning curve analysis
#' @param model_func Model training function
#' @param data Training data
#' @param sizes Vector of sample sizes to test
create_learning_curve <- function(model_func, data, sizes = seq(10, nrow(data), by = 10)) {
  aucs <- numeric(length(sizes))
  
  for (i in seq_along(sizes)) {
    if (sizes[i] <= nrow(data)) {
      sub_data <- data[sample(nrow(data), sizes[i]), ]
      
      tryCatch({
        set.seed(42)
        model <- model_func(sub_data)
        pred_prob <- predict(model, sub_data, type = "prob")[, 2]
        roc_obj <- roc(sub_data$Group, pred_prob)
        aucs[i] <- auc(roc_obj)
      }, error = function(e) {
        aucs[i] <<- NA
      })
    } else {
      aucs[i] <- NA
    }
  }
  
  return(data.frame(sample_size = sizes, auc = aucs))
}

#' Print session information for reproducibility
print_session_info <- function() {
  cat("\n=== SESSION INFORMATION ===\n")
  print(sessionInfo())
  
  # Save session info to file
  capture.output(sessionInfo(), file = "outputs/session_info.txt")
}

# =============================================================================
# EXECUTE ANALYSIS
# =============================================================================

# Run the complete analysis pipeline
if (interactive()) {
  results <- main_analysis_pipeline()
  print_session_info()
} else {
  cat("Script loaded. Run main_analysis_pipeline() to execute the complete analysis.\n")
}

# =============================================================================
# END OF SCRIPT
# =============================================================================