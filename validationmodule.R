# =============================================================================
# IMPUTATION QUALITY ASSESSMENT WITHOUT ORIGINAL DATA
# =============================================================================
# 
# Methods to evaluate imputation quality when original complete data is unavailable
#
# =============================================================================

library(mice)
library(VIM)
library(dplyr)
library(ggplot2)
library(corrplot)

# =============================================================================
# MAIN VALIDATION FUNCTION
# =============================================================================

validate_imputation_no_original <- function(original_with_missing, mice_object, completed_data) {
  
  cat("=============================================================================\n")
  cat("IMPUTATION QUALITY ASSESSMENT (NO ORIGINAL DATA REQUIRED)\n")
  cat("=============================================================================\n\n")
  
  results <- list()
  
  # 1. CONVERGENCE DIAGNOSTICS
  cat("1. CONVERGENCE DIAGNOSTICS\n")
  cat("-------------------------\n")
  
  convergence_results <- assess_convergence(mice_object)
  results$convergence <- convergence_results
  
  # 2. DISTRIBUTION PRESERVATION
  cat("\n2. DISTRIBUTION PRESERVATION\n")
  cat("----------------------------\n")
  
  distribution_results <- assess_distribution_preservation(original_with_missing, completed_data)
  results$distribution <- distribution_results
  
  # 3. CORRELATION STRUCTURE PRESERVATION
  cat("\n3. CORRELATION STRUCTURE PRESERVATION\n")
  cat("-------------------------------------\n")
  
  correlation_results <- assess_correlation_preservation(original_with_missing, completed_data)
  results$correlation <- correlation_results
  
  # 4. PLAUSIBILITY CHECKS
  cat("\n4. PLAUSIBILITY CHECKS\n")
  cat("----------------------\n")
  
  plausibility_results <- assess_plausibility(original_with_missing, completed_data)
  results$plausibility <- plausibility_results
  
  # 5. BETWEEN-IMPUTATION VARIABILITY
  cat("\n5. BETWEEN-IMPUTATION VARIABILITY\n")
  cat("---------------------------------\n")
  
  variability_results <- assess_between_imputation_variability(mice_object)
  results$variability <- variability_results
  
  # 6. CROSS-VALIDATION ASSESSMENT
  cat("\n6. CROSS-VALIDATION ASSESSMENT\n")
  cat("------------------------------\n")
  
  cv_results <- cross_validation_assessment(original_with_missing)
  results$cross_validation <- cv_results
  
  # 7. OVERALL QUALITY SCORE
  cat("\n7. OVERALL QUALITY ASSESSMENT\n")
  cat("-----------------------------\n")
  
  overall_score <- calculate_overall_quality_score(results)
  results$overall_score <- overall_score
  
  # 8. SAVE DETAILED REPORT
  save_validation_report(results, "imputation_quality_report.txt")
  
  cat("\n=============================================================================\n")
  cat("VALIDATION COMPLETED - See 'imputation_quality_report.txt' for full details\n")
  cat("=============================================================================\n")
  
  return(results)
}

# =============================================================================
# 1. CONVERGENCE DIAGNOSTICS
# =============================================================================

assess_convergence <- function(mice_object) {
  
  cat("Analyzing convergence patterns...\n")
  
  results <- list()
  
  # Check if chains have converged
  try({
    # Plot convergence (saved to file)
    png("convergence_plots.png", width = 1200, height = 800)
    plot(mice_object, layout = c(2, 3))
    dev.off()
    cat("✓ Convergence plots saved to 'convergence_plots.png'\n")
    
    # Gelman-Rubin diagnostic approximation
    # Check variance between chains vs within chains
    convergence_stats <- list()
    
    for(var in names(mice_object$method)[mice_object$method != ""]) {
      if(var %in% names(mice_object$chainMean)) {
        chain_means <- mice_object$chainMean[[var]]
        if(!is.null(chain_means) && length(chain_means) > 1) {
          between_chain_var <- var(chain_means, na.rm = TRUE)
          # Approximate within-chain variance
          within_chain_var <- mean(apply(mice_object$imp[[var]], 2, var, na.rm = TRUE), na.rm = TRUE)
          
          if(within_chain_var > 0) {
            psrf_approx <- sqrt((between_chain_var + within_chain_var) / within_chain_var)
            convergence_stats[[var]] <- psrf_approx
            
            status <- ifelse(psrf_approx < 1.1, "✓ Converged", 
                             ifelse(psrf_approx < 1.2, "⚠ Marginal", "✗ Poor"))
            cat(paste("  ", var, ": PSRF ≈", round(psrf_approx, 3), status, "\n"))
          }
        }
      }
    }
    
    results$convergence_stats <- convergence_stats
    results$overall_convergence <- ifelse(mean(unlist(convergence_stats), na.rm = TRUE) < 1.1, 
                                          "Good", "Needs more iterations")
    
  }, silent = TRUE)
  
  cat("Convergence assessment:", results$overall_convergence %||% "Unable to assess", "\n")
  
  return(results)
}

# =============================================================================
# 2. DISTRIBUTION PRESERVATION
# =============================================================================

assess_distribution_preservation <- function(original_missing, completed) {
  
  cat("Comparing observed vs imputed value distributions...\n")
  
  results <- list()
  
  for(var in names(original_missing)) {
    if(sum(is.na(original_missing[[var]])) > 0) {
      
      observed_vals <- original_missing[[var]][!is.na(original_missing[[var]])]
      imputed_vals <- completed[[var]][is.na(original_missing[[var]])]
      
      if(is.numeric(observed_vals)) {
        # Numeric variables
        ks_test <- ks.test(observed_vals, imputed_vals)
        
        # Compare moments
        obs_mean <- mean(observed_vals, na.rm = TRUE)
        imp_mean <- mean(imputed_vals, na.rm = TRUE)
        obs_sd <- sd(observed_vals, na.rm = TRUE)
        imp_sd <- sd(imputed_vals, na.rm = TRUE)
        
        mean_ratio <- imp_mean / obs_mean
        sd_ratio <- imp_sd / obs_sd
        
        results[[var]] <- list(
          type = "numeric",
          ks_p_value = ks_test$p.value,
          mean_ratio = mean_ratio,
          sd_ratio = sd_ratio,
          quality = ifelse(ks_test$p.value > 0.05 & abs(mean_ratio - 1) < 0.1 & abs(sd_ratio - 1) < 0.2, 
                           "Good", "Fair")
        )
        
        cat(paste("  ", var, "(numeric): KS p =", round(ks_test$p.value, 4), 
                  ", Mean ratio =", round(mean_ratio, 3), 
                  ", SD ratio =", round(sd_ratio, 3), 
                  "[", results[[var]]$quality, "]\n"))
        
      } else {
        # Categorical variables
        obs_props <- prop.table(table(observed_vals))
        imp_props <- prop.table(table(imputed_vals))
        
        # Chi-square test if possible
        chi_p <- NA
        try({
          common_levels <- intersect(names(obs_props), names(imp_props))
          if(length(common_levels) > 1) {
            chi_test <- chisq.test(table(imputed_vals), p = obs_props[common_levels])
            chi_p <- chi_test$p.value
          }
        }, silent = TRUE)
        
        results[[var]] <- list(
          type = "categorical",
          chi_p_value = chi_p,
          observed_props = obs_props,
          imputed_props = imp_props,
          quality = ifelse(!is.na(chi_p) && chi_p > 0.05, "Good", "Fair")
        )
        
        cat(paste("  ", var, "(categorical): Chi-square p =", 
                  ifelse(is.na(chi_p), "N/A", round(chi_p, 4)), 
                  "[", results[[var]]$quality, "]\n"))
      }
    }
  }
  
  return(results)
}

# =============================================================================
# 3. CORRELATION STRUCTURE PRESERVATION
# =============================================================================

assess_correlation_preservation <- function(original_missing, completed) {
  
  cat("Analyzing correlation structure preservation...\n")
  
  # Get numeric variables only
  numeric_vars <- names(original_missing)[sapply(original_missing, is.numeric)]
  
  if(length(numeric_vars) < 2) {
    cat("  Insufficient numeric variables for correlation analysis\n")
    return(list(quality = "N/A"))
  }
  
  # Calculate correlations for observed data
  observed_data <- original_missing[numeric_vars]
  observed_complete_cases <- observed_data[complete.cases(observed_data), ]
  
  if(nrow(observed_complete_cases) < 10) {
    cat("  Insufficient complete cases for correlation analysis\n")
    return(list(quality = "N/A"))
  }
  
  # Correlations
  cor_observed <- cor(observed_complete_cases, use = "complete.obs")
  cor_imputed <- cor(completed[numeric_vars], use = "complete.obs")
  
  # Compare correlation matrices
  cor_diff <- abs(cor_observed - cor_imputed)
  mean_cor_diff <- mean(cor_diff[upper.tri(cor_diff)], na.rm = TRUE)
  
  # Save correlation plots
  png("correlation_comparison.png", width = 1200, height = 600)
  par(mfrow = c(1, 2))
  corrplot(cor_observed, method = "color", title = "Observed Data Correlations", mar = c(0,0,1,0))
  corrplot(cor_imputed, method = "color", title = "After Imputation Correlations", mar = c(0,0,1,0))
  dev.off()
  
  quality <- ifelse(mean_cor_diff < 0.1, "Excellent",
                    ifelse(mean_cor_diff < 0.2, "Good", 
                           ifelse(mean_cor_diff < 0.3, "Fair", "Poor")))
  
  cat("  Mean absolute correlation difference:", round(mean_cor_diff, 4), "[", quality, "]\n")
  cat("✓ Correlation comparison plots saved to 'correlation_comparison.png'\n")
  
  return(list(
    cor_observed = cor_observed,
    cor_imputed = cor_imputed,
    mean_difference = mean_cor_diff,
    quality = quality
  ))
}

# =============================================================================
# 4. PLAUSIBILITY CHECKS
# =============================================================================

assess_plausibility <- function(original_missing, completed) {
  
  cat("Performing plausibility checks...\n")
  
  results <- list()
  
  for(var in names(original_missing)) {
    if(sum(is.na(original_missing[[var]])) > 0) {
      
      observed_vals <- original_missing[[var]][!is.na(original_missing[[var]])]
      imputed_vals <- completed[[var]][is.na(original_missing[[var]])]
      
      if(is.numeric(observed_vals)) {
        # Check for outliers in imputed values
        obs_range <- range(observed_vals, na.rm = TRUE)
        outliers <- sum(imputed_vals < obs_range[1] | imputed_vals > obs_range[2], na.rm = TRUE)
        outlier_pct <- outliers / length(imputed_vals) * 100
        
        # Check for impossible values (negative where shouldn't be, etc.)
        impossible_count <- 0
        if(min(observed_vals, na.rm = TRUE) >= 0) {
          impossible_count <- sum(imputed_vals < 0, na.rm = TRUE)
        }
        
        results[[var]] <- list(
          type = "numeric",
          outliers_pct = outlier_pct,
          impossible_values = impossible_count,
          quality = ifelse(outlier_pct < 5 && impossible_count == 0, "Good", "Fair")
        )
        
        cat(paste("  ", var, ": Outliers =", round(outlier_pct, 1), "%, Impossible =", 
                  impossible_count, "[", results[[var]]$quality, "]\n"))
        
      } else {
        # Categorical variables - check for new levels
        obs_levels <- levels(as.factor(observed_vals))
        imp_levels <- levels(as.factor(imputed_vals))
        new_levels <- setdiff(imp_levels, obs_levels)
        
        results[[var]] <- list(
          type = "categorical",
          new_levels = new_levels,
          quality = ifelse(length(new_levels) == 0, "Good", "Poor")
        )
        
        cat(paste("  ", var, ": New levels =", length(new_levels), 
                  "[", results[[var]]$quality, "]\n"))
      }
    }
  }
  
  return(results)
}

# =============================================================================
# 5. BETWEEN-IMPUTATION VARIABILITY
# =============================================================================

assess_between_imputation_variability <- function(mice_object) {
  
  cat("Analyzing between-imputation variability...\n")
  
  results <- list()
  
  for(var in names(mice_object$method)[mice_object$method != ""]) {
    if(var %in% names(mice_object$imp)) {
      
      imp_data <- mice_object$imp[[var]]
      
      if(!is.null(imp_data) && ncol(imp_data) > 1) {
        # Calculate coefficient of variation across imputations
        row_means <- rowMeans(imp_data, na.rm = TRUE)
        row_sds <- apply(imp_data, 1, sd, na.rm = TRUE)
        cv <- mean(row_sds / abs(row_means), na.rm = TRUE)
        
        # Lower CV indicates more stable imputations
        stability <- ifelse(cv < 0.1, "Excellent",
                            ifelse(cv < 0.2, "Good",
                                   ifelse(cv < 0.3, "Fair", "Poor")))
        
        results[[var]] <- list(
          coefficient_of_variation = cv,
          stability = stability
        )
        
        cat(paste("  ", var, ": CV =", round(cv, 4), "[", stability, "]\n"))
      }
    }
  }
  
  return(results)
}

# =============================================================================
# 6. CROSS-VALIDATION ASSESSMENT
# =============================================================================

cross_validation_assessment <- function(data_with_missing, k_folds = 5) {
  
  cat("Performing cross-validation assessment...\n")
  
  # This creates additional missing data and tests imputation accuracy
  numeric_vars <- names(data_with_missing)[sapply(data_with_missing, is.numeric)]
  
  if(length(numeric_vars) == 0) {
    cat("  No numeric variables for CV assessment\n")
    return(list(quality = "N/A"))
  }
  
  cv_results <- list()
  
  for(var in numeric_vars) {
    if(sum(!is.na(data_with_missing[[var]])) > 20) {  # Need sufficient observed data
      
      # Get complete cases for this variable
      complete_indices <- which(!is.na(data_with_missing[[var]]))
      
      if(length(complete_indices) > k_folds) {
        
        # Create k-fold splits
        fold_size <- floor(length(complete_indices) / k_folds)
        mae_scores <- numeric(k_folds)
        
        for(fold in 1:k_folds) {
          # Create test set
          test_start <- (fold - 1) * fold_size + 1
          test_end <- min(fold * fold_size, length(complete_indices))
          test_indices <- complete_indices[test_start:test_end]
          
          # Create training data with additional missing values
          train_data <- data_with_missing
          true_values <- train_data[[var]][test_indices]
          train_data[[var]][test_indices] <- NA
          
          # Impute
          try({
            mice_cv <- mice(train_data, m = 1, maxit = 5, printFlag = FALSE)
            completed_cv <- complete(mice_cv)
            predicted_values <- completed_cv[[var]][test_indices]
            
            # Calculate MAE
            mae_scores[fold] <- mean(abs(true_values - predicted_values), na.rm = TRUE)
          }, silent = TRUE)
        }
        
        # Average MAE across folds
        avg_mae <- mean(mae_scores[mae_scores > 0], na.rm = TRUE)
        
        # Normalize by variable scale
        var_range <- diff(range(data_with_missing[[var]], na.rm = TRUE))
        normalized_mae <- avg_mae / var_range
        
        quality <- ifelse(normalized_mae < 0.1, "Excellent",
                          ifelse(normalized_mae < 0.2, "Good",
                                 ifelse(normalized_mae < 0.3, "Fair", "Poor")))
        
        cv_results[[var]] <- list(
          mae = avg_mae,
          normalized_mae = normalized_mae,
          quality = quality
        )
        
        cat(paste("  ", var, ": Normalized MAE =", round(normalized_mae, 4), 
                  "[", quality, "]\n"))
      }
    }
  }
  
  return(cv_results)
}

# =============================================================================
# 7. OVERALL QUALITY SCORE
# =============================================================================

calculate_overall_quality_score <- function(results) {
  
  scores <- list()
  
  # Convergence score
  if(!is.null(results$convergence$overall_convergence)) {
    scores$convergence <- ifelse(results$convergence$overall_convergence == "Good", 1, 0.5)
  }
  
  # Distribution preservation score
  if(!is.null(results$distribution)) {
    dist_scores <- sapply(results$distribution, function(x) {
      ifelse(x$quality == "Good", 1, ifelse(x$quality == "Fair", 0.5, 0))
    })
    scores$distribution <- mean(dist_scores, na.rm = TRUE)
  }
  
  # Correlation preservation score
  if(!is.null(results$correlation$quality)) {
    scores$correlation <- switch(results$correlation$quality,
                                 "Excellent" = 1, "Good" = 0.8, "Fair" = 0.5, "Poor" = 0.2, 0)
  }
  
  # Plausibility score
  if(!is.null(results$plausibility)) {
    plaus_scores <- sapply(results$plausibility, function(x) {
      ifelse(x$quality == "Good", 1, 0.5)
    })
    scores$plausibility <- mean(plaus_scores, na.rm = TRUE)
  }
  
  # Calculate overall score
  overall_score <- mean(unlist(scores), na.rm = TRUE)
  
  grade <- ifelse(overall_score >= 0.8, "A (Excellent)",
                  ifelse(overall_score >= 0.7, "B (Good)",
                         ifelse(overall_score >= 0.6, "C (Fair)", "D (Poor)")))
  
  cat("Overall Quality Score:", round(overall_score, 3), "- Grade:", grade, "\n")
  
  return(list(
    individual_scores = scores,
    overall_score = overall_score,
    grade = grade
  ))
}

# =============================================================================
# 8. SAVE DETAILED REPORT
# =============================================================================

save_validation_report <- function(results, filename) {
  
  sink(filename)
  
  cat("IMPUTATION QUALITY ASSESSMENT REPORT\n")
  cat("=====================================\n")
  cat("Generated on:", Sys.time(), "\n\n")
  
  # Summary
  cat("OVERALL ASSESSMENT\n")
  cat("------------------\n")
  cat("Quality Score:", results$overall_score$overall_score, "\n")
  cat("Grade:", results$overall_score$grade, "\n\n")
  
  # Detailed results for each component
  cat("DETAILED RESULTS\n")
  cat("----------------\n\n")
  
  # Add all detailed results here...
  str(results)
  
  cat("\nRECOMMENDations\n")
  cat("---------------\n")
  
  if(results$overall_score$overall_score >= 0.8) {
    cat("✓ Imputation quality is excellent. Results can be used confidently.\n")
  } else if(results$overall_score$overall_score >= 0.7) {
    cat("✓ Imputation quality is good. Results are reliable for most analyses.\n")
  } else if(results$overall_score$overall_score >= 0.6) {
    cat("⚠ Imputation quality is fair. Consider additional validation or alternative methods.\n")
  } else {
    cat("✗ Imputation quality is poor. Consider different imputation methods or more data.\n")
  }
  
  sink()
  
  cat("✓ Detailed report saved to:", filename, "\n")
}

# =============================================================================
# UTILITY FUNCTIONS
# =============================================================================

`%||%` <- function(x, y) if(is.null(x)) y else x

# =============================================================================
# MAIN WRAPPER FUNCTION
# =============================================================================

assess_imputation_quality <- function(data_with_missing_file, mice_result = NULL, completed_data_file = NULL) {
  
  # Load data
  data_with_missing <- read.csv(data_with_missing_file)
  data_with_missing[data_with_missing == ""] <- NA
  
  # If mice result not provided, perform basic imputation
  if(is.null(mice_result)) {
    cat("No MICE result provided. Performing basic imputation for assessment...\n")
    mice_result <- mice(data_with_missing, m = 5, maxit = 10, printFlag = FALSE)
  }
  
  # Get completed data
  if(is.null(completed_data_file)) {
    completed_data <- complete(mice_result)
  } else {
    completed_data <- read.csv(completed_data_file)
  }
  
  # Run validation
  validation_results <- validate_imputation_no_original(data_with_missing, mice_result, completed_data)
  
  return(validation_results)
}

cat("=============================================================================\n")
cat("IMPUTATION VALIDATION MODULE (NO ORIGINAL DATA) LOADED\n")
cat("=============================================================================\n")
cat("Usage: results <- assess_imputation_quality('data_with_missing.csv')\n")
cat("=============================================================================\n")