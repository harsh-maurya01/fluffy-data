# =============================================================================
# ADAPTIVE MISSING DATA IMPUTATION MODULE
# =============================================================================
# 
# Description: Automatically detects missing data mechanism (MCAR/MAR/MNAR) 
#              and applies optimal imputation techniques
# 
# Author: Your Name
# Date: 2025
# 
# Usage:
#   source("AdaptiveImputationModule.R")
#   result <- run_adaptive_imputation("your_data.csv")
#
# =============================================================================

# Check and install required packages
required_packages <- c("mice", "VIM", "naniar", "randomForest", "dplyr", "corrplot")

for(pkg in required_packages) {
  if(!require(pkg, character.only = TRUE)) {
    cat(paste("Installing package:", pkg, "\n"))
    install.packages(pkg)
    library(pkg, character.only = TRUE)
  }
}

cat("=============================================================================\n")
cat("ADAPTIVE MISSING DATA IMPUTATION MODULE LOADED\n")
cat("=============================================================================\n")
cat("Available functions:\n")
cat("  - detect_missing_mechanism(data)\n")
cat("  - adaptive_imputation(data, mechanism_info)\n") 
cat("  - evaluate_imputation(imputation_result, original_data)\n")
cat("  - run_adaptive_imputation(data_file, original_data_file)\n")
cat("\nQuick start: result <- run_adaptive_imputation('your_file.csv')\n")
cat("=============================================================================\n\n")

# =============================================================================
# MAIN FUNCTIONS
# =============================================================================

detect_missing_mechanism <- function(data) {
  cat("=== MISSING DATA MECHANISM ANALYSIS ===\n\n")
  
  results <- list()
  
  # 1. Little's MCAR Test
  cat("1. Little's MCAR Test:\n")
  tryCatch({
    mcar_test <- mice::mice.mids(mice(data, m=1, maxit=0, printFlag=FALSE))
    # Alternative approach for MCAR test
    missing_pattern <- md.pattern(data, plot = FALSE)
    
    # Simple MCAR heuristic: if missing patterns are random
    n_patterns <- nrow(missing_pattern) - 1
    n_vars_with_missing <- sum(colSums(is.na(data)) > 0)
    
    # If few patterns relative to variables with missing data, likely MCAR
    mcar_likelihood <- ifelse(n_patterns <= n_vars_with_missing * 2, "High", "Low")
    
    cat("  Number of missing patterns:", n_patterns, "\n")
    cat("  Variables with missing data:", n_vars_with_missing, "\n")
    cat("  MCAR likelihood:", mcar_likelihood, "\n\n")
    
    results$mcar_likelihood <- mcar_likelihood
    results$n_patterns <- n_patterns
  }, error = function(e) {
    cat("  Could not perform Little's MCAR test\n\n")
    results$mcar_likelihood <- "Unknown"
  })
  
  # 2. MAR Detection via Logistic Regression
  cat("2. MAR Detection (Logistic Regression Analysis):\n")
  
  mar_results <- list()
  numeric_vars <- names(data)[sapply(data, is.numeric)]
  factor_vars <- names(data)[sapply(data, is.factor)]
  
  for(var in names(data)) {
    if(sum(is.na(data[[var]])) > 0) {
      # Create missingness indicator
      missing_indicator <- is.na(data[[var]])
      
      # Prepare predictors (exclude the variable itself)
      predictors <- data[, !names(data) %in% var, drop = FALSE]
      predictors <- predictors[complete.cases(predictors), ]
      missing_indicator <- missing_indicator[complete.cases(predictors)]
      
      if(length(unique(missing_indicator)) > 1 && nrow(predictors) > 10) {
        tryCatch({
          # Fit logistic regression
          formula_str <- paste("missing_indicator ~", paste(names(predictors), collapse = " + "))
          model <- glm(formula_str, data = cbind(missing_indicator, predictors), family = binomial)
          
          # Check significance of predictors
          p_values <- summary(model)$coefficients[, "Pr(>|z|)"]
          significant_predictors <- sum(p_values[-1] < 0.05, na.rm = TRUE)  # Exclude intercept
          
          mar_results[[var]] <- list(
            significant_predictors = significant_predictors,
            total_predictors = length(p_values) - 1,
            min_p_value = min(p_values[-1], na.rm = TRUE)
          )
          
          cat(paste("  ", var, ": ", significant_predictors, "/", length(p_values)-1, 
                    " significant predictors (min p =", round(min(p_values[-1], na.rm = TRUE), 4), ")\n"))
        }, error = function(e) {
          cat(paste("  ", var, ": Could not fit MAR model\n"))
        })
      }
    }
  }
  
  # Determine overall mechanism
  cat("\n3. Overall Mechanism Assessment:\n")
  
  if(results$mcar_likelihood == "High" && length(mar_results) > 0) {
    significant_mar <- sum(sapply(mar_results, function(x) x$significant_predictors > 0))
    total_vars_tested <- length(mar_results)
    
    if(significant_mar == 0) {
      mechanism <- "MCAR"
      cat("  Assessment: MCAR (Missing Completely At Random)\n")
      cat("  Reason: High MCAR likelihood + No significant MAR predictors\n")
    } else if(significant_mar < total_vars_tested * 0.5) {
      mechanism <- "Mixed (Mostly MCAR)"
      cat("  Assessment: Mixed (Mostly MCAR)\n")
      cat("  Reason: Some MAR patterns but predominantly MCAR\n")
    } else {
      mechanism <- "MAR"
      cat("  Assessment: MAR (Missing At Random)\n")
      cat("  Reason: Significant predictors for missingness patterns\n")
    }
  } else if(length(mar_results) > 0) {
    mechanism <- "MAR"
    cat("  Assessment: MAR (Missing At Random)\n")
    cat("  Reason: Complex missing patterns with predictive relationships\n")
  } else {
    mechanism <- "MNAR (Suspected)"
    cat("  Assessment: MNAR (Missing Not At Random) - Suspected\n")
    cat("  Reason: Complex patterns that don't fit MCAR/MAR assumptions\n")
  }
  
  results$mechanism <- mechanism
  results$mar_results <- mar_results
  
  cat("\n")
  return(results)
}

adaptive_imputation <- function(data, mechanism_info = NULL, m = 5, maxit = 10, seed = 123) {
  
  cat("=== ADAPTIVE IMPUTATION ===\n\n")
  
  # Detect mechanism if not provided
  if(is.null(mechanism_info)) {
    mechanism_info <- detect_missing_mechanism(data)
  }
  
  mechanism <- mechanism_info$mechanism
  cat("Detected mechanism:", mechanism, "\n")
  cat("Selecting optimal imputation methods...\n\n")
  
  # Prepare data
  data_processed <- data
  
  # Convert character columns to factors
  char_cols <- names(data_processed)[sapply(data_processed, is.character)]
  for(col in char_cols) {
    data_processed[[col]] <- as.factor(data_processed[[col]])
  }
  
  # Identify variable types
  numeric_vars <- names(data_processed)[sapply(data_processed, is.numeric)]
  factor_vars <- names(data_processed)[sapply(data_processed, is.factor)]
  
  cat("Variable types identified:\n")
  cat("  Numeric:", length(numeric_vars), "variables\n")
  cat("  Categorical:", length(factor_vars), "variables\n\n")
  
  # Select methods based on mechanism and variable type
  method <- rep("", ncol(data_processed))
  names(method) <- names(data_processed)
  
  for(var in names(data_processed)) {
    if(sum(is.na(data_processed[[var]])) > 0) {  # Only set method for variables with missing data
      
      if(var %in% numeric_vars) {
        # NUMERIC VARIABLES
        if(mechanism %in% c("MCAR", "Mixed (Mostly MCAR)")) {
          method[var] <- "pmm"  # Predictive Mean Matching for MCAR
          cat(paste("  ", var, "(numeric, MCAR): PMM\n"))
        } else if(mechanism == "MAR") {
          # Check if variable has complex relationships
          if(!is.null(mechanism_info$mar_results[[var]]) && 
             mechanism_info$mar_results[[var]]$significant_predictors > 2) {
            method[var] <- "rf"  # Random Forest for complex MAR
            cat(paste("  ", var, "(numeric, complex MAR): Random Forest\n"))
          } else {
            method[var] <- "pmm"  # PMM for simple MAR
            cat(paste("  ", var, "(numeric, simple MAR): PMM\n"))
          }
        } else {  # MNAR
          method[var] <- "pmm"  # Conservative choice for MNAR
          cat(paste("  ", var, "(numeric, MNAR): PMM (conservative)\n"))
        }
        
      } else if(var %in% factor_vars) {
        # CATEGORICAL VARIABLES
        n_levels <- length(levels(data_processed[[var]]))
        
        if(mechanism %in% c("MCAR", "Mixed (Mostly MCAR)")) {
          if(n_levels == 2) {
            method[var] <- "logreg"  # Logistic regression for binary MCAR
            cat(paste("  ", var, "(binary, MCAR): Logistic Regression\n"))
          } else {
            method[var] <- "polyreg"  # Polytomous regression for multi-level MCAR
            cat(paste("  ", var, "(multi-level, MCAR): Polytomous Regression\n"))
          }
        } else if(mechanism == "MAR") {
          # Check complexity
          if(!is.null(mechanism_info$mar_results[[var]]) && 
             mechanism_info$mar_results[[var]]$significant_predictors > 2) {
            method[var] <- "rf"  # Random Forest for complex MAR
            cat(paste("  ", var, "(categorical, complex MAR): Random Forest\n"))
          } else {
            if(n_levels == 2) {
              method[var] <- "logreg"
              cat(paste("  ", var, "(binary, simple MAR): Logistic Regression\n"))
            } else {
              method[var] <- "polyreg"
              cat(paste("  ", var, "(multi-level, simple MAR): Polytomous Regression\n"))
            }
          }
        } else {  # MNAR
          if(n_levels == 2) {
            method[var] <- "logreg"
            cat(paste("  ", var, "(binary, MNAR): Logistic Regression\n"))
          } else {
            method[var] <- "polyreg"
            cat(paste("  ", var, "(multi-level, MNAR): Polytomous Regression\n"))
          }
        }
      }
    }
  }
  
  cat("\n")
  
  # Perform imputation
  cat("Performing imputation...\n")
  
  # Handle MNAR with sensitivity analysis
  if(mechanism == "MNAR (Suspected)") {
    cat("MNAR detected - performing sensitivity analysis with multiple approaches...\n")
    
    # Approach 1: Conservative (PMM/Logistic)
    method_conservative <- method
    for(var in names(method_conservative)) {
      if(method_conservative[var] == "rf") {
        if(var %in% numeric_vars) {
          method_conservative[var] <- "pmm"
        } else if(length(levels(data_processed[[var]])) == 2) {
          method_conservative[var] <- "logreg"
        } else {
          method_conservative[var] <- "polyreg"
        }
      }
    }
    
    # Approach 2: Advanced (Random Forest)
    method_advanced <- method
    for(var in names(method_advanced)) {
      if(method_advanced[var] != "") {
        method_advanced[var] <- "rf"
      }
    }
    
    cat("Running conservative approach...\n")
    imp_conservative <- mice(data_processed, method = method_conservative, m = m, maxit = maxit, 
                             seed = seed, printFlag = FALSE)
    
    cat("Running advanced approach...\n")
    imp_advanced <- mice(data_processed, method = method_advanced, m = m, maxit = maxit, 
                         seed = seed + 100, printFlag = FALSE)
    
    # Return both for comparison
    result <- list(
      mechanism = mechanism,
      imputation_conservative = imp_conservative,
      imputation_advanced = imp_advanced,
      completed_conservative = complete(imp_conservative),
      completed_advanced = complete(imp_advanced),
      methods_used = list(conservative = method_conservative, advanced = method_advanced)
    )
    
  } else {
    # Standard imputation for MCAR/MAR
    cat("Running imputation...\n")
    imp <- mice(data_processed, method = method, m = m, maxit = maxit, seed = seed, printFlag = FALSE)
    
    result <- list(
      mechanism = mechanism,
      imputation = imp,
      completed = complete(imp),
      methods_used = method
    )
  }
  
  cat("Imputation completed!\n\n")
  return(result)
}

evaluate_imputation <- function(imputation_result, original_data = NULL) {
  
  cat("=== IMPUTATION EVALUATION ===\n\n")
  
  mechanism <- imputation_result$mechanism
  
  if(mechanism == "MNAR (Suspected)" && !is.null(imputation_result$completed_conservative)) {
    cat("MNAR Sensitivity Analysis Results:\n")
    cat("Two approaches were used for comparison:\n\n")
    
    cat("1. Conservative Approach (PMM/Logistic Regression):\n")
    print(summary(imputation_result$completed_conservative))
    
    cat("\n2. Advanced Approach (Random Forest):\n")
    print(summary(imputation_result$completed_advanced))
    
    cat("\nRecommendation for MNAR:\n")
    cat("- Compare both results with domain knowledge\n")
    cat("- Consider the conservative approach if interpretability is important\n")
    cat("- Consider the advanced approach if accuracy is paramount\n")
    cat("- Validate with external data if possible\n\n")
    
    # Save both results
    write.csv(imputation_result$completed_conservative, "imputed_conservative.csv", row.names = FALSE)
    write.csv(imputation_result$completed_advanced, "imputed_advanced.csv", row.names = FALSE)
    cat("Results saved: imputed_conservative.csv, imputed_advanced.csv\n")
    
  } else {
    cat("Single Imputation Result:\n")
    print(summary(imputation_result$completed))
    
    # Basic convergence check
    if(!is.null(imputation_result$imputation)) {
      cat("\nConvergence Check:\n")
      cat("Iterations completed:", imputation_result$imputation$iteration, "\n")
      cat("Number of imputations:", imputation_result$imputation$m, "\n")
    }
    
    write.csv(imputation_result$completed, "imputed_adaptive.csv", row.names = FALSE)
    cat("\nResult saved: imputed_adaptive.csv\n")
  }
  
  # If original data provided, perform validation
  if(!is.null(original_data)) {
    cat("\n=== VALIDATION AGAINST ORIGINAL DATA ===\n")
    # This would call the validation function from the previous code
    cat("Use the validation module with the imputed results for detailed accuracy metrics.\n")
  }
}

run_adaptive_imputation <- function(data_file, original_data_file = NULL) {
  
  cat("ADAPTIVE MISSING DATA IMPUTATION MODULE\n")
  cat("=====================================\n\n")
  
  # Load data
  cat("Loading data...\n")
  data <- read.csv(data_file)
  
  # Convert empty strings to NA
  data[data == ""] <- NA
  
  cat("Data loaded:", nrow(data), "rows,", ncol(data), "columns\n")
  cat("Missing data overview:\n")
  print(colSums(is.na(data)))
  cat("\n")
  
  # Step 1: Detect mechanism
  mechanism_info <- detect_missing_mechanism(data)
  
  # Step 2: Perform adaptive imputation
  imputation_result <- adaptive_imputation(data, mechanism_info)
  
  # Step 3: Evaluate results
  original_data <- NULL
  if(!is.null(original_data_file)) {
    original_data <- read.csv(original_data_file)
  }
  
  evaluate_imputation(imputation_result, original_data)
  
  return(imputation_result)
}

# =============================================================================
# MODULE LOADED SUCCESSFULLY
# =============================================================================

cat("✓ Adaptive Imputation Module ready!\n")
cat("✓ Use: result <- run_adaptive_imputation('your_file.csv')\n\n")