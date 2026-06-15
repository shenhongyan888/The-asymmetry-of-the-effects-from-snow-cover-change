# Variance-Covariance Matrix-Based Multilevel Meta-Analysis with Subgroup Comparisons
# Method name: Variance-covariance matrix-based multilevel meta-analysis with subgroup comparisons

# 1. Install and load required packages --------------------------------------------------------
cat("Step 1: Load required packages...\n")

required_packages <- c("metafor", "readxl", "multcomp", "openxlsx", "dplyr", "purrr")
new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages)

library(metafor)
library(readxl)
library(multcomp)
library(openxlsx)
library(dplyr)
library(purrr)

cat("✓ All packages loaded\n\n")

# 2. Data import and preprocessing ---------------------------------------------------------
cat("Step 2: Data import and preprocessing...\n")

# Read data
d2 <- read_excel("C:/Users/Administrator/Desktop/25meta/Meta_effect_sizes_r_2025.12.2_version_0-1.xlsx")

# Function for cleaning trait names
clean_trait_name <- function(trait) {
  trait <- gsub("/", "_", trait)
  trait <- gsub("-", "_", trait)
  trait <- gsub("'", "", trait)
  trait <- gsub("\\s+", "_", trait)
  trait <- gsub("__", "_", trait)
  return(trait)
}

# Define grouping variables
eco_types <- c("Forest", "Wetland", "Grassland", "Tundra")
seasons <- c("Annual average", "Freeze", "Freeze thawing", "Growing season")
frozen_types <- c("P", "NP")

# Updated study-duration grouping: change 0.5-1 years to 1 year
stand_age_groups <- c("1", "2-5", "6-10", "11-15", "16-25")

# Soil-depth grouping  
soil_depth_groups <- c("0-10", "11-30", "31-60", "61-100")

# New: convert continuous variables into grouped variables
cat("  Convert continuous variables into grouped variables...\n")

# Safe data-type conversion
d2$StandAge_year <- as.numeric(as.character(d2$StandAge_year))
d2$Soil_depth <- as.numeric(as.character(d2$Soil_depth))

# Study-duration grouping: use manual grouping to ensure accuracy
d2$StandAge_group <- NA
d2$StandAge_group[d2$StandAge_year >= 0 & d2$StandAge_year <= 1] <- "1"
d2$StandAge_group[d2$StandAge_year > 1 & d2$StandAge_year <= 5] <- "2-5"
d2$StandAge_group[d2$StandAge_year > 5 & d2$StandAge_year <= 10] <- "6-10"
d2$StandAge_group[d2$StandAge_year > 10 & d2$StandAge_year <= 15] <- "11-15"
d2$StandAge_group[d2$StandAge_year > 15 & d2$StandAge_year <= 25] <- "16-25"
d2$StandAge_group[d2$StandAge_year > 25] <- ">25"

# Soil-depth grouping
d2$SoilDepth_group <- cut(d2$Soil_depth,
                          breaks = c(0, 10, 30, 60, 100, Inf),
                          labels = c("0-10", "11-30", "31-60", "61-100", ">100"),
                          include.lowest = TRUE,
                          right = FALSE)

# Check grouping distribution
cat("  Study-duration group distribution:\n")
print(table(d2$StandAge_group, useNA = "always"))

cat("  Soil-depth group distribution:\n")
print(table(d2$SoilDepth_group, useNA = "always"))

# Verify the sample size of the 1-year group
one_year_count <- sum(d2$StandAge_group == "1", na.rm = TRUE)
cat("  Sample size of the 1-year group:", one_year_count, "\n")

cat("✓ Data import and preprocessing completed\n")
cat("  Total data size:", nrow(d2), "rows\n")
cat("  Number of traits:", length(unique(d2$Trait)), "traits\n")
cat("  Grouping types: ecosystem(", length(eco_types), "), season(", length(seasons), "), frozen soil(", length(frozen_types), ")\n")
cat("           study duration(", length(stand_age_groups), "), soil depth(", length(soil_depth_groups), ")\n\n")

# 3. Multiple-comparison function ------------------------------------------------------------
cat("Step 3: Define analysis functions...\n")

run_multiple_comparisons <- function(model, group_type) {
  tryCatch({
    k <- length(coef(model))
    n_comparisons <- k * (k - 1) / 2
    
    if (n_comparisons == 0) return(NULL)
    
    cont_matrix <- matrix(0, nrow = n_comparisons, ncol = k)
    row_index <- 1
    
    for (i in 1:(k-1)) {
      for (j in (i+1):k) {
        cont_matrix[row_index, i] <- 1
        cont_matrix[row_index, j] <- -1
        row_index <- row_index + 1
      }
    }
    
    mc_result <- summary(glht(model, linfct = cont_matrix), 
                         test = adjusted("holm"))
    return(mc_result)
    
  }, error = function(e) {
    return(NULL)
  })
}

# 4. Main analysis function ------------------------------------------------------------
run_comprehensive_meta_analysis <- function(data, trait_name, group_name = "overall") {
  results <- list()
  
  tryCatch({
    if(nrow(data) < 5) return(NULL)
    
    treatment_levels <- unique(na.omit(data$Treatment_1))
    has_both_treatments <- length(treatment_levels) >= 2
    
    # 1. Single-treatment models
    sa_data <- subset(data, `Treatment_1` == "Snow addition")
    if(nrow(sa_data) >= 3) {
      suppressWarnings({
        results$sa <- metafor::rma.mv(yi, V, data = sa_data, 
                                      random = list(~1 | Study_id / Row_id), method = "REML")
      })
    }
    
    sr_data <- subset(data, `Treatment_1` == "Snow remove")
    if(nrow(sr_data) >= 3) {
      suppressWarnings({
        results$sr <- metafor::rma.mv(yi, V, data = sr_data, 
                                      random = list(~1 | Study_id / Row_id), method = "REML")
      })
    }
    
    # 2. Key models: Treatment_1 model for testing differences and Treatment_1 - 1 model for extracting effect sizes
    if(has_both_treatments && nrow(data) >= 5) {
      # Treatment_1 model - test differences
      suppressWarnings({
        results$treatment_model <- metafor::rma.mv(yi, V, 
                                                   mods = ~ `Treatment_1`, 
                                                   data = data,
                                                   random = list(~1 | Study_id / Row_id), 
                                                   method = "REML")
      })
      
      # Treatment_1 - 1 model - extract effect sizes for plotting
      suppressWarnings({
        results$treatment_model_nointercept <- metafor::rma.mv(yi, V, 
                                                               mods = ~ `Treatment_1` - 1, 
                                                               data = data,
                                                               random = list(~1 | Study_id / Row_id), 
                                                               method = "REML")
      })
    }
    
    # 3. Interaction models for multiple comparisons
    if(has_both_treatments && nrow(data) >= 10) {
      
      # Ecosystem interaction
      eco_levels <- unique(na.omit(data$Eco_1))
      if(length(eco_levels) >= 2) {
        suppressWarnings({
          results$eco_interaction <- metafor::rma.mv(yi, V, 
                                                     mods = ~ Treatment_1:Eco_1 - 1,
                                                     data = data,
                                                     random = list(~1 | Study_id / Row_id), 
                                                     method = "REML")
        })
        
        if(!is.null(results$eco_interaction)) {
          results$eco_multcomp <- run_multiple_comparisons(results$eco_interaction, "ecosystem")
        }
      }
      
      # Sampling-season interaction
      season_levels <- unique(na.omit(data$Sampling_season))
      if(length(season_levels) >= 2) {
        suppressWarnings({
          results$season_interaction <- metafor::rma.mv(yi, V, 
                                                        mods = ~ Treatment_1:Sampling_season - 1,
                                                        data = data,
                                                        random = list(~1 | Study_id / Row_id), 
                                                        method = "REML")
        })
        
        if(!is.null(results$season_interaction)) {
          results$season_multcomp <- run_multiple_comparisons(results$season_interaction, "season")
        }
      }
      
      # Frozen-soil-type interaction
      frozen_levels <- unique(na.omit(data$Types_of_frozen_soil))
      frozen_levels <- frozen_levels[!is.na(frozen_levels) & frozen_levels != "NA"]
      if(length(frozen_levels) >= 2) {
        valid_frozen_data <- data[!is.na(data$Types_of_frozen_soil) & data$Types_of_frozen_soil != "NA", ]
        suppressWarnings({
          results$frozen_interaction <- metafor::rma.mv(yi, V, 
                                                        mods = ~ Treatment_1:Types_of_frozen_soil - 1,
                                                        data = valid_frozen_data,
                                                        random = list(~1 | Study_id / Row_id), 
                                                        method = "REML")
        })
        
        if(!is.null(results$frozen_interaction)) {
          results$frozen_multcomp <- run_multiple_comparisons(results$frozen_interaction, "frozen")
        }
      }
    }
    
    if(length(results) > 0) return(results)
    else return(NULL)
    
  }, error = function(e) {
    return(NULL)
  })
}

cat("✓ Analysis functions defined\n\n")

# 5. Run meta-analysis --------------------------------------------------------------
cat("Step 4: Run meta-analysis...\n")

# Get all traits
traits <- unique(d2$Trait)
clean_traits <- sapply(traits, clean_trait_name)

cat("Start analyzing", length(traits), "traits...\n")

# Create a list for storing results
all_results <- list()
success_count <- 0

# Create progress bar
total <- length(traits)
pb <- txtProgressBar(min = 0, max = total, style = 3)

# Main loop
suppressWarnings({
  for(i in seq_along(traits)) {
    trait <- traits[i]
    clean_trait <- clean_trait_name(trait)
    
    setTxtProgressBar(pb, i)
    
    # Filter data for the current trait
    d_current <- subset(d2, Trait == trait)
    
    if(nrow(d_current) < 10) next  # Insufficient data, skip
    
    # Overall analysis
    overall_results <- run_comprehensive_meta_analysis(d_current, clean_trait, "overall")
    if(!is.null(overall_results)) {
      all_results[[clean_trait]] <- list(overall = overall_results)
      success_count <- success_count + 1
    }
  }
})

close(pb)

cat("✓ Meta-analysis completed\n")
cat("  Successfully analyzed:", success_count, "/", length(traits), "traits\n\n")

# 5.5 Group analysis: separate analyses by ecosystem, season, frozen-soil type, study duration, and soil depth
cat("Step 5.5: Run group analysis...\n")

# Revised group-analysis function: now supports all grouping types
run_group_analysis <- function(data, trait_name, group_var, group_values) {
  group_results <- list()
  
  for(group_val in group_values) {
    # Handle NA values
    if(is.na(group_val) || group_val == "NA") {
      group_data <- data[is.na(data[[group_var]]), ]
    } else {
      group_data <- data[data[[group_var]] == group_val & !is.na(data[[group_var]]), ]
    }
    
    if(nrow(group_data) >= 5) {
      # Run the full meta-analysis for each group
      group_analysis <- run_comprehensive_meta_analysis(group_data, trait_name, group_val)
      if(!is.null(group_analysis)) {
        group_results[[as.character(group_val)]] <- group_analysis
      }
    }
  }
  
  if(length(group_results) > 0) return(group_results)
  else return(NULL)
}

# New: Cross-group analysis function
run_cross_group_analysis <- function(data, trait_name, group_var1, group_values1, group_var2, group_values2) {
  cross_results <- list()
  
  for(group1 in group_values1) {
    for(group2 in group_values2) {
      # Filter by the first grouping condition
      if(is.na(group1) || group1 == "NA") {
        data_group1 <- data[is.na(data[[group_var1]]), ]
      } else {
        data_group1 <- data[data[[group_var1]] == group1 & !is.na(data[[group_var1]]), ]
      }
      
      # Filter by the second grouping condition within the first group
      if(is.na(group2) || group2 == "NA") {
        cross_data <- data_group1[is.na(data_group1[[group_var2]]), ]
      } else {
        cross_data <- data_group1[data_group1[[group_var2]] == group2 & !is.na(data_group1[[group_var2]]), ]
      }
      
      if(nrow(cross_data) >= 3) {  # The data requirement for cross-group analysis can be moderately relaxed
        cross_name <- paste(group1, group2, sep = "_")
        cross_analysis <- run_comprehensive_meta_analysis(cross_data, trait_name, cross_name)
        if(!is.null(cross_analysis)) {
          cross_results[[cross_name]] <- cross_analysis
        }
      }
    }
  }
  
  if(length(cross_results) > 0) return(cross_results)
  else return(NULL)
}

# Run group analyses for all traits
cat("Start group analysis...\n")
total_indicators <- length(names(all_results))
pb_groups <- txtProgressBar(min = 0, max = total_indicators, style = 3)

for(i in seq_along(names(all_results))) {
  trait <- names(all_results)[i]
  setTxtProgressBar(pb_groups, i)
  
  # Fix: correctly retrieve the original trait name
  original_trait <- traits[clean_traits == trait][1]
  d_current <- subset(d2, Trait == original_trait)
  
  if(nrow(d_current) >= 10) {
    # Basic group analysis
    eco_results <- run_group_analysis(d_current, trait, "Eco_1", eco_types)
    if(!is.null(eco_results)) all_results[[trait]]$ecosystem_groups <- eco_results
    
    season_results <- run_group_analysis(d_current, trait, "Sampling_season", seasons)
    if(!is.null(season_results)) all_results[[trait]]$season_groups <- season_results
    
    frozen_results <- run_group_analysis(d_current, trait, "Types_of_frozen_soil", frozen_types)
    if(!is.null(frozen_results)) all_results[[trait]]$frozen_groups <- frozen_results
    
    stand_age_results <- run_group_analysis(d_current, trait, "StandAge_group", stand_age_groups)
    if(!is.null(stand_age_results)) all_results[[trait]]$stand_age_groups <- stand_age_results
    
    soil_depth_results <- run_group_analysis(d_current, trait, "SoilDepth_group", soil_depth_groups)
    if(!is.null(soil_depth_results)) all_results[[trait]]$soil_depth_groups <- soil_depth_results
    
    # New: Cross-group analysis
    # Ecosystem x study duration
    eco_stand_age_results <- run_cross_group_analysis(d_current, trait, "Eco_1", eco_types, "StandAge_group", stand_age_groups)
    if(!is.null(eco_stand_age_results)) all_results[[trait]]$eco_stand_age_groups <- eco_stand_age_results
    
    # Ecosystem x soil depth
    eco_soil_depth_results <- run_cross_group_analysis(d_current, trait, "Eco_1", eco_types, "SoilDepth_group", soil_depth_groups)
    if(!is.null(eco_soil_depth_results)) all_results[[trait]]$eco_soil_depth_groups <- eco_soil_depth_results
    
    # Season x study duration
    season_stand_age_results <- run_cross_group_analysis(d_current, trait, "Sampling_season", seasons, "StandAge_group", stand_age_groups)
    if(!is.null(season_stand_age_results)) all_results[[trait]]$season_stand_age_groups <- season_stand_age_results
    
    # Season x soil depth
    season_soil_depth_results <- run_cross_group_analysis(d_current, trait, "Sampling_season", seasons, "SoilDepth_group", soil_depth_groups)
    if(!is.null(season_soil_depth_results)) all_results[[trait]]$season_soil_depth_groups <- season_soil_depth_results
    
    # Frozen soil x study duration
    frozen_stand_age_results <- run_cross_group_analysis(d_current, trait, "Types_of_frozen_soil", frozen_types, "StandAge_group", stand_age_groups)
    if(!is.null(frozen_stand_age_results)) all_results[[trait]]$frozen_stand_age_groups <- frozen_stand_age_results
    
    # Frozen soil x soil depth
    frozen_soil_depth_results <- run_cross_group_analysis(d_current, trait, "Types_of_frozen_soil", frozen_types, "SoilDepth_group", soil_depth_groups)
    if(!is.null(frozen_soil_depth_results)) all_results[[trait]]$frozen_soil_depth_groups <- frozen_soil_depth_results
  }
}

close(pb_groups)
cat("✓ Group analysis completed\n\n")

# 6. Export results ----------------------------------------------------------------
cat("Step 6: Export results...\n")

# Revised data-statistics function
get_analysis_stats <- function(data, trait_name, analysis_type, group_info = NULL) {
  if(is.null(data) || nrow(data) == 0) return(NULL)
  
  # Calculate the actual number of Study_id values included in the analysis after removing duplicates
  study_count <- length(unique(data$Study_id))
  
  # Calculate the actual number of Row_id values included in the analysis
  row_count <- nrow(data)
  
  # Construct statistical information
  stats_info <- data.frame(
    Trait = trait_name,
    Analysis_Type = analysis_type,
    Study_id_Count = study_count,
    Row_id_Count = row_count,
    stringsAsFactors = FALSE
  )
  
  # Add grouping information if available
  if(!is.null(group_info)) {
    stats_info$Group_Type <- group_info$type
    stats_info$Group_Value <- group_info$value
  }
  
  return(stats_info)
}

# ==================== Modified extract_model_summary function ====================
# Modification: redefine the function and add the third argument, stats
extract_model_summary <- function(model, model_name, stats = NULL) {
    if(is.null(model)) return(NULL)
    
    treatment_p <- NA
    if(model_name == "Treatment Model" && !is.null(model$pval) && length(model$pval) > 1) {
        treatment_p <- model$pval[2]
    }
    
    # Create the basic data frame
    result <- data.frame(
        Model = model_name,
        k = model$k,
        Estimate = ifelse(!is.null(model$beta), model$beta[1], NA),
        SE = ifelse(!is.null(model$se), model$se[1], NA),
        z = ifelse(!is.null(model$zval), model$zval[1], NA),
        p_value = ifelse(!is.null(model$pval), model$pval[1], NA),
        Treatment_p_value = treatment_p,
        CI_lower = ifelse(!is.null(model$ci.lb), model$ci.lb[1], NA),
        CI_upper = ifelse(!is.null(model$ci.ub), model$ci.ub[1], NA),
        QE = ifelse(!is.null(model$QE), model$QE, NA),
        QEp = ifelse(!is.null(model$QEp), model$QEp, NA),
        stringsAsFactors = FALSE
    )
    
    # Modification: add statistical information to the results
    if(!is.null(stats)) {
        result$Sample_Size <- ifelse(!is.null(stats$Study_id_Count), stats$Study_id_Count, NA)
        result$Row_Count <- ifelse(!is.null(stats$Row_id_Count), stats$Row_id_Count, NA)
        # Other statistical information fields can be added if needed
    }
    
    return(result)
}

# Function: get a unique worksheet name
get_unique_sheet_name <- function(base_name, existing_names) {
  if (!base_name %in% existing_names) {
    return(base_name)
  }
  
  # If the name already exists, add a numeric suffix
  i <- 1
  new_name <- paste0(base_name, "_", i)
  while (new_name %in% existing_names) {
    i <- i + 1
    new_name <- paste0(base_name, "_", i)
  }
  return(new_name)
}

# Store used worksheet names
used_sheet_names <- c()

# Create Excel workbook
wb <- createWorkbook()

# ==================== Worksheet1: Basic_Model_Summary ====================
cat("  Generate basic model summary...\n")
basic_summary_list <- list()

for(trait in names(all_results)) {
  # Get the original data
  original_trait <- traits[clean_traits == trait][1]
  d_current <- subset(d2, Trait == original_trait)
  
  overall <- all_results[[trait]]$overall
  
  # Snow Addition model
  if(!is.null(overall$sa)) {
    # Get the data actually used for the Snow Addition analysis
    sa_data <- subset(d_current, `Treatment_1` == "Snow addition")
    sa_stats <- get_analysis_stats(sa_data, trait, "Snow Addition")
    
    # Modification: Now the third argument, sa_stats, can be passed correctly
    basic_summary_list[[paste0(trait, "_SnowAddition")]] <- extract_model_summary(
      overall$sa, "Snow Addition", sa_stats)
  }
  
  # Snow Remove model
  if(!is.null(overall$sr)) {
    # Get the data actually used for the Snow Remove analysis
    sr_data <- subset(d_current, `Treatment_1` == "Snow remove")
    sr_stats <- get_analysis_stats(sr_data, trait, "Snow Remove")
    
    # Modification: Now the third argument, sr_stats, can be passed correctly
    basic_summary_list[[paste0(trait, "_SnowRemove")]] <- extract_model_summary(
      overall$sr, "Snow Remove", sr_stats)
  }
  
  # Treatment model(use all data)
  if(!is.null(overall$treatment_model)) {
    treatment_stats <- get_analysis_stats(d_current, trait, "Treatment Comparison")
    
    # Modification: Now the third argument, treatment_stats, can be passed correctly
    basic_summary_list[[paste0(trait, "_TreatmentModel")]] <- extract_model_summary(
      overall$treatment_model, "Treatment Model", treatment_stats)
  }
  
  # Treatment No Intercept model
  if(!is.null(overall$treatment_model_nointercept)) {
    treatment_noint_stats <- get_analysis_stats(d_current, trait, "Treatment No Intercept")
    
    # Modification: Now the third argument, treatment_noint_stats, can be passed correctly
    basic_summary_list[[paste0(trait, "_TreatmentNoIntercept")]] <- extract_model_summary(
      overall$treatment_model_nointercept, "Treatment No Intercept", treatment_noint_stats)
  }
}

if(length(basic_summary_list) > 0) {
  basic_summary_df <- do.call(rbind, lapply(basic_summary_list, function(x) {
    if(!is.null(x)) return(x)
  }))
  
  if(!is.null(basic_summary_df)) {
    basic_summary_df$Trait <- gsub("_.*", "", rownames(basic_summary_df))
    basic_summary_df$ModelType <- gsub(".*_", "", rownames(basic_summary_df))
    rownames(basic_summary_df) <- NULL
    
    sheet_name <- get_unique_sheet_name("Basic_Model_Summary", used_sheet_names)
    used_sheet_names <- c(used_sheet_names, sheet_name)
    addWorksheet(wb, sheet_name)
    writeData(wb, sheet_name, basic_summary_df)
    cat("  Generated worksheet:", sheet_name, "\n")
  }
}

# ==================== Worksheet2: Treatment_Effect_Sizes ====================
cat("  Generate Treatment effect sizes...\n")
treatment_coefs_list <- list()

for(trait in names(all_results)) {
  # Get the original data for statistics
  original_trait <- traits[clean_traits == trait][1]
  d_current <- subset(d2, Trait == original_trait)
  overall_stats <- get_analysis_stats(d_current, trait, "Treatment No Intercept")
  
  overall <- all_results[[trait]]$overall
  
  if(!is.null(overall$treatment_model_nointercept)) {
    coefs <- coef(overall$treatment_model_nointercept)
    ses <- sqrt(diag(vcov(overall$treatment_model_nointercept)))
    ci_lb <- overall$treatment_model_nointercept$ci.lb
    ci_ub <- overall$treatment_model_nointercept$ci.ub
    
    for(i in 1:length(coefs)) {
      treatment_coefs_list[[paste0(trait, "_", names(coefs)[i])]] <- data.frame(
        Trait = trait,
        Coefficient = names(coefs)[i],
        Estimate = coefs[i],
        SE = ses[i],
        CI_lower = ci_lb[i],
        CI_upper = ci_ub[i],
        z_value = coefs[i]/ses[i],
        p_value = 2*pnorm(abs(coefs[i]/ses[i]), lower.tail = FALSE),
        Study_id_Count = overall_stats$Study_id_Count,
        Row_id_Count = overall_stats$Row_id_Count,
        stringsAsFactors = FALSE
      )
    }
  }
}

if(length(treatment_coefs_list) > 0) {
  treatment_coefs_df <- do.call(rbind, treatment_coefs_list)
  rownames(treatment_coefs_df) <- NULL
  
  sheet_name <- get_unique_sheet_name("Treatment_Effect_Sizes", used_sheet_names)
  used_sheet_names <- c(used_sheet_names, sheet_name)
  addWorksheet(wb, sheet_name)
  writeData(wb, sheet_name, treatment_coefs_df)
  cat("  Generated worksheet:", sheet_name, "\n")
}

# ==================== Worksheet3: Interaction_Coefficients ====================
cat("  Generate interaction coefficients...\n")
interaction_list <- list()

for(trait in names(all_results)) {
  # Get the original data for statistics
  original_trait <- traits[clean_traits == trait][1]
  d_current <- subset(d2, Trait == original_trait)
  
  overall <- all_results[[trait]]$overall
  
  if(!is.null(overall$eco_interaction)) {
    # Get the actual data for the interaction model
    eco_int_stats <- get_analysis_stats(d_current, trait, "Ecosystem Interaction")
    
    coefs <- coef(overall$eco_interaction)
    ses <- sqrt(diag(vcov(overall$eco_interaction)))
    
    for(i in 1:length(coefs)) {
      interaction_list[[paste0(trait, "_Eco_", i)]] <- data.frame(
        Trait = trait,
        Model = "Ecosystem Interaction",
        Coefficient = names(coefs)[i],
        Estimate = coefs[i],
        SE = ses[i],
        z_value = coefs[i]/ses[i],
        p_value = 2*pnorm(abs(coefs[i]/ses[i]), lower.tail = FALSE),
        Study_id_Count = eco_int_stats$Study_id_Count,
        Row_id_Count = eco_int_stats$Row_id_Count,
        stringsAsFactors = FALSE
      )
    }
  }
  
  if(!is.null(overall$season_interaction)) {
    # Get the actual data for the interaction model
    season_int_stats <- get_analysis_stats(d_current, trait, "Season Interaction")
    
    coefs <- coef(overall$season_interaction)
    ses <- sqrt(diag(vcov(overall$season_interaction)))
    
    for(i in 1:length(coefs)) {
      interaction_list[[paste0(trait, "_Season_", i)]] <- data.frame(
        Trait = trait,
        Model = "Season Interaction",
        Coefficient = names(coefs)[i],
        Estimate = coefs[i],
        SE = ses[i],
        z_value = coefs[i]/ses[i],
        p_value = 2*pnorm(abs(coefs[i]/ses[i]), lower.tail = FALSE),
        Study_id_Count = season_int_stats$Study_id_Count,
        Row_id_Count = season_int_stats$Row_id_Count,
        stringsAsFactors = FALSE
      )
    }
  }
}

if(length(interaction_list) > 0) {
  interaction_df <- do.call(rbind, interaction_list)
  rownames(interaction_df) <- NULL
  
  sheet_name <- get_unique_sheet_name("Interaction_Coefficients", used_sheet_names)
  used_sheet_names <- c(used_sheet_names, sheet_name)
  addWorksheet(wb, sheet_name)
  writeData(wb, sheet_name, interaction_df)
  cat("  Generated worksheet:", sheet_name, "\n")
}

# ==================== Worksheet4: Multiple_Comparison_Results ====================
cat("  Generate multiple-comparison results...\n")
multcomp_list <- list()

for(trait in names(all_results)) {
  # Get the original data for statistics
  original_trait <- traits[clean_traits == trait][1]
  d_current <- subset(d2, Trait == original_trait)
  
  overall <- all_results[[trait]]$overall
  
  if(!is.null(overall$eco_multcomp)) {
    # Get the actual data for multiple comparisons
    eco_multcomp_stats <- get_analysis_stats(d_current, trait, "Ecosystem Multiple Comparison")
    
    tryCatch({
      coefs <- overall$eco_multcomp$test$coefficients
      ses <- overall$eco_multcomp$test$sigma
      tstats <- overall$eco_multcomp$test$tstat
      pvals <- overall$eco_multcomp$test$pvalues
      
      for(i in 1:length(coefs)) {
        multcomp_list[[paste0(trait, "_Eco_", i)]] <- data.frame(
          Trait = trait,
          Comparison = "Ecosystem",
          Contrast = ifelse(!is.null(names(coefs)[i]), names(coefs)[i], paste0("Contrast_", i)),
          Estimate = coefs[i],
          SE = ses[i],
          z_value = tstats[i],
          p_value = pvals[i],
          Study_id_Count = eco_multcomp_stats$Study_id_Count,
          Row_id_Count = eco_multcomp_stats$Row_id_Count,
          stringsAsFactors = FALSE
        )
      }
    }, error = function(e) {
      cat("  Processing", trait, "'s ecosystem multiple comparison failed:", e$message, "\n")
    })
  }
}

if(length(multcomp_list) > 0) {
  multcomp_df <- do.call(rbind, multcomp_list)
  rownames(multcomp_df) <- NULL
  
  sheet_name <- get_unique_sheet_name("Multiple_Comparisons", used_sheet_names)
  used_sheet_names <- c(used_sheet_names, sheet_name)
  addWorksheet(wb, sheet_name)
  writeData(wb, sheet_name, multcomp_df)
  cat("  Generated worksheet:", sheet_name, "\n")
}

# ==================== Worksheet5: Ecosystem_Group_Effect_Sizes ====================
cat("  Generate ecosystem-group effect sizes...\n")
eco_group_list <- list()

for(trait in names(all_results)) {
  if(!is.null(all_results[[trait]]$ecosystem_groups)) {
    # Get the original data
    original_trait <- traits[clean_traits == trait][1]
    d_current <- subset(d2, Trait == original_trait)
    
    for(eco in names(all_results[[trait]]$ecosystem_groups)) {
      # Get the actual data for this ecosystem
      eco_data <- d_current[d_current$Eco_1 == eco & !is.na(d_current$Eco_1), ]
      
      group_results <- all_results[[trait]]$ecosystem_groups[[eco]]
      
      # Extract Treatment effect sizes for this ecosystem
      if(!is.null(group_results$treatment_model_nointercept)) {
        # Get the actual data for the Treatment model in this group
        eco_treatment_stats <- get_analysis_stats(eco_data, trait, "Treatment No Intercept", 
                                                 list(type = "Ecosystem", value = eco))
        
        coefs <- coef(group_results$treatment_model_nointercept)
        ses <- sqrt(diag(vcov(group_results$treatment_model_nointercept)))
        ci_lb <- group_results$treatment_model_nointercept$ci.lb
        ci_ub <- group_results$treatment_model_nointercept$ci.ub
        
        for(i in 1:length(coefs)) {
          eco_group_list[[paste0(trait, "_", eco, "_", i)]] <- data.frame(
            Trait = trait,
            Ecosystem = eco,
            Coefficient = names(coefs)[i],
            Estimate = coefs[i],
            SE = ses[i],
            CI_lower = ci_lb[i],
            CI_upper = ci_ub[i],
            z_value = coefs[i]/ses[i],
            p_value = 2*pnorm(abs(coefs[i]/ses[i]), lower.tail = FALSE),
            Study_id_Count = eco_treatment_stats$Study_id_Count,
            Row_id_Count = eco_treatment_stats$Row_id_Count,
            stringsAsFactors = FALSE
          )
        }
      }
    }
  }
}

if(length(eco_group_list) > 0) {
  eco_group_df <- do.call(rbind, eco_group_list)
  rownames(eco_group_df) <- NULL
  
  sheet_name <- get_unique_sheet_name("Ecosystem_Group_Effect_Sizes", used_sheet_names)
  used_sheet_names <- c(used_sheet_names, sheet_name)
  addWorksheet(wb, sheet_name)
  writeData(wb, sheet_name, eco_group_df)
  cat("  Generated worksheet:", sheet_name, "\n")
}

# ==================== Worksheet6: Season_Group_Effect_Sizes ====================
cat("  Generate season-group effect sizes...\n")
season_group_list <- list()

for(trait in names(all_results)) {
  if(!is.null(all_results[[trait]]$season_groups)) {
    # Get the original data
    original_trait <- traits[clean_traits == trait][1]
    d_current <- subset(d2, Trait == original_trait)
    
    for(season in names(all_results[[trait]]$season_groups)) {
      # Get the actual data for this season
      season_data <- d_current[d_current$Sampling_season == season & !is.na(d_current$Sampling_season), ]
      
      group_results <- all_results[[trait]]$season_groups[[season]]
      
      # Extract Treatment effect sizes for this season
      if(!is.null(group_results$treatment_model_nointercept)) {
        # Get the actual data for the Treatment model in this group
        season_treatment_stats <- get_analysis_stats(season_data, trait, "Treatment No Intercept", 
                                                    list(type = "Season", value = season))
        
        coefs <- coef(group_results$treatment_model_nointercept)
        ses <- sqrt(diag(vcov(group_results$treatment_model_nointercept)))
        ci_lb <- group_results$treatment_model_nointercept$ci.lb
        ci_ub <- group_results$treatment_model_nointercept$ci.ub
        
        for(i in 1:length(coefs)) {
          season_group_list[[paste0(trait, "_", season, "_", i)]] <- data.frame(
            Trait = trait,
            Season = season,
            Coefficient = names(coefs)[i],
            Estimate = coefs[i],
            SE = ses[i],
            CI_lower = ci_lb[i],
            CI_upper = ci_ub[i],
            z_value = coefs[i]/ses[i],
            p_value = 2*pnorm(abs(coefs[i]/ses[i]), lower.tail = FALSE),
            Study_id_Count = season_treatment_stats$Study_id_Count,
            Row_id_Count = season_treatment_stats$Row_id_Count,
            stringsAsFactors = FALSE
          )
        }
      }
    }
  }
}

if(length(season_group_list) > 0) {
  season_group_df <- do.call(rbind, season_group_list)
  rownames(season_group_df) <- NULL
  
  sheet_name <- get_unique_sheet_name("Season_Group_Effect_Sizes", used_sheet_names)
  used_sheet_names <- c(used_sheet_names, sheet_name)
  addWorksheet(wb, sheet_name)
  writeData(wb, sheet_name, season_group_df)
  cat("  Generated worksheet:", sheet_name, "\n")
}

# ==================== Worksheet7: Frozen_Soil_Group_Effect_Sizes ====================
cat("  Generate frozen-soil-group effect sizes...\n")
frozen_group_list <- list()

for(trait in names(all_results)) {
  if(!is.null(all_results[[trait]]$frozen_groups)) {
    # Get the original data
    original_trait <- traits[clean_traits == trait][1]
    d_current <- subset(d2, Trait == original_trait)
    
    for(frozen in names(all_results[[trait]]$frozen_groups)) {
      # Get the actual data for this frozen-soil type
      frozen_data <- d_current[d_current$Types_of_frozen_soil == frozen & !is.na(d_current$Types_of_frozen_soil), ]
      frozen_name <- ifelse(frozen == "P", "Permafrost", "Non_permafrost")
      
      group_results <- all_results[[trait]]$frozen_groups[[frozen]]
      
      # Extract Treatment effect sizes for this frozen-soil type
      if(!is.null(group_results$treatment_model_nointercept)) {
        # Get the actual data for the Treatment model in this group
        frozen_treatment_stats <- get_analysis_stats(frozen_data, trait, "Treatment No Intercept", 
                                                    list(type = "Frozen_Soil", value = frozen_name))
        
        coefs <- coef(group_results$treatment_model_nointercept)
        ses <- sqrt(diag(vcov(group_results$treatment_model_nointercept)))
        ci_lb <- group_results$treatment_model_nointercept$ci.lb
        ci_ub <- group_results$treatment_model_nointercept$ci.ub
        
        for(i in 1:length(coefs)) {
          frozen_group_list[[paste0(trait, "_", frozen, "_", i)]] <- data.frame(
            Trait = trait,
            Frozen_Type = frozen_name,
            Coefficient = names(coefs)[i],
            Estimate = coefs[i],
            SE = ses[i],
            CI_lower = ci_lb[i],
            CI_upper = ci_ub[i],
            z_value = coefs[i]/ses[i],
            p_value = 2*pnorm(abs(coefs[i]/ses[i]), lower.tail = FALSE),
            Study_id_Count = frozen_treatment_stats$Study_id_Count,
            Row_id_Count = frozen_treatment_stats$Row_id_Count,
            stringsAsFactors = FALSE
          )
        }
      }
    }
  }
}

if(length(frozen_group_list) > 0) {
  frozen_group_df <- do.call(rbind, frozen_group_list)
  rownames(frozen_group_df) <- NULL
  
  sheet_name <- get_unique_sheet_name("Frozen_Soil_Group_Effect_Sizes", used_sheet_names)
  used_sheet_names <- c(used_sheet_names, sheet_name)
  addWorksheet(wb, sheet_name)
  writeData(wb, sheet_name, frozen_group_df)
  cat("  Generated worksheet:", sheet_name, "\n")
}

# ==================== Worksheet8: Study_Duration_Group_Effect_Sizes ====================
cat("  Generate study-duration-group effect sizes...\n")
stand_age_group_list <- list()

for(trait in names(all_results)) {
  if(!is.null(all_results[[trait]]$stand_age_groups)) {
    # Get the original data
    original_trait <- traits[clean_traits == trait][1]
    d_current <- subset(d2, Trait == original_trait)
    
    for(age_group in names(all_results[[trait]]$stand_age_groups)) {
      # Get the actual data for this group
      age_data <- d_current[d_current$StandAge_group == age_group & !is.na(d_current$StandAge_group), ]
      
      group_results <- all_results[[trait]]$stand_age_groups[[age_group]]
      
      # Extract Treatment effect sizes for this group
      if(!is.null(group_results$treatment_model_nointercept)) {
        # Get the actual data for the Treatment model in this group
        age_treatment_stats <- get_analysis_stats(age_data, trait, "Treatment No Intercept", 
                                                 list(type = "StandAge", value = paste0(age_group, " years")))
        
        coefs <- coef(group_results$treatment_model_nointercept)
        ses <- sqrt(diag(vcov(group_results$treatment_model_nointercept)))
        ci_lb <- group_results$treatment_model_nointercept$ci.lb
        ci_ub <- group_results$treatment_model_nointercept$ci.ub
        
        for(i in 1:length(coefs)) {
          stand_age_group_list[[paste0(trait, "_", age_group, "_", i)]] <- data.frame(
            Trait = trait,
            StandAge_Group = paste0(age_group, " years"),
            Coefficient = names(coefs)[i],
            Estimate = coefs[i],
            SE = ses[i],
            CI_lower = ci_lb[i],
            CI_upper = ci_ub[i],
            z_value = coefs[i]/ses[i],
            p_value = 2*pnorm(abs(coefs[i]/ses[i]), lower.tail = FALSE),
            Study_id_Count = age_treatment_stats$Study_id_Count,
            Row_id_Count = age_treatment_stats$Row_id_Count,
            stringsAsFactors = FALSE
          )
        }
      }
    }
  }
}

if(length(stand_age_group_list) > 0) {
  stand_age_group_df <- do.call(rbind, stand_age_group_list)
  rownames(stand_age_group_df) <- NULL
  
  sheet_name <- get_unique_sheet_name("Study_Duration_Group_Effect_Sizes", used_sheet_names)
  used_sheet_names <- c(used_sheet_names, sheet_name)
  addWorksheet(wb, sheet_name)
  writeData(wb, sheet_name, stand_age_group_df)
  cat("  Generated worksheet:", sheet_name, "\n")
}

# ==================== Worksheet9: Soil_Depth_Group_Effect_Sizes ====================
cat("  Generate soil-depth-group effect sizes...\n")
soil_depth_group_list <- list()

for(trait in names(all_results)) {
  if(!is.null(all_results[[trait]]$soil_depth_groups)) {
    # Get the original data
    original_trait <- traits[clean_traits == trait][1]
    d_current <- subset(d2, Trait == original_trait)
    
    for(depth_group in names(all_results[[trait]]$soil_depth_groups)) {
      # Get the actual data for this group
      depth_data <- d_current[d_current$SoilDepth_group == depth_group & !is.na(d_current$SoilDepth_group), ]
      
      group_results <- all_results[[trait]]$soil_depth_groups[[depth_group]]
      
      # Extract Treatment effect sizes for this group
      if(!is.null(group_results$treatment_model_nointercept)) {
        # Get the actual data for the Treatment model in this group
        depth_treatment_stats <- get_analysis_stats(depth_data, trait, "Treatment No Intercept", 
                                                   list(type = "SoilDepth", value = paste0(depth_group, " cm")))
        
        coefs <- coef(group_results$treatment_model_nointercept)
        ses <- sqrt(diag(vcov(group_results$treatment_model_nointercept)))
        ci_lb <- group_results$treatment_model_nointercept$ci.lb
        ci_ub <- group_results$treatment_model_nointercept$ci.ub
        
        for(i in 1:length(coefs)) {
          soil_depth_group_list[[paste0(trait, "_", depth_group, "_", i)]] <- data.frame(
            Trait = trait,
            SoilDepth_Group = paste0(depth_group, " cm"),
            Coefficient = names(coefs)[i],
            Estimate = coefs[i],
            SE = ses[i],
            CI_lower = ci_lb[i],
            CI_upper = ci_ub[i],
            z_value = coefs[i]/ses[i],
            p_value = 2*pnorm(abs(coefs[i]/ses[i]), lower.tail = FALSE),
            Study_id_Count = depth_treatment_stats$Study_id_Count,
            Row_id_Count = depth_treatment_stats$Row_id_Count,
            stringsAsFactors = FALSE
          )
        }
      }
    }
  }
}

if(length(soil_depth_group_list) > 0) {
  soil_depth_group_df <- do.call(rbind, soil_depth_group_list)
  rownames(soil_depth_group_df) <- NULL
  
  sheet_name <- get_unique_sheet_name("Soil_Depth_Group_Effect_Sizes", used_sheet_names)
  used_sheet_names <- c(used_sheet_names, sheet_name)
  addWorksheet(wb, sheet_name)
  writeData(wb, sheet_name, soil_depth_group_df)
  cat("  Generated worksheet:", sheet_name, "\n")
}

# ==================== Worksheet10: Group comparison summaries in multiple worksheets ====================
# Define grouping variables
eco_types <- c("Forest", "Wetland", "Grassland", "Tundra")
seasons <- c("Annual average", "Freeze", "Freeze thawing", "Growing season")
frozen_types <- c("P", "NP")
stand_age_groups <- c("1", "2-5", "6-10", "11-15", "15-25")
soil_depth_groups <- c("0-10", "11-30", "31-60", "61-100")
# Define a function for extracting group-specific treatment differences with two levels
extract_group_comparison <- function(trait, trait_results, d_current, group_type, 
                                     group_var, group_values, groups_list_name = NULL) {
  comparison_list <- list()
  
  # ==================== Level 1: overall treatment difference ====================
  if(group_type == "Overall") {
    if(!is.null(trait_results$overall$treatment_model)) {
      overall_model <- trait_results$overall$treatment_model
      overall_stats <- get_analysis_stats(d_current, trait, "Treatment Comparison")
      
      comparison_list[[paste0(trait, "_Overall")]] <- data.frame(
        Trait = trait,
        Group = "Overall",
        Group_Type = "Overall",
        Analysis_Level = "Overall",
        QM = overall_model$QM,
        df1 = overall_model$m,
        p_value = ifelse(length(overall_model$pval) > 1, overall_model$pval[2], overall_model$pval[1]),
        k = overall_model$k,
        Study_id_Count = overall_stats$Study_id_Count,
        Row_id_Count = overall_stats$Row_id_Count,
        stringsAsFactors = FALSE
      )
    }
  } 
  # ==================== Level 1: QM value with the grouping variable as a moderator ====================
  else if(group_type != "Overall") {
    # Run a meta-regression model using the grouping variable as a moderator
    # That is: yi ~ Treatment_1 * Group_Variable
    tryCatch({
      # Ensure that the grouping variable has sufficient data
      valid_data <- d_current[!is.na(d_current[[group_var]]), ]
      unique_groups <- unique(valid_data[[group_var]])
      
      if(length(unique_groups) >= 2 && nrow(valid_data) >= 10) {
        # Run the interaction model
        suppressWarnings({
          group_interaction_model <- metafor::rma.mv(
            yi, V,
            mods = ~ `Treatment_1` * factor(get(group_var)) - 1,
            data = valid_data,
            random = list(~1 | Study_id / Row_id),
            method = "REML"
          )
        })
        
        if(!is.null(group_interaction_model)) {
          group_stats <- get_analysis_stats(valid_data, trait, paste0(group_type, " Interaction"))
          
          # Extract the QM value for the overall test of treatment effects
          comparison_list[[paste0(trait, "_", group_type, "_GroupEffect")]] <- data.frame(
            Trait = trait,
            Group = paste0("All ", group_type, " groups"),
            Group_Type = group_type,
            Analysis_Level = "Group_Effect",
            QM = group_interaction_model$QM,
            df1 = group_interaction_model$m,
            p_value = ifelse(length(group_interaction_model$pval) > 1, 
                           group_interaction_model$pval[2], 
                           group_interaction_model$pval[1]),
            k = group_interaction_model$k,
            Study_id_Count = group_stats$Study_id_Count,
            Row_id_Count = group_stats$Row_id_Count,
            stringsAsFactors = FALSE
          )
        }
      }
    }, error = function(e) {
      # If the model fails, skip
    })
  }
  
  # ==================== Level 2: QM value within each subgroup ====================
  if(group_type != "Overall" && !is.null(groups_list_name) && !is.null(trait_results[[groups_list_name]])) {
    actual_groups <- names(trait_results[[groups_list_name]])
    
    for(group_val in actual_groups) {
      if(is.na(group_val) || group_val == "" || group_val == "NA") next
      
      # Get the data for this group
      if(group_type == "Ecosystem") {
        group_data <- d_current[d_current$Eco_1 == group_val & !is.na(d_current$Eco_1), ]
      } else if(group_type == "Season") {
        group_data <- d_current[d_current$Sampling_season == group_val & !is.na(d_current$Sampling_season), ]
      } else if(group_type == "Frozen_Soil") {
        group_data <- d_current[d_current$Types_of_frozen_soil == group_val & !is.na(d_current$Types_of_frozen_soil), ]
      } else if(group_type == "StandAge") {
        group_data <- d_current[d_current$StandAge_group == group_val & !is.na(d_current$StandAge_group), ]
      } else if(group_type == "SoilDepth") {
        group_data <- d_current[d_current$SoilDepth_group == group_val & !is.na(d_current$SoilDepth_group), ]
      } else {
        group_data <- d_current[d_current[[group_var]] == group_val & !is.na(d_current[[group_var]]), ]
      }
      
      if(nrow(group_data) < 3) next
      
      group_stats <- get_analysis_stats(group_data, trait, "Treatment Comparison", 
                                        list(type = group_type, value = group_val))
      
      group_results <- trait_results[[groups_list_name]][[group_val]]
      
      if(!is.null(group_results$treatment_model)) {
        model <- group_results$treatment_model
        comparison_list[[paste0(trait, "_", group_type, "_Subgroup_", group_val)]] <- data.frame(
          Trait = trait,
          Group = group_val,
          Group_Type = group_type,
          Analysis_Level = "Within_Subgroup",
          QM = model$QM,
          df1 = model$m,
          p_value = ifelse(length(model$pval) > 1, model$pval[2], model$pval[1]),
          k = model$k,
          Study_id_Count = group_stats$Study_id_Count,
          Row_id_Count = group_stats$Row_id_Count,
          stringsAsFactors = FALSE
        )
      }
    }
  }
  
  return(comparison_list)}
# Define group-type configurations
group_configs <- list(
  list(
    type = "Overall", 
    var = NULL, 
    name = "Overall_Comparison_Summary",
    groups_list = NULL,
    values = NULL
  ),
  list(
    type = "Ecosystem", 
    var = "Eco_1", 
    name = "Ecosystem_Group_Comparison_Summary",
    groups_list = "ecosystem_groups",
    values = eco_types
  ),
  list(
    type = "Season", 
    var = "Sampling_season", 
    name = "Sampling_Season_Group_Comparison_Summary",
    groups_list = "season_groups",
    values = seasons
  ),
  list(
    type = "SoilDepth", 
    var = "SoilDepth_group", 
    name = "Soil_Layer_Depth_Group_Comparison_Summary",
    groups_list = "soil_depth_groups",
    values = soil_depth_groups
  ),
  list(
    type = "Frozen_Soil", 
    var = "Types_of_frozen_soil", 
    name = "Frozen_Soil_Type_Group_Comparison_Summary",
    groups_list = "frozen_groups",
    values = frozen_types
  ),
  list(
    type = "StandAge", 
    var = "StandAge_group", 
    name = "Study_Duration_Group_Comparison_Summary",
    groups_list = "stand_age_groups",
    values = stand_age_groups
  ))
# Create a separate worksheet for each group type
for(config in group_configs) {
  cat("  Generate", config$name, "...\n")
  
  group_comparison_list <- list()
  
  for(trait in names(all_results)) {
    original_trait <- traits[clean_traits == trait][1]
    d_current <- subset(d2, Trait == original_trait)
    
    comparisons <- extract_group_comparison(
      trait = trait,
      trait_results = all_results[[trait]],
      d_current = d_current,
      group_type = config$type,
      group_var = config$var,
      group_values = config$values,
      groups_list_name = config$groups_list
    )
    
    if(length(comparisons) > 0) {
      group_comparison_list <- c(group_comparison_list, comparisons)
    }
  }
  
  if(length(group_comparison_list) > 0) {
    group_comparison_df <- do.call(rbind, group_comparison_list)
    rownames(group_comparison_df) <- NULL
    
    # Sort by analysis level first and then by group
    group_comparison_df <- group_comparison_df[
      order(
        group_comparison_df$Trait,
        factor(group_comparison_df$Analysis_Level, 
               levels = c("Overall", "Group_Effect", "Within_Subgroup")),
        group_comparison_df$Group
      ), 
    ]
    
    sheet_name <- get_unique_sheet_name(config$name, used_sheet_names)
    used_sheet_names <- c(used_sheet_names, sheet_name)
    addWorksheet(wb, sheet_name)
    writeData(wb, sheet_name, group_comparison_df)
    cat("  Generated worksheet:", sheet_name, "(", nrow(group_comparison_df), "rows)\n")
  } else {
    cat("  Warning:", config$name, "has no valid results, skip\n")
  }
}
# ==================== Summarize QM values for all group effects ====================
cat("  Generate group-effect summary worksheet...\n")
group_effect_summary_list <- list()
for(trait in names(all_results)) {
  original_trait <- traits[clean_traits == trait][1]
  d_current <- subset(d2, Trait == original_trait)
  
  # Collect effect QM values for all group types
  for(config in group_configs) {
    if(config$type == "Overall") next
    
    # Try to run the group interaction model
    tryCatch({
      valid_data <- d_current[!is.na(d_current[[config$var]]), ]
      unique_groups <- unique(valid_data[[config$var]])
      
      if(length(unique_groups) >= 2 && nrow(valid_data) >= 10) {
        suppressWarnings({
          group_interaction_model <- metafor::rma.mv(
            yi, V,
            mods = ~ `Treatment_1` * factor(get(config$var)) - 1,
            data = valid_data,
            random = list(~1 | Study_id / Row_id),
            method = "REML"
          )
        })
        
        if(!is.null(group_interaction_model)) {
          group_stats <- get_analysis_stats(valid_data, trait, paste0(config$type, " Interaction"))
          
          group_effect_summary_list[[paste0(trait, "_", config$type)]] <- data.frame(
            Trait = trait,
            Group_Type = config$type,
            QM = group_interaction_model$QM,
            df1 = group_interaction_model$m,
            p_value = ifelse(length(group_interaction_model$pval) > 1, 
                           group_interaction_model$pval[2], 
                           group_interaction_model$pval[1]),
            k = group_interaction_model$k,
            n_Groups = length(unique_groups),
            Groups = paste(unique_groups, collapse = ", "),
            Study_id_Count = group_stats$Study_id_Count,
            Row_id_Count = group_stats$Row_id_Count,
            stringsAsFactors = FALSE
          )
        }
      }
    }, error = function(e) {
      # Skip errors
    })
  }
}
# Create the group-effect summary worksheet
if(length(group_effect_summary_list) > 0) {
  group_effect_summary_df <- do.call(rbind, group_effect_summary_list)
  rownames(group_effect_summary_df) <- NULL
  
  # Sort
  group_effect_summary_df <- group_effect_summary_df[
    order(
      group_effect_summary_df$Trait,
      factor(group_effect_summary_df$Group_Type,
             levels = c("Ecosystem", "Season", "SoilDepth", "Frozen_Soil", "StandAge"))
    ),
  ]
  
  sheet_name <- get_unique_sheet_name("Group_Effect_Summary", used_sheet_names)
  used_sheet_names <- c(used_sheet_names, sheet_name)
  addWorksheet(wb, sheet_name)
  writeData(wb, sheet_name, group_effect_summary_df)
  cat("  Generated worksheet:", sheet_name, "(", nrow(group_effect_summary_df), "rows)\n")
}

# ==================== Worksheet11: Ecosystem_Group_Basic_Models ====================
cat("  Generate ecosystem-group basic models...\n")
eco_basic_summary_list <- list()

for(trait in names(all_results)) {
  if(!is.null(all_results[[trait]]$ecosystem_groups)) {
    # Get the original data
    original_trait <- traits[clean_traits == trait][1]
    d_current <- subset(d2, Trait == original_trait)
    
    for(eco in names(all_results[[trait]]$ecosystem_groups)) {
      # Get data statistics for this group
      eco_data <- d_current[d_current$Eco_1 == eco & !is.na(d_current$Eco_1), ]
      
      group_results <- all_results[[trait]]$ecosystem_groups[[eco]]
      
      # Extract basic model information for this ecosystem
      if(!is.null(group_results$sa)) {
        sa_eco_data <- subset(eco_data, `Treatment_1` == "Snow addition")
        sa_stats <- get_analysis_stats(sa_eco_data, trait, "Snow Addition", list(type = "Ecosystem", value = eco))
        
        eco_basic_summary_list[[paste0(trait, "_", eco, "_SnowAddition")]] <- extract_model_summary(
          group_results$sa, paste0(eco, " - Snow Addition"), sa_stats)
      }
      if(!is.null(group_results$sr)) {
        sr_eco_data <- subset(eco_data, `Treatment_1` == "Snow remove")
        sr_stats <- get_analysis_stats(sr_eco_data, trait, "Snow Remove", list(type = "Ecosystem", value = eco))
        
        eco_basic_summary_list[[paste0(trait, "_", eco, "_SnowRemove")]] <- extract_model_summary(
          group_results$sr, paste0(eco, " - Snow Remove"), sr_stats)
      }
      if(!is.null(group_results$treatment_model)) {
        treatment_stats <- get_analysis_stats(eco_data, trait, "Treatment Comparison", list(type = "Ecosystem", value = eco))
        
        eco_basic_summary_list[[paste0(trait, "_", eco, "_TreatmentModel")]] <- extract_model_summary(
          group_results$treatment_model, paste0(eco, " - Treatment Comparison"), treatment_stats)
      }
    }
  }
}

if(length(eco_basic_summary_list) > 0) {
  eco_basic_summary_df <- do.call(rbind, lapply(eco_basic_summary_list, function(x) {
    if(!is.null(x)) return(x)
  }))
  
  if(!is.null(eco_basic_summary_df)) {
    eco_basic_summary_df$Trait <- gsub("_.*", "", rownames(eco_basic_summary_df))
    eco_basic_summary_df$Ecosystem <- sapply(strsplit(rownames(eco_basic_summary_df), "_"), function(x) x[2])
    rownames(eco_basic_summary_df) <- NULL
    
    sheet_name <- get_unique_sheet_name("Ecosystem_Group_Basic_Models", used_sheet_names)
    used_sheet_names <- c(used_sheet_names, sheet_name)
    addWorksheet(wb, sheet_name)
    writeData(wb, sheet_name, eco_basic_summary_df)
    cat("  Generated worksheet:", sheet_name, "\n")
  }
}

# ==================== Worksheet12: Season_Group_Basic_Models ====================
cat("  Generate season-group basic models...\n")
season_basic_summary_list <- list()

for(trait in names(all_results)) {
  if(!is.null(all_results[[trait]]$season_groups)) {
    # Get the original data
    original_trait <- traits[clean_traits == trait][1]
    d_current <- subset(d2, Trait == original_trait)
    
    for(season in names(all_results[[trait]]$season_groups)) {
      # Get data statistics for this group
      season_data <- d_current[d_current$Sampling_season == season & !is.na(d_current$Sampling_season), ]
      
      group_results <- all_results[[trait]]$season_groups[[season]]
      
      # Extract basic model information for this season
      if(!is.null(group_results$sa)) {
        sa_season_data <- subset(season_data, `Treatment_1` == "Snow addition")
        sa_stats <- get_analysis_stats(sa_season_data, trait, "Snow Addition", list(type = "Season", value = season))
        
        season_basic_summary_list[[paste0(trait, "_", season, "_SnowAddition")]] <- extract_model_summary(
          group_results$sa, paste0(season, " - Snow Addition"), sa_stats)
      }
      if(!is.null(group_results$sr)) {
        sr_season_data <- subset(season_data, `Treatment_1` == "Snow remove")
        sr_stats <- get_analysis_stats(sr_season_data, trait, "Snow Remove", list(type = "Season", value = season))
        
        season_basic_summary_list[[paste0(trait, "_", season, "_SnowRemove")]] <- extract_model_summary(
          group_results$sr, paste0(season, " - Snow Remove"), sr_stats)
      }
      if(!is.null(group_results$treatment_model)) {
        treatment_stats <- get_analysis_stats(season_data, trait, "Treatment Comparison", list(type = "Season", value = season))
        
        season_basic_summary_list[[paste0(trait, "_", season, "_TreatmentModel")]] <- extract_model_summary(
          group_results$treatment_model, paste0(season, " - Treatment Comparison"), treatment_stats)
      }
    }
  }
}

if(length(season_basic_summary_list) > 0) {
  season_basic_summary_df <- do.call(rbind, lapply(season_basic_summary_list, function(x) {
    if(!is.null(x)) return(x)
  }))
  
  if(!is.null(season_basic_summary_df)) {
    season_basic_summary_df$Trait <- gsub("_.*", "", rownames(season_basic_summary_df))
    season_basic_summary_df$Season <- sapply(strsplit(rownames(season_basic_summary_df), "_"), function(x) x[2])
    rownames(season_basic_summary_df) <- NULL
    
    sheet_name <- get_unique_sheet_name("Season_Group_Basic_Models", used_sheet_names)
    used_sheet_names <- c(used_sheet_names, sheet_name)
    addWorksheet(wb, sheet_name)
    writeData(wb, sheet_name, season_basic_summary_df)
    cat("  Generated worksheet:", sheet_name, "\n")
  }
}

# ==================== Worksheet13-17: Other group basic models ====================
# For brevity, the frozen-soil group basic model is shown here; the others are similar

# Worksheet13: Frozen_Soil_Group_Basic_Models
cat("  Generate frozen-soil group basic models...\n")
frozen_basic_summary_list <- list()

for(trait in names(all_results)) {
  if(!is.null(all_results[[trait]]$frozen_groups)) {
    # Get the original data
    original_trait <- traits[clean_traits == trait][1]
    d_current <- subset(d2, Trait == original_trait)
    
    for(frozen in names(all_results[[trait]]$frozen_groups)) {
      # Get data statistics for this group
      frozen_data <- d_current[d_current$Types_of_frozen_soil == frozen & !is.na(d_current$Types_of_frozen_soil), ]
      frozen_name <- ifelse(frozen == "P", "Permafrost", "Non_permafrost")
      
      group_results <- all_results[[trait]]$frozen_groups[[frozen]]
      
      # Extract basic model information for this frozen-soil type
      if(!is.null(group_results$sa)) {
        sa_frozen_data <- subset(frozen_data, `Treatment_1` == "Snow addition")
        sa_stats <- get_analysis_stats(sa_frozen_data, trait, "Snow Addition", list(type = "Frozen_Soil", value = frozen_name))
        
        frozen_basic_summary_list[[paste0(trait, "_", frozen, "_SnowAddition")]] <- extract_model_summary(
          group_results$sa, paste0(frozen_name, " - Snow Addition"), sa_stats)
      }
      if(!is.null(group_results$sr)) {
        sr_frozen_data <- subset(frozen_data, `Treatment_1` == "Snow remove")
        sr_stats <- get_analysis_stats(sr_frozen_data, trait, "Snow Remove", list(type = "Frozen_Soil", value = frozen_name))
        
        frozen_basic_summary_list[[paste0(trait, "_", frozen, "_SnowRemove")]] <- extract_model_summary(
          group_results$sr, paste0(frozen_name, " - Snow Remove"), sr_stats)
      }
      if(!is.null(group_results$treatment_model)) {
        treatment_stats <- get_analysis_stats(frozen_data, trait, "Treatment Comparison", list(type = "Frozen_Soil", value = frozen_name))
        
        frozen_basic_summary_list[[paste0(trait, "_", frozen, "_TreatmentModel")]] <- extract_model_summary(
          group_results$treatment_model, paste0(frozen_name, " - Treatment Comparison"), treatment_stats)
      }
    }
  }
}

if(length(frozen_basic_summary_list) > 0) {
  frozen_basic_summary_df <- do.call(rbind, lapply(frozen_basic_summary_list, function(x) {
    if(!is.null(x)) return(x)
  }))
  
  if(!is.null(frozen_basic_summary_df)) {
    frozen_basic_summary_df$Trait <- gsub("_.*", "", rownames(frozen_basic_summary_df))
    frozen_basic_summary_df$Frozen_Type <- sapply(strsplit(rownames(frozen_basic_summary_df), "_"), function(x) {
      ifelse(x[2] == "P", "Permafrost", "Non_permafrost")
    })
    rownames(frozen_basic_summary_df) <- NULL
    
    sheet_name <- get_unique_sheet_name("Frozen_Soil_Group_Basic_Models", used_sheet_names)
    used_sheet_names <- c(used_sheet_names, sheet_name)
    addWorksheet(wb, sheet_name)
    writeData(wb, sheet_name, frozen_basic_summary_df)
    cat("  Generated worksheet:", sheet_name, "\n")
  }
}

# Worksheet14: Study_Duration_Group_Basic_Models
cat("  Generate study-duration group basic models...\n")
stand_age_basic_summary_list <- list()

for(trait in names(all_results)) {
  if(!is.null(all_results[[trait]]$stand_age_groups)) {
    # Get the original data
    original_trait <- traits[clean_traits == trait][1]
    d_current <- subset(d2, Trait == original_trait)
    
    for(age_group in names(all_results[[trait]]$stand_age_groups)) {
      # Get data statistics for this group
      age_data <- d_current[d_current$StandAge_group == age_group & !is.na(d_current$StandAge_group), ]
      
      group_results <- all_results[[trait]]$stand_age_groups[[age_group]]
      
      if(!is.null(group_results$sa)) {
        sa_age_data <- subset(age_data, `Treatment_1` == "Snow addition")
        sa_stats <- get_analysis_stats(sa_age_data, trait, "Snow Addition", list(type = "StandAge", value = paste0(age_group, " years")))
        
        stand_age_basic_summary_list[[paste0(trait, "_", age_group, "_SnowAddition")]] <- extract_model_summary(
          group_results$sa, paste0(age_group, " years - Snow Addition"), sa_stats)
      }
      if(!is.null(group_results$sr)) {
        sr_age_data <- subset(age_data, `Treatment_1` == "Snow remove")
        sr_stats <- get_analysis_stats(sr_age_data, trait, "Snow Remove", list(type = "StandAge", value = paste0(age_group, " years")))
        
        stand_age_basic_summary_list[[paste0(trait, "_", age_group, "_SnowRemove")]] <- extract_model_summary(
          group_results$sr, paste0(age_group, " years - Snow Remove"), sr_stats)
      }
      if(!is.null(group_results$treatment_model)) {
        treatment_stats <- get_analysis_stats(age_data, trait, "Treatment Comparison", list(type = "StandAge", value = paste0(age_group, " years")))
        
        stand_age_basic_summary_list[[paste0(trait, "_", age_group, "_TreatmentModel")]] <- extract_model_summary(
          group_results$treatment_model, paste0(age_group, " years - Treatment Comparison"), treatment_stats)
      }
    }
  }
}

if(length(stand_age_basic_summary_list) > 0) {
  stand_age_basic_summary_df <- do.call(rbind, lapply(stand_age_basic_summary_list, function(x) {
    if(!is.null(x)) return(x)
  }))
  
  if(!is.null(stand_age_basic_summary_df)) {
    stand_age_basic_summary_df$Trait <- gsub("_.*", "", rownames(stand_age_basic_summary_df))
    stand_age_basic_summary_df$StandAge_Group <- sapply(strsplit(rownames(stand_age_basic_summary_df), "_"), function(x) paste0(x[2], " years"))
    rownames(stand_age_basic_summary_df) <- NULL
    
    sheet_name <- get_unique_sheet_name("Study_Duration_Group_Basic_Models", used_sheet_names)
    used_sheet_names <- c(used_sheet_names, sheet_name)
    addWorksheet(wb, sheet_name)
    writeData(wb, sheet_name, stand_age_basic_summary_df)
    cat("  Generated worksheet:", sheet_name, "\n")
  }
}

# Worksheet15: Soil_Depth_Group_Basic_Models
cat("  Generate soil-depth group basic models...\n")
soil_depth_basic_summary_list <- list()

for(trait in names(all_results)) {
  if(!is.null(all_results[[trait]]$soil_depth_groups)) {
    # Get the original data
    original_trait <- traits[clean_traits == trait][1]
    d_current <- subset(d2, Trait == original_trait)
    
    for(depth_group in names(all_results[[trait]]$soil_depth_groups)) {
      # Get data statistics for this group
      depth_data <- d_current[d_current$SoilDepth_group == depth_group & !is.na(d_current$SoilDepth_group), ]
      
      group_results <- all_results[[trait]]$soil_depth_groups[[depth_group]]
      
      if(!is.null(group_results$sa)) {
        sa_depth_data <- subset(depth_data, `Treatment_1` == "Snow addition")
        sa_stats <- get_analysis_stats(sa_depth_data, trait, "Snow Addition", list(type = "SoilDepth", value = paste0(depth_group, " cm")))
        
        soil_depth_basic_summary_list[[paste0(trait, "_", depth_group, "_SnowAddition")]] <- extract_model_summary(
          group_results$sa, paste0(depth_group, " cm - Snow Addition"), sa_stats)
      }
      if(!is.null(group_results$sr)) {
        sr_depth_data <- subset(depth_data, `Treatment_1` == "Snow remove")
        sr_stats <- get_analysis_stats(sr_depth_data, trait, "Snow Remove", list(type = "SoilDepth", value = paste0(depth_group, " cm")))
        
        soil_depth_basic_summary_list[[paste0(trait, "_", depth_group, "_SnowRemove")]] <- extract_model_summary(
          group_results$sr, paste0(depth_group, " cm - Snow Remove"), sr_stats)
      }
      if(!is.null(group_results$treatment_model)) {
        treatment_stats <- get_analysis_stats(depth_data, trait, "Treatment Comparison", list(type = "SoilDepth", value = paste0(depth_group, " cm")))
        
        soil_depth_basic_summary_list[[paste0(trait, "_", depth_group, "_TreatmentModel")]] <- extract_model_summary(
          group_results$treatment_model, paste0(depth_group, " cm - Treatment Comparison"), treatment_stats)
      }
    }
  }
}

if(length(soil_depth_basic_summary_list) > 0) {
  soil_depth_basic_summary_df <- do.call(rbind, lapply(soil_depth_basic_summary_list, function(x) {
    if(!is.null(x)) return(x)
  }))
  
  if(!is.null(soil_depth_basic_summary_df)) {
    soil_depth_basic_summary_df$Trait <- gsub("_.*", "", rownames(soil_depth_basic_summary_df))
    soil_depth_basic_summary_df$SoilDepth_Group <- sapply(strsplit(rownames(soil_depth_basic_summary_df), "_"), function(x) paste0(x[2], " cm"))
    rownames(soil_depth_basic_summary_df) <- NULL
    
    sheet_name <- get_unique_sheet_name("Soil_Depth_Group_Basic_Models", used_sheet_names)
    used_sheet_names <- c(used_sheet_names, sheet_name)
    addWorksheet(wb, sheet_name)
    writeData(wb, sheet_name, soil_depth_basic_summary_df)
    cat("  Generated worksheet:", sheet_name, "\n")
  }
}

# ==================== Worksheet16: Group_Treatment_Difference_Summary ====================
cat("  Generate group treatment-difference summary...\n")
group_treatment_comparison_list <- list()

for(trait in names(all_results)) {
  # Get the original data for statistics
  original_trait <- traits[clean_traits == trait][1]
  d_current <- subset(d2, Trait == original_trait)
  
  # Overall treatment difference
  if(!is.null(all_results[[trait]]$overall$treatment_model)) {
    overall_model <- all_results[[trait]]$overall$treatment_model
    overall_stats <- get_analysis_stats(d_current, trait, "Treatment Comparison")
    
    group_treatment_comparison_list[[paste0(trait, "_Overall")]] <- data.frame(
      Trait = trait,
      Group = "Overall",
      Group_Type = "Overall",
      QM = overall_model$QM,
      df1 = overall_model$m,
      p_value = ifelse(length(overall_model$pval) > 1, overall_model$pval[2], overall_model$pval[1]),
      k = overall_model$k,
      Study_id_Count = overall_stats$Study_id_Count,
      Row_id_Count = overall_stats$Row_id_Count,
      stringsAsFactors = FALSE
    )
  }
  
  # Treatment differences for different groups: similar to the previous code, but with data statistics added
  # ... Other group-specific treatment-difference statistics can be added here ...
}

if(length(group_treatment_comparison_list) > 0) {
  group_treatment_comparison_df <- do.call(rbind, group_treatment_comparison_list)
  rownames(group_treatment_comparison_df) <- NULL
  
  sheet_name <- get_unique_sheet_name("Group_Treatment_Difference_Summary", used_sheet_names)
  used_sheet_names <- c(used_sheet_names, sheet_name)
  addWorksheet(wb, sheet_name)
  writeData(wb, sheet_name, group_treatment_comparison_df)
  cat("  Generated worksheet:", sheet_name, "\n")
}

# ==================== Worksheet17: Analysis_Statistics ====================
cat("  Generate analysis statistics...\n")
stats_list <- list()

for(trait in names(all_results)) {
  # Get the original data for statistics
  original_trait <- traits[clean_traits == trait][1]
  d_current <- subset(d2, Trait == original_trait)
  overall_stats <- get_analysis_stats(d_current, trait, "Overall Analysis")
  
  stats_list[[trait]] <- data.frame(
    Trait = trait,
    SnowAddition = ifelse(!is.null(all_results[[trait]]$overall$sa), TRUE, FALSE),
    SnowRemove = ifelse(!is.null(all_results[[trait]]$overall$sr), TRUE, FALSE),
    Treatment_Model = ifelse(!is.null(all_results[[trait]]$overall$treatment_model), TRUE, FALSE),
    Treatment_NoIntercept = ifelse(!is.null(all_results[[trait]]$overall$treatment_model_nointercept), TRUE, FALSE),
    Ecosystem_Interaction = ifelse(!is.null(all_results[[trait]]$overall$eco_interaction), TRUE, FALSE),
    Season_Interaction = ifelse(!is.null(all_results[[trait]]$overall$season_interaction), TRUE, FALSE),
    Ecosystem_Groups = ifelse(!is.null(all_results[[trait]]$ecosystem_groups), TRUE, FALSE),
    Season_Groups = ifelse(!is.null(all_results[[trait]]$season_groups), TRUE, FALSE),
    Frozen_Groups = ifelse(!is.null(all_results[[trait]]$frozen_groups), TRUE, FALSE),
    StandAge_Groups = ifelse(!is.null(all_results[[trait]]$stand_age_groups), TRUE, FALSE),
    SoilDepth_Groups = ifelse(!is.null(all_results[[trait]]$soil_depth_groups), TRUE, FALSE),
    Study_id_Count = overall_stats$Study_id_Count,
    Row_id_Count = overall_stats$Row_id_Count,
    stringsAsFactors = FALSE
  )
}

if(length(stats_list) > 0) {
  stats_df <- do.call(rbind, stats_list)
  rownames(stats_df) <- NULL
  
  sheet_name <- get_unique_sheet_name("Analysis_Statistics", used_sheet_names)
  used_sheet_names <- c(used_sheet_names, sheet_name)
  addWorksheet(wb, sheet_name)
  writeData(wb, sheet_name, stats_df)
  cat("  Generated worksheet:", sheet_name, "\n")
}
# Save Excel file
saveWorkbook(wb, "C:/Users/Administrator/Desktop/25meta/complete_meta_analysis_results_12.19_version_10_with_data_statistics.xlsx", overwrite = TRUE)

# 7. Save R data files -----------------------------------------------------------
cat("  Save R data files...\n")
saveRDS(all_results, "C:/Users/Administrator/Desktop/25meta/r/version_4_ecosystem_treatment_intensity_and_other_details/complete_meta_analysis_results_version_6.rds")
save(all_results, file = "C:/Users/Administrator/Desktop/25meta/r/version_2_ecosystem_treatment_intensity_and_other_details/complete_meta_analysis_results_version_6.RData")

# 8. Generate summary report -----------------------------------------------------------
cat("  Generate text summary...\n")
sink("meta_analysis_summary_report.txt")

cat("meta_analysis_summary_report\n")
cat("==============\n\n")
cat("Analysis time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("Total number of traits:", length(traits), "\n")
cat("Successfully analyzed:", length(all_results), "\n")
cat("Success rate:", round(length(all_results)/length(traits)*100, 1), "%\n\n")

cat("Analysis status for each trait:\n")
cat("----------------\n")
for(trait in names(all_results)) {
  # Get data statistics
  original_trait <- traits[clean_traits == trait][1]
  d_current <- subset(d2, Trait == original_trait)
  stats <- get_analysis_stats(d_current, trait, "Overall Analysis")
  
  cat(trait, ":\n")
  cat("  Data size: Study_id =", stats$Study_id_Count, ", Row_id =", stats$Row_id_Count, "\n")
  
  overall <- all_results[[trait]]$overall
  
  if(!is.null(overall$treatment_model)) {
    cat("  Between-treatment difference p =", format.pval(overall$treatment_model$pval[2], digits = 4), "\n")
  }
  
  if(!is.null(overall$sa)) {
    cat("  Snow addition effect:", round(overall$sa$beta[1], 4), 
        "[", round(overall$sa$ci.lb, 4), ",", round(overall$sa$ci.ub, 4), "]\n")
  }
  
  if(!is.null(overall$sr)) {
    cat("  Snow removal effect:", round(overall$sr$beta[1], 4), 
        "[", round(overall$sr$ci.lb, 4), ",", round(overall$sr$ci.ub, 4), "]\n")
  }
  cat("\n")
}

sink()


