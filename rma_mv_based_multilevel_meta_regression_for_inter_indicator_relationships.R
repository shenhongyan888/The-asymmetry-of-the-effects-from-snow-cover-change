# =============================================================================
# ST Regression Analysis - Final Revised Version
# =============================================================================

# 1. Load packages -----------------------------------------------------------
cat("Step 1: Loading required packages...\n")
required_packages <- c("metafor", "readxl", "multcomp", "openxlsx", "dplyr", "purrr", "ggplot2")
new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages)

library(metafor)
library(readxl)
library(multcomp)
library(openxlsx)
library(dplyr)
library(purrr)
library(ggplot2)
cat("All packages were loaded successfully.\n\n")

# 2. Import and preprocess data ----------------------------------------------
cat("Step 2: Importing and preprocessing data...\n")
d2 <- read_excel("E:/\u535a\u58eb\u671f\u95f4\u6587\u6863/25meta/\u56de\u5f52\u5206\u6790/ST\u56de\u5f52\u6570\u636e\u6548\u5e94\u503c.xlsx")

cat("\n=== Data type conversion ===\n")
if("yi" %in% names(d2)) {
  d2$yi <- as.numeric(as.character(d2$yi))
  cat("Non-missing yi values:", sum(!is.na(d2$yi)), "/", nrow(d2), "\n")
} else {
  stop("The yi column does not exist.")
}
if("V" %in% names(d2)) {
  d2$V <- as.numeric(as.character(d2$V))
  cat("Non-missing V values:", sum(!is.na(d2$V)), "/", nrow(d2), "\n")
} else {
  stop("The V column does not exist.")
}
if("vi" %in% names(d2)) {
  d2$vi <- as.numeric(as.character(d2$vi))
} else {
  cat("The vi column does not exist; the V column will be used instead.\n")
  d2$vi <- d2$V
}
if("Treatment_1" %in% names(d2)) {
  cat("\nDistribution of Treatment_1:\n")
  print(table(d2$Treatment_1, useNA = "always"))
} else {
  stop("The Treatment_1 column does not exist.")
}

# Check random-effects columns
required_random_cols <- c("Study_id", "Row_id", "Sampling_season", "StandAge_year", "Soil_depth")
cat("\n=== Checking random-effects columns ===\n")
for(col in required_random_cols) {
  if(col %in% names(d2)) {
    cat(sprintf("%s exists.\n", col))
  } else {
    cat(sprintf("%s does not exist; a default value will be created.\n", col))
    d2[[col]] <- "Default"
  }
}

# 3. Define the variable list and column mapping -----------------------------
cat("\nStep 3: Defining variable mapping...\n")
variables_list <- unique(c(
  "Soil temperature",
  "Soil water content",
  "Snowmelt date",
  "Number of frozen days",
  "Frost depth",
  "Number of freeze-thaw cycles",
  "Soil water content",
  "pH",
  "Cation exchange capacity",
  "CO2 flux",
  "N2O flux",
  "CH4 flux",
  "Total carbon",
  "Dissolved organic carbon",
  "Soil organic carbon",
  "Total nitrogen",
  "Dissolved organic nitrogen",
  "Nitrate nitrogen",
  "Ammonium nitrogen",
  "Inorganic nitrogen",
  "Nitrogen limitation",
  "Nitrogen availability",
  "Carbon to nitrogen ratio",
  "Nitrogen to phosphorus ratio",
  "Ammonification",
  "Mineralization rate",
  "Nitrification",
  "Total phosphorus",
  "Available phosphorus",
  "Phosphate",
  "Soil Ca",
  "Soil Mg",
  "Soil K",
  "Leaf carbon",
  "Leaf nitrogen",
  "Leaf phosphorus",
  "Total PLFA",
  "Fungal PLFA",
  "Bacterial PLFA",
  "Microbial biomass carbon",
  "Microbial biomass nitrogen",
  "Microbial biomass phosphorus",
  "MBC/MBN",
  "MBN/MBP",
  "Bacterial Chao1 index",
  "Bacterial Shannon",
  "Fungal Shannon",
  "Bacterial Simpson",
  "β-glucosidase",
  "β-n-acetylglucosaminidase",
  "Nitrate reductase",
  "Nitrite reductase",
  "Urease",
  "Invertase",
  "β-xylosidase",
  "Peroxidase",
  "Polyphenol oxidase",
  "Phosphate",
  "Acid phosphatase",
  "Aboveground biomass",
  "Belowground biomass",
  "Plant cover",
  "Plant branch length",
  "Plant root length",
  "Plant first budding date",
  "Gross primary production",
  "Net ecosystem exchange",
  "Litter mass",
  "Litter decomposition rate",
  "Carbon remaining in litter",
  "Nitrogen remaining in litter",
  "Phosphorus remaining in litter"
))

# Each display name is linked to one or more possible source-column names.
# The first matching column present in the input dataset will be used.
column_mapping <- list(
  "Soil temperature" = c("ST", "Soil temperature"),
  "Soil water content" = c("SWC", "Soil water content"),
  "Snowmelt date" = c("Thaw_date", "Snowmelt_date", "Snowmelt date"),
  "Number of frozen days" = c("Frost_days", "Number_of_frozen_days", "Number of frozen days"),
  "Frost depth" = c("Frost_depth", "Frost depth"),
  "Number of freeze-thaw cycles" = c(
    "Freeze_thaw_cycles",
    "Number_of_freeze_thaw_cycles",
    "Number of freeze-thaw cycles"
  ),
  "pH" = c("pH"),
  "Cation exchange capacity" = c("CEC", "Cation exchange capacity"),
  "CO2 flux" = c("CO2", "CO2_flux", "CO2 flux"),
  "N2O flux" = c("N2O", "N2O_flux", "N2O flux"),
  "CH4 flux" = c("CH4", "CH4_flux", "CH4 flux"),
  "Total carbon" = c("TC", "Total carbon"),
  "Dissolved organic carbon" = c("DOC", "Dissolved organic carbon"),
  "Soil organic carbon" = c("SOC", "Soil organic carbon"),
  "Total nitrogen" = c("TN", "Total nitrogen"),
  "Dissolved organic nitrogen" = c("DON", "Dissolved organic nitrogen"),
  "Nitrate nitrogen" = c("NO3-", "Nitrate nitrogen"),
  "Ammonium nitrogen" = c("NH4+", "Ammonium nitrogen"),
  "Inorganic nitrogen" = c("Inorganic nitrogen", "DIN"),
  "Nitrogen limitation" = c("Nitrogen limitation", "N_limitation", "N limitation"),
  "Nitrogen availability" = c("Nitrogen availability", "N_availability", "N availability"),
  "Carbon to nitrogen ratio" = c("C/N", "Carbon_to_nitrogen_ratio", "Carbon to nitrogen ratio"),
  "Nitrogen to phosphorus ratio" = c("N/P", "Nitrogen_to_phosphorus_ratio", "Nitrogen to phosphorus ratio"),
  "Ammonification" = c("Ammonification"),
  "Mineralization rate" = c("Mineralization_rate", "Mineralization rate"),
  "Nitrification" = c("Nitrification"),
  "Total phosphorus" = c("TP", "Total phosphorus"),
  "Available phosphorus" = c("AP", "Available phosphorus"),
  "Phosphate" = c("PO43-", "Phosphate"),
  "Soil Ca" = c("Soil_Ca", "Soil Ca"),
  "Soil Mg" = c("Soil_Mg", "Soil Mg"),
  "Soil K" = c("Soil_K", "Soil K"),
  "Leaf carbon" = c("Leaf_C", "Leaf carbon"),
  "Leaf nitrogen" = c("Leaf_N", "Leaf nitrogen"),
  "Leaf phosphorus" = c("Leaf_P", "Leaf phosphorus"),
  "Total PLFA" = c("Totel_PLFA", "Total_PLFA", "Total PLFA"),
  "Fungal PLFA" = c("F_PLFA", "Fungal_PLFA", "Fungal PLFA"),
  "Bacterial PLFA" = c("B_PLFA", "Bacterial_PLFA", "Bacterial PLFA"),
  "Microbial biomass carbon" = c("MBC", "Microbial biomass carbon"),
  "Microbial biomass nitrogen" = c("MBN", "Microbial biomass nitrogen"),
  "Microbial biomass phosphorus" = c("MBP", "Microbial biomass phosphorus"),
  "MBC/MBN" = c("MBC/MBN"),
  "MBN/MBP" = c("MBN/MBP"),
  "Bacterial Chao1 index" = c("B_Chao1", "Bacterial_Chao1_index", "Bacterial Chao1 index"),
  "Bacterial Shannon" = c("B_Shannon", "Bacterial_Shannon", "Bacterial Shannon"),
  "Fungal Shannon" = c("F_Shannon", "Fungal_Shannon", "Fungal Shannon"),
  "Bacterial Simpson" = c("B_Simpson's", "B_Simpson", "Bacterial_Simpson", "Bacterial Simpson"),
  "β-glucosidase" = c("βG", "BG", "Beta_glucosidase", "β-glucosidase"),
  "β-n-acetylglucosaminidase" = c(
    "NAG",
    "Beta_n_acetylglucosaminidase",
    "β-n-acetylglucosaminidase"
  ),
  "Nitrate reductase" = c("NaRS", "Nitrate_reductase", "Nitrate reductase"),
  "Nitrite reductase" = c("NiRS", "Nitrite_reductase", "Nitrite reductase"),
  "Urease" = c("Urease"),
  "Invertase" = c("Invertase"),
  "β-xylosidase" = c("βX", "BX", "Beta_xylosidase", "β-xylosidase"),
  "Peroxidase" = c("Peroxidase", "POD"),
  "Polyphenol oxidase" = c("Polyphenol_oxidase", "Polyphenol oxidase", "PPO"),
  "Acid phosphatase" = c("Acid_phosphatase", "Acid phosphatase", "ACP"),
  "Aboveground biomass" = c("AGB", "Aboveground biomass"),
  "Belowground biomass" = c("BGB", "Belowground biomass"),
  "Plant cover" = c("Plant_cover", "Plant cover"),
  "Plant branch length" = c("Plant_branch_length", "Plant branch length"),
  "Plant root length" = c("Plant_root_length", "Plant root length"),
  "Plant first budding date" = c(
    "Plant_the_first_budding_deat",
    "Plant_the_first_budding_date",
    "Plant_first_budding_date",
    "Plant first budding date"
  ),
  "Gross primary production" = c("GPP", "Gross_primary_production", "Gross primary production"),
  "Net ecosystem exchange" = c("NEE", "Net_ecosystem_exchange", "Net ecosystem exchange"),
  "Litter mass" = c("Litter_mass", "Litter mass"),
  "Litter decomposition rate" = c("Litter_decomposition_rate", "Litter decomposition rate"),
  "Carbon remaining in litter" = c(
    "Litter_C",
    "Carbon_remaining_in_litter",
    "Carbon remaining in litter"
  ),
  "Nitrogen remaining in litter" = c(
    "Litter_N",
    "Nitrogen_remaining_in_litter",
    "Nitrogen remaining in litter"
  ),
  "Phosphorus remaining in litter" = c(
    "Litter_P",
    "Phosphorus_remaining_in_litter",
    "Phosphorus remaining in litter"
  )
)

# 4. Create output directories -----------------------------------------------
main_output_path <- "E:/\u535a\u58eb\u671f\u95f4\u6587\u6863/25meta/\u56de\u5f52\u5206\u6790/2026.6.28/ST"
if (!dir.exists(main_output_path)) dir.create(main_output_path, recursive = TRUE)
combo_path <- file.path(main_output_path, "combined_plots")
if (!dir.exists(combo_path)) dir.create(combo_path, recursive = TRUE)

all_results_summary <- data.frame()
regression_coefficients <- data.frame()

# ============================
# Helper functions
# ============================
get_column_name <- function(variable_name) {
  if(variable_name %in% names(column_mapping)) {
    candidate_names <- unname(column_mapping[[variable_name]])
    matched_names <- candidate_names[candidate_names %in% names(d2)]
    if(length(matched_names) > 0) return(matched_names[1])
  }
  return(NULL)
}

fit_simple_model <- function(data, col_name) {
  if(nrow(data) < 3) return(NULL)
  data$yi <- as.numeric(as.character(data$yi))
  data$V <- as.numeric(as.character(data$V))
  data$vi <- as.numeric(as.character(data$vi))
  data[[col_name]] <- as.numeric(as.character(data[[col_name]]))
  data$Study_id <- as.factor(data$Study_id)
  data$Row_id <- as.factor(data$Row_id)
  data$Sampling_season <- as.factor(data$Sampling_season)
  data$StandAge_year <- as.factor(data$StandAge_year)
  data$Soil_depth <- as.factor(data$Soil_depth)
  data <- data[!is.na(data$yi) & !is.na(data$V) & !is.na(data[[col_name]]), ]
  if(nrow(data) < 3) return(NULL)
  tryCatch({
    model <- rma.mv(
      yi = data$yi, V = data$V, mods = ~ data[[col_name]],
      random = list(~ 1 | Study_id/Row_id, ~ 1 | Sampling_season,
                    ~ 1 | StandAge_year, ~ 1 | Soil_depth),
      data = data, method = "REML", test = "t", control = list(maxiter = 1000)
    )
    return(model)
  }, error = function(e) {
    cat("  The full model failed; attempting a simplified model...\n")
    tryCatch({
      model_simple <- rma.mv(
        yi = data$yi, V = data$V, mods = ~ data[[col_name]],
        random = list(~ 1 | Study_id/Row_id),
        data = data, method = "REML", test = "t", control = list(maxiter = 1000)
      )
      return(model_simple)
    }, error = function(e2) {
      cat("  The simplified model also failed:", e2$message, "\n")
      return(NULL)
    })
  })
}

get_ticks <- function(min_val, max_val) {
  if(min_val == max_val) { min_val <- min_val - 0.5; max_val <- max_val + 0.5 }
  min_int <- floor(min_val); max_int <- ceiling(max_val)
  if(max_int - min_int < 3) max_int <- min_int + 3
  step <- (max_int - min_int) / 3
  if(step != round(step)) max_int <- min_int + 3 * ceiling(step)
  ticks <- seq(min_int, max_int, length.out = 4)
  return(list(ticks = ticks, range = c(min_int, max_int)))
}

# Regression y-axis: retain four equally spaced ticks and include all observations and confidence intervals
get_regression_y_info <- function(values, margin_ratio = 0.05) {
  values <- values[is.finite(values)]
  if(length(values) == 0) return(get_ticks(-1, 1))
  min_val <- min(values)
  max_val <- max(values)
  value_range <- max_val - min_val
  if(!is.finite(value_range) || value_range <= 0) {
    value_range <- max(abs(c(min_val, max_val)), 1)
  }
  margin <- margin_ratio * value_range
  return(get_ticks(min_val - margin, max_val + margin))
}

# Bar-chart y-axis: use three equally spaced ticks with equal numerical intervals and maximize vertical bar coverage
get_bar_axis <- function(beta_vals, margin_ratio = 0.08) {
  beta_vals <- beta_vals[is.finite(beta_vals)]
  if(length(beta_vals) == 0) {
    return(list(ticks = c(-0.1, 0.0, 0.1), range = c(-0.1, 0.1)))
  }
  
  beta_min <- min(beta_vals)
  beta_max <- max(beta_vals)
  
  # All values are nonnegative: start the lower axis at 0.0 with no lower margin
  if(beta_min >= 0) {
    target_upper <- beta_max * (1 + margin_ratio)
    upper <- ceiling(target_upper / 0.2) * 0.2
    if(upper <= beta_max + 1e-10) upper <- upper + 0.2
    upper <- max(upper, 0.2)
    ticks <- c(0, upper / 2, upper)
    return(list(ticks = round(ticks, 10), range = c(0, upper)))
  }
  
  # All values are nonpositive: end the upper axis at 0.0 with no upper margin
  if(beta_max <= 0) {
    target_lower_abs <- abs(beta_min) * (1 + margin_ratio)
    lower_abs <- ceiling(target_lower_abs / 0.2) * 0.2
    if(lower_abs <= abs(beta_min) + 1e-10) lower_abs <- lower_abs + 0.2
    lower_abs <- max(lower_abs, 0.2)
    ticks <- c(-lower_abs, -lower_abs / 2, 0)
    return(list(ticks = round(ticks, 10), range = c(-lower_abs, 0)))
  }
  
  # Values have opposite signs: use 0.0 as the middle tick with no additional end margins
  max_abs <- max(abs(beta_vals))
  target_limit <- max_abs * (1 + margin_ratio)
  limit <- ceiling(target_limit / 0.1) * 0.1
  if(limit <= max_abs + 1e-10) limit <- limit + 0.1
  limit <- max(limit, 0.1)
  ticks <- c(-limit, 0, limit)
  return(list(ticks = round(ticks, 10), range = c(-limit, limit)))
}

# ============================
# Main loop
# ============================
cat("\n", strrep("=", 80), "\n")
cat("Starting the analysis of all variables using rma.mv models...\n")
cat("Total number of variables:", length(variables_list), "\n")
cat(strrep("=", 80), "\n\n")

success_count <- 0

for (i in seq_along(variables_list)) {
  variable_name <- variables_list[i]
  cat("\n", strrep("-", 60), "\n")
  cat("Analyzing variable:", i, "/", length(variables_list), "-", variable_name, "\n")
  cat(strrep("-", 60), "\n")
  
  col_name <- get_column_name(variable_name)
  if(is.null(col_name)) { cat("The corresponding column could not be found; skipping this variable.\n"); next }
  cat("Column used:", col_name, "\n")
  
  safe_name <- gsub("[^[:alnum:]_]", "_", variable_name)
  variable_output_path <- file.path(main_output_path, safe_name)
  if (!dir.exists(variable_output_path)) dir.create(variable_output_path, recursive = TRUE)
  
  setwd(variable_output_path)
  log_file <- file.path(variable_output_path, paste0("analysis_log_", safe_name, ".txt"))
  sink(log_file, split = TRUE)
  
  cat("Variable:", variable_name, "\n")
  cat("Corresponding column:", col_name, "\n")
  cat("Analysis time:", Sys.time(), "\n\n")
  
  # Data preprocessing
  cat("=== Data preprocessing ===\n")
  d2[[col_name]] <- as.numeric(as.character(d2[[col_name]]))
  n_valid <- sum(!is.na(d2[[col_name]]) & !is.na(d2$yi) & !is.na(d2$V) & !is.na(d2$vi))
  cat("Total valid observations:", n_valid, "\n")
  if (n_valid < 3) { cat("Insufficient valid observations.\n"); sink(); next }
  
  snow_add <- subset(d2, Treatment_1 == "Snow addition" & 
                       !is.na(d2[[col_name]]) & !is.na(yi) & !is.na(V) & !is.na(vi))
  snow_remove <- subset(d2, Treatment_1 == "Snow remove" & 
                          !is.na(d2[[col_name]]) & !is.na(yi) & !is.na(V) & !is.na(vi))
  cat("Snow addition:", nrow(snow_add), "observations\n")
  cat("Snow removal:", nrow(snow_remove), "observations\n")
  
  # Fit models
  cat("\n=== Model fitting using rma.mv ===\n")
  model_add <- NULL; model_remove <- NULL
  if(nrow(snow_add) >= 3) {
    cat("Fitting the snow-addition model...\n")
    model_add <- fit_simple_model(snow_add, col_name)
    if(!is.null(model_add)) saveRDS(model_add, file.path(variable_output_path, paste0("model_snow_add_", safe_name, ".rds")))
  } else { cat("Insufficient snow-addition data (n < 3).\n") }
  if(nrow(snow_remove) >= 3) {
    cat("\nFitting the snow-removal model...\n")
    model_remove <- fit_simple_model(snow_remove, col_name)
    if(!is.null(model_remove)) saveRDS(model_remove, file.path(variable_output_path, paste0("model_snow_remove_", safe_name, ".rds")))
  } else { cat("Insufficient snow-removal data (n < 3).\n") }
  
  # Summarize results
  cat("\n=== Results summary ===\n")
  variable_results <- data.frame()
  if(!is.null(model_add)) {
    results_add <- data.frame(Variable = variable_name, Treatment = "Snow addition",
                              Beta = model_add$b[2], SE = model_add$se[2],
                              CI_lower = model_add$ci.lb[2], CI_upper = model_add$ci.ub[2],
                              p_value = model_add$pval[2], QM = model_add$QM, QM_p = model_add$QMp,
                              n_obs = nrow(snow_add), stringsAsFactors = FALSE)
    variable_results <- rbind(variable_results, results_add)
    success_count <- success_count + 1
    regression_coefficients <- rbind(regression_coefficients,
                                     data.frame(Variable = variable_name, Treatment = "Increased snowpack thickness",
                                                Beta = round(model_add$b[2], 4),
                                                SE = round(model_add$se[2], 4),
                                                CI_lower = model_add$ci.lb[2],
                                                CI_upper = model_add$ci.ub[2],
                                                p_value = round(model_add$pval[2], 4),
                                                n = nrow(snow_add), stringsAsFactors = FALSE))
  }
  if(!is.null(model_remove)) {
    results_remove <- data.frame(Variable = variable_name, Treatment = "Snow remove",
                                 Beta = model_remove$b[2], SE = model_remove$se[2],
                                 CI_lower = model_remove$ci.lb[2], CI_upper = model_remove$ci.ub[2],
                                 p_value = model_remove$pval[2], QM = model_remove$QM, QM_p = model_remove$QMp,
                                 n_obs = nrow(snow_remove), stringsAsFactors = FALSE)
    variable_results <- rbind(variable_results, results_remove)
    success_count <- success_count + 1
    regression_coefficients <- rbind(regression_coefficients,
                                     data.frame(Variable = variable_name, Treatment = "Decreased snowpack thickness",
                                                Beta = round(model_remove$b[2], 4),
                                                SE = round(model_remove$se[2], 4),
                                                CI_lower = model_remove$ci.lb[2],
                                                CI_upper = model_remove$ci.ub[2],
                                                p_value = round(model_remove$pval[2], 4),
                                                n = nrow(snow_remove), stringsAsFactors = FALSE))
  }
  if(nrow(variable_results) > 0) {
    all_results_summary <- rbind(all_results_summary, variable_results)
    write.csv(variable_results, file.path(variable_output_path, paste0("results_", safe_name, ".csv")), row.names = FALSE)
    cat("Results were saved successfully.\n")
  } else { cat("No valid model results were available.\n") }
  
  # =====================================================================
  # Plotting section
  # =====================================================================
  if(!is.null(model_add) || !is.null(model_remove)) {
    cat("\n=== Generating figures ===\n")
    if(sink.number() > 0) sink()
    
    # ---- 1. Separate plots ----------------------------------------------------
    try({
      img_width <- 1890; img_height <- 1890
      par_original <- par(no.readonly = TRUE)
      
      # Increased-snowpack-thickness plot
      if(!is.null(model_add) && nrow(snow_add) >= 3) {
        png_file_add <- file.path(variable_output_path, paste0("increased_snowpack_thickness_", safe_name, ".png"))
        png(png_file_add, width = img_width, height = img_height, units = "px", res = 600)
        par(family = "Arial", oma = c(1.5, 1.5, 0, 0), mar = c(1.5, 2.5, 1.0, 0.5),
            mgp = c(1.5, 0.2, 0), tcl = -0.02, cex.axis = 1.0, las = 1, lwd = 1.2,
            font = 1, pty = "s")
        x_data <- snow_add[[col_name]]; y_data <- snow_add$yi; vi_data <- snow_add$vi
        if(length(vi_data) > 0 && all(!is.na(vi_data)) && all(vi_data > 0)) {
          inv_vi <- 1/vi_data
          if(length(inv_vi) > 1 && var(inv_vi) > 0) {
            point_sizes <- 0.5 + 2.0 * (inv_vi - min(inv_vi)) / (max(inv_vi) - min(inv_vi))
          } else { point_sizes <- rep(1.5, length(x_data)) }
        } else { point_sizes <- rep(1.5, length(x_data)) }
        x_info <- get_ticks(min(x_data), max(x_data)); x_ticks <- x_info$ticks; xlim <- x_info$range
        pred_range <- seq(xlim[1], xlim[2], length = 100)
        preds <- predict(model_add, newmods = pred_range)
        y_info <- get_regression_y_info(c(y_data, preds$ci.lb, preds$ci.ub), margin_ratio = 0.05)
        y_ticks <- y_info$ticks
        ylim <- y_info$range
        plot(NA, xlim = xlim, ylim = ylim, xlab = "", ylab = "", axes = FALSE, lwd = 1.2, yaxs = "i")
        polygon(c(pred_range, rev(pred_range)), c(preds$ci.lb, rev(preds$ci.ub)),
                col = rgb(254,241,231,maxColorValue=255), border = NA)
        points(x_data, y_data, pch = 21, bg = rgb(255,192,127,maxColorValue=255),
               col = gray(0.2), cex = point_sizes, lwd = 0.5)
        lines(pred_range, preds$pred, col = gray(0.1), lwd = 1.2)
        lines(pred_range, preds$ci.lb, col = gray(0.6), lty = 2, lwd = 0.8)
        lines(pred_range, preds$ci.ub, col = gray(0.6), lty = 2, lwd = 0.8)
        abline(h = 0, lty = 3, col = "gray40", lwd = 0.8)
        axis(1, at = x_ticks, tck = -0.02, mgp = c(1.5, 0.2, 0), cex.axis = 1.0, lwd = 1.2, font = 1)
        axis(2, at = y_ticks, tck = -0.02, mgp = c(1.5, 0.4, 0), cex.axis = 1.0, las = 1, lwd = 1.2, font = 1)
        box(lwd = 1.2)
        beta_val <- round(model_add$b[2], 3)
        p_val <- model_add$pval[2]
        n_val <- nrow(snow_add)
        p_text <- if(p_val < 0.001) "<0.001" else paste0("=", round(p_val, 3))
        lab <- bquote(italic(β) == .(beta_val) * ", " * italic(p) * .(p_text) * ", n=" * .(n_val))
        text(x = mean(xlim), y = ylim[2] - 0.05*(ylim[2]-ylim[1]),
             labels = lab, cex = 1.1, font = 2, adj = c(0.5, 1), family = "Arial")
        dev.off()
        cat("  Increased-snowpack-thickness plot saved:", png_file_add, "\n")
      }
      # Decreased-snowpack-thickness plot
      if(!is.null(model_remove) && nrow(snow_remove) >= 3) {
        png_file_remove <- file.path(variable_output_path, paste0("decreased_snowpack_thickness_", safe_name, ".png"))
        png(png_file_remove, width = img_width, height = img_height, units = "px", res = 600)
        par(family = "Arial", oma = c(1.5, 1.5, 0, 0), mar = c(1.5, 2.5, 1.0, 0.5),
            mgp = c(1.5, 0.2, 0), tcl = -0.02, cex.axis = 1.0, las = 1, lwd = 1.2,
            font = 1, pty = "s")
        x_data <- snow_remove[[col_name]]; y_data <- snow_remove$yi; vi_data <- snow_remove$vi
        if(length(vi_data) > 0 && all(!is.na(vi_data)) && all(vi_data > 0)) {
          inv_vi <- 1/vi_data
          if(length(inv_vi) > 1 && var(inv_vi) > 0) {
            point_sizes <- 0.5 + 2.0 * (inv_vi - min(inv_vi)) / (max(inv_vi) - min(inv_vi))
          } else { point_sizes <- rep(1.5, length(x_data)) }
        } else { point_sizes <- rep(1.5, length(x_data)) }
        x_info <- get_ticks(min(x_data), max(x_data)); x_ticks <- x_info$ticks; xlim <- x_info$range
        pred_range <- seq(xlim[1], xlim[2], length = 100)
        preds <- predict(model_remove, newmods = pred_range)
        y_info <- get_regression_y_info(c(y_data, preds$ci.lb, preds$ci.ub), margin_ratio = 0.05)
        y_ticks <- y_info$ticks
        ylim <- y_info$range
        plot(NA, xlim = xlim, ylim = ylim, xlab = "", ylab = "", axes = FALSE, lwd = 1.2, yaxs = "i")
        polygon(c(pred_range, rev(pred_range)), c(preds$ci.lb, rev(preds$ci.ub)),
                col = rgb(233,243,249,maxColorValue=255), border = NA)
        points(x_data, y_data, pch = 21, bg = rgb(143,196,222,maxColorValue=255),
               col = gray(0.2), cex = point_sizes, lwd = 0.5)
        lines(pred_range, preds$pred, col = gray(0.1), lwd = 1.2)
        lines(pred_range, preds$ci.lb, col = gray(0.6), lty = 2, lwd = 0.8)
        lines(pred_range, preds$ci.ub, col = gray(0.6), lty = 2, lwd = 0.8)
        abline(h = 0, lty = 3, col = "gray40", lwd = 0.8)
        axis(1, at = x_ticks, tck = -0.02, mgp = c(1.5, 0.2, 0), cex.axis = 1.0, lwd = 1.2, font = 1)
        axis(2, at = y_ticks, tck = -0.02, mgp = c(1.5, 0.4, 0), cex.axis = 1.0, las = 1, lwd = 1.2, font = 1)
        box(lwd = 1.2)
        beta_val <- round(model_remove$b[2], 3)
        p_val <- model_remove$pval[2]
        n_val <- nrow(snow_remove)
        p_text <- if(p_val < 0.001) "<0.001" else paste0("=", round(p_val, 3))
        lab <- bquote(italic(β) == .(beta_val) * ", " * italic(p) * .(p_text) * ", n=" * .(n_val))
        text(x = mean(xlim), y = ylim[2] - 0.05*(ylim[2]-ylim[1]),
             labels = lab, cex = 1.1, font = 2, adj = c(0.5, 1), family = "Arial")
        dev.off()
        cat("  Decreased-snowpack-thickness plot saved:", png_file_remove, "\n")
      }
      par(par_original)
    }, silent = FALSE)
    
    # ---- 2. Combined plot -----------------------------------------------------
    try({
      bar_data <- subset(regression_coefficients, Variable == variable_name)
      if(nrow(bar_data) == 0) { has_bar <- FALSE } else {
        has_bar <- TRUE
        bar_data <- bar_data[order(bar_data$Treatment == "Increased snowpack thickness", decreasing = TRUE), ]
        beta_vals <- bar_data$Beta
      }
      
      # Common y-axis range for regression plots: calculate confidence intervals across the complete plotted x-axis ranges
      x_add <- if(!is.null(model_add)) snow_add[[col_name]] else numeric(0)
      x_rem <- if(!is.null(model_remove)) snow_remove[[col_name]] else numeric(0)
      y_add <- if(!is.null(model_add)) snow_add$yi else numeric(0)
      y_rem <- if(!is.null(model_remove)) snow_remove$yi else numeric(0)
      all_y <- c(y_add, y_rem)
      all_pred_ci <- numeric(0)
      xlim_add_common <- NULL
      xlim_rem_common <- NULL
      if(!is.null(model_add)) {
        x_info_add_common <- get_ticks(min(x_add), max(x_add))
        xlim_add_common <- x_info_add_common$range
        pred_x_add <- seq(xlim_add_common[1], xlim_add_common[2], length = 100)
        p_add <- predict(model_add, newmods = pred_x_add)
        all_pred_ci <- c(all_pred_ci, p_add$ci.lb, p_add$ci.ub)
      }
      if(!is.null(model_remove)) {
        x_info_rem_common <- get_ticks(min(x_rem), max(x_rem))
        xlim_rem_common <- x_info_rem_common$range
        pred_x_rem <- seq(xlim_rem_common[1], xlim_rem_common[2], length = 100)
        p_rem <- predict(model_remove, newmods = pred_x_rem)
        all_pred_ci <- c(all_pred_ci, p_rem$ci.lb, p_rem$ci.ub)
      }
      y_info <- get_regression_y_info(c(all_y, all_pred_ci), margin_ratio = 0.05)
      y_ticks <- y_info$ticks
      ylim <- y_info$range
      
      # Bar-chart y-axis: three equally spaced ticks with equal numerical intervals and one decimal place
      if(has_bar) {
        bar_axis <- get_bar_axis(beta_vals, margin_ratio = 0.08)
        bar_y_ticks <- bar_axis$ticks
        bar_ylim <- bar_axis$range
      } else {
        bar_y_ticks <- c(-0.1, 0.0, 0.1)
        bar_ylim <- c(-0.1, 0.1)
      }
      
      png_file_combo <- file.path(combo_path, paste0("regression_and_beta_bar_", safe_name, ".png"))
      png(png_file_combo, width = 3600, height = 1440, units = "px", res = 600)
      par(family = "Arial", oma = c(1.5, 2.5, 0.6, 0.8))
      layout(matrix(1:3, nrow = 1), widths = c(1, 1, 0.5))
      
      # Panel 1: Increased-snowpack-thickness regression
      par(mar = c(1.2, 0.6, 0, 0), mgp = c(1.5, 0.6, 0),
          tcl = -0.02, cex.axis = 1.2, las = 1, lwd = 1.2, font = 1)
      if(!is.null(model_add) && nrow(snow_add) >= 3) {
        x_data <- snow_add[[col_name]]; y_data <- snow_add$yi
        x_info <- get_ticks(min(x_data), max(x_data)); x_ticks <- x_info$ticks; xlim_add <- x_info$range
        pred_x <- seq(xlim_add[1], xlim_add[2], length = 100)
        preds <- predict(model_add, newmods = pred_x)
        plot(NA, xlim = xlim_add, ylim = ylim, xlab = "", ylab = "", axes = FALSE, lwd = 1.2, yaxs = "i")
        polygon(c(pred_x, rev(pred_x)), c(preds$ci.lb, rev(preds$ci.ub)),
                col = rgb(254,241,231,maxColorValue=255), border = NA)
        lines(pred_x, preds$pred, col = gray(0.1), lwd = 1.2)
        lines(pred_x, preds$ci.lb, col = gray(0.6), lty = 2, lwd = 0.8)
        lines(pred_x, preds$ci.ub, col = gray(0.6), lty = 2, lwd = 0.8)
        vi <- snow_add$vi
        if(length(vi) > 0 && all(!is.na(vi)) && all(vi > 0)) {
          inv_vi <- 1/vi
          if(length(inv_vi) > 1 && var(inv_vi) > 0) {
            point_sizes <- 1 + 3 * (inv_vi - min(inv_vi)) / (max(inv_vi) - min(inv_vi))
          } else { point_sizes <- rep(2.5, length(x_data)) }
        } else { point_sizes <- rep(2.5, length(x_data)) }
        points(x_data, y_data, pch = 21, bg = rgb(255,192,127,maxColorValue=255),
               col = gray(0.2), cex = point_sizes, lwd = 0.5)
        abline(h = 0, lty = 3, col = "gray40", lwd = 0.8)
        axis(1, at = x_ticks, tck = -0.02, mgp = c(1.5, 0.4, 0), cex.axis = 1.2, lwd = 1.2, font = 1)
        axis(2, at = y_ticks, tck = -0.02, mgp = c(1.5, 0.6, 0), cex.axis = 1.2, las = 1, lwd = 1.2, font = 1)
        box(lwd = 1.2)
        beta_val <- round(model_add$b[2], 3); p_val <- model_add$pval[2]; n_val <- nrow(snow_add)
        if(p_val < 0.001) {
          lab <- bquote(italic(β) == .(beta_val) ~ "," ~ italic(p) < 0.001 ~ "," ~ n == .(n_val))
        } else {
          lab <- bquote(italic(β) == .(beta_val) ~ "," ~ italic(p) == .(round(p_val, 3)) ~ "," ~ n == .(n_val))
        }
        text(x = mean(xlim_add), y = ylim[2] - 0.02*(ylim[2]-ylim[1]),
             labels = lab, cex = 1.3, font = 1, adj = c(0.5, 1), family = "Arial")
      } else {
        plot(NA, xlim = c(0,1), ylim = ylim, axes = FALSE, xlab = "", ylab = "", yaxs = "i")
        box(); text(0.5, 0.5, "No increased-snowpack-thickness data", cex = 1.6, font = 1, family = "Arial")
      }
      
      # Panel 2: Decreased-snowpack-thickness regression
      par(mar = c(1.2, 0, 0, 0.6), mgp = c(1.5, 0.6, 0),
          tcl = -0.02, cex.axis = 1.2, las = 1, lwd = 1.2, font = 1)
      if(!is.null(model_remove) && nrow(snow_remove) >= 3) {
        x_data <- snow_remove[[col_name]]; y_data <- snow_remove$yi
        x_info <- get_ticks(min(x_data), max(x_data)); x_ticks <- x_info$ticks; xlim_rem <- x_info$range
        pred_x <- seq(xlim_rem[1], xlim_rem[2], length = 100)
        preds <- predict(model_remove, newmods = pred_x)
        plot(NA, xlim = xlim_rem, ylim = ylim, xlab = "", ylab = "", axes = FALSE, lwd = 1.2, yaxs = "i")
        polygon(c(pred_x, rev(pred_x)), c(preds$ci.lb, rev(preds$ci.ub)),
                col = rgb(233,243,249,maxColorValue=255), border = NA)
        lines(pred_x, preds$pred, col = gray(0.1), lwd = 1.2)
        lines(pred_x, preds$ci.lb, col = gray(0.6), lty = 2, lwd = 0.8)
        lines(pred_x, preds$ci.ub, col = gray(0.6), lty = 2, lwd = 0.8)
        vi <- snow_remove$vi
        if(length(vi) > 0 && all(!is.na(vi)) && all(vi > 0)) {
          inv_vi <- 1/vi
          if(length(inv_vi) > 1 && var(inv_vi) > 0) {
            point_sizes <- 1 + 3 * (inv_vi - min(inv_vi)) / (max(inv_vi) - min(inv_vi))
          } else { point_sizes <- rep(2.5, length(x_data)) }
        } else { point_sizes <- rep(2.5, length(x_data)) }
        points(x_data, y_data, pch = 21, bg = rgb(143,196,222,maxColorValue=255),
               col = gray(0.2), cex = point_sizes, lwd = 0.5)
        abline(h = 0, lty = 3, col = "gray40", lwd = 0.8)
        axis(1, at = x_ticks, tck = -0.02, mgp = c(1.5, 0.4, 0), cex.axis = 1.2, lwd = 1.2, font = 1)
        box(lwd = 1.2)
        beta_val <- round(model_remove$b[2], 3); p_val <- model_remove$pval[2]; n_val <- nrow(snow_remove)
        if(p_val < 0.001) {
          lab <- bquote(italic(β) == .(beta_val) ~ "," ~ italic(p) < 0.001 ~ "," ~ n == .(n_val))
        } else {
          lab <- bquote(italic(β) == .(beta_val) ~ "," ~ italic(p) == .(round(p_val, 3)) ~ "," ~ n == .(n_val))
        }
        text(x = mean(xlim_rem), y = ylim[2] - 0.02*(ylim[2]-ylim[1]),
             labels = lab, cex = 1.3, font = 1, adj = c(0.5, 1), family = "Arial")
      } else {
        plot(NA, xlim = c(0,1), ylim = ylim, axes = FALSE, xlab = "", ylab = "", yaxs = "i")
        box(); text(0.5, 0.5, "No decreased-snowpack-thickness data", cex = 1.6, font = 1, family = "Arial")
      }
      
      # Panel 3: Beta bar chart
      par(mar = c(1.2, 0, 0, 3.0), mgp = c(1.5, 0.6, 0),
          tcl = -0.02, cex.axis = 1.2, las = 1, lwd = 1.2, font = 1)
      if(has_bar) {
        n_bars <- nrow(bar_data)
        if(n_bars == 1) { x_vals <- 1.5; bar_width <- 0.6 } else { x_vals <- c(1,2); bar_width <- 0.6 }
        plot(NA, xlim = c(0.5, 2.5), ylim = bar_ylim, xlab = "", ylab = "", axes = FALSE, lwd = 1.2, yaxs = "i")
        for(j in seq_along(x_vals)) {
          x_pos <- x_vals[j]; beta <- beta_vals[j]
          col_bar <- if(bar_data$Treatment[j] == "Increased snowpack thickness") rgb(255,192,127,maxColorValue=255) else rgb(143,196,222,maxColorValue=255)
          rect(x_pos - bar_width/2, 0, x_pos + bar_width/2, beta, col = col_bar, border = gray(0.2), lwd = 0.5)
        }
        abline(h = 0, lty = 3, col = "gray40", lwd = 0.8)
        axis(4, at = bar_y_ticks, labels = sprintf("%.1f", bar_y_ticks),
             tck = -0.02, mgp = c(1.5, 0.6, 0), cex.axis = 1.2, las = 1, lwd = 1.2, font = 1)
        box(lwd = 1.2)
      } else {
        plot(NA, xlim = c(0,1), ylim = bar_ylim, axes = FALSE, xlab = "", ylab = "", yaxs = "i")
        box(); text(0.5, 0.5, "No beta data", cex = 1.6, font = 1, family = "Arial")
      }
      par(oma = c(0, 0, 0, 0))
      dev.off()
      cat("  Combined regression and beta-bar plot saved:", png_file_combo, "\n")
    }, silent = FALSE)
    
    # Resume log output
    sink(log_file, split = TRUE, append = TRUE)
  }  # end if model exists
  
  sink()
  cat("\nVariable", variable_name, "was analyzed successfully.\n")
}  # end for loop

# ============================
# Save global results and create the summary beta bar chart
# ============================
if(nrow(regression_coefficients) > 0) {
  write.csv(regression_coefficients, file.path(main_output_path, "regression_coefficients_summary.csv"), row.names = FALSE)
  if(require(openxlsx)) {
    wb <- createWorkbook()
    addWorksheet(wb, "Coefficients")
    writeData(wb, "Coefficients", regression_coefficients)
    saveWorkbook(wb, file.path(main_output_path, "regression_coefficients_summary.xlsx"), overwrite = TRUE)
    cat("\nThe regression coefficient summary was saved in Excel format.\n")
  }
  if(require(ggplot2)) {
    df_bar <- regression_coefficients
    df_bar$Variable <- factor(df_bar$Variable, levels = variables_list[variables_list %in% df_bar$Variable])
    p_bar <- ggplot(df_bar, aes(x = Variable, y = Beta, fill = Treatment)) +
      geom_bar(stat = "identity", position = position_dodge(0.9), width = 0.7) +
      geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper),
                    position = position_dodge(0.9), width = 0.2, size = 0.5) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", size = 0.5) +
      scale_fill_manual(values = c("Increased snowpack thickness" = "#FFC07F", "Decreased snowpack thickness" = "#8FC4DE")) +
      labs(x = "Variable", y = expression("Regression coefficient (" * beta * ")"),
           title = "Regression coefficients for the effects of each variable on soil temperature") +
      theme_minimal(base_family = "Arial") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9, face = "bold"),
            axis.text.y = element_text(size = 9, face = "bold"),
            axis.title = element_text(size = 10, face = "bold"),
            legend.title = element_text(size = 9, face = "bold"),
            legend.text = element_text(size = 8, face = "bold"),
            plot.title = element_text(hjust = 0.5, size = 11, face = "bold"))
    ggsave(file.path(main_output_path, "beta_barplot.png"), p_bar, width = 14, height = 8, dpi = 300)
    cat("The beta summary bar chart was saved: beta_barplot.png\n")
  }
}
if(nrow(all_results_summary) > 0) {
  write.csv(all_results_summary, file.path(main_output_path, "all_results.csv"), row.names = FALSE)
}

cat("\n========== Analysis completed ==========\n")
cat("Number of successfully fitted models:", success_count, "\n")
cat("Separate plots were saved in the corresponding variable subdirectories.\n")
cat("Combined regression and beta-bar plots were saved in:", combo_path, "\n")
cat("The beta summary bar chart was saved in the main output directory.\n")
save.image(file.path(main_output_path, "analysis_workspace.RData"))
cat("The R workspace was saved.\n")