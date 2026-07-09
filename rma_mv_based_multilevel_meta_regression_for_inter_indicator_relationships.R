# rma.mv-Based Multilevel Meta-Regression for Inter-Indicator Relationships
# Complete revised analysis code - SOC regression analysis using the rma.mv model

# 1. Install and load required packages --------------------------------------------------------
cat("Step 1: Loading required packages...\n")

required_packages <- c("metafor", "readxl", "multcomp", "openxlsx", "dplyr", "purrr", "ggplot2")
new_packages <- required_packages[!(required_packages %in% installed.packages()[, "Package"])]
if(length(new_packages)) install.packages(new_packages)

library(metafor)
library(readxl)
library(multcomp)
library(openxlsx)
library(dplyr)
library(purrr)
library(ggplot2)

cat("All packages loaded successfully.\n\n")

# 2. Data import and preprocessing ---------------------------------------------------------
cat("Step 2: Data import and preprocessing...\n")

# Read data - updated input file path
d2 <- read_excel("C:/Users/Administrator/Desktop/25meta/regression_analysis/SOC_regression_analysis_effect_sizes.xlsx")

# Key step: convert data types
cat("\n=== Data Type Conversion ===\n")

# Check and convert the yi column, which is the effect size of SOC
if("yi" %in% names(d2)) {
  cat("Converting the yi column to numeric type...\n")
  d2$yi <- as.numeric(as.character(d2$yi))
  cat("Non-NA yi values:", sum(!is.na(d2$yi)), "/", nrow(d2), "\n")
  cat("Summary statistics for yi:\n")
  print(summary(d2$yi))
} else {
  stop("The yi column does not exist!")
}

# Check and convert the V column
if("V" %in% names(d2)) {
  cat("\nConverting the V column to numeric type...\n")
  d2$V <- as.numeric(as.character(d2$V))
  cat("Non-NA V values:", sum(!is.na(d2$V)), "/", nrow(d2), "\n")
  cat("Summary statistics for V:\n")
  print(summary(d2$V))
} else {
  stop("The V column does not exist!")
}

# Check and convert the vi column
if("vi" %in% names(d2)) {
  cat("\nConverting the vi column to numeric type...\n")
  d2$vi <- as.numeric(as.character(d2$vi))
  cat("Non-NA vi values:", sum(!is.na(d2$vi)), "/", nrow(d2), "\n")
  cat("Summary statistics for vi:\n")
  print(summary(d2$vi))
} else {
  cat("The vi column does not exist; the V column will be used instead.\n")
  d2$vi <- d2$V
}

# Check the Treatment_1 column
if("Treatment_1" %in% names(d2)) {
  cat("\nDistribution of Treatment_1:\n")
  treatment_counts <- table(d2$Treatment_1, useNA = "always")
  print(treatment_counts)
} else {
  stop("The Treatment_1 column does not exist!")
}

# Check required random-effect columns
required_random_cols <- c("Study_id", "Row_id", "Sampling_season", "StandAge_year", "Soil_depth")
cat("\n=== Checking Random-Effect Columns ===\n")
for(col in required_random_cols) {
  if(col %in% names(d2)) {
    cat(sprintf("%s exists\n", col))
  } else {
    cat(sprintf("%s does not exist; a default value will be created\n", col))
    d2[[col]] <- "Default"
  }
}

# 3. Define the variable list and column-name mapping ---------------------------------------------------
cat("\nStep 3: Defining variable mapping...\n")

# Define the list of indicators to be analyzed: x-axis variables, i.e., factors affecting SOC
variables_list <- c(
  "SOC/TN", "Dissolved_organic_carbon", "Total_carbon", "Total_nitrogen", 
  "C/N", "N/P", "Total_phosphorus", "Microbial_biomass_carbon", 
  "Microbial_biomass_nitrogen", "Microbial_biomass_phosphorus", 
  "Soil_water_content", "pH", "Soil_temperature", "Dissolved_organic_nitrogen", 
  "Nitrate_nitrogen", "Ammonium_nitrogen", "Belowground_biomass", 
  "Available_phosphorus", "CO2", "β-glucosidase", "NAG", "PER", "PPO", "ACP", 
  "β-xylosidase", "Urease", "Invertase", "Active_layer_thickness",
  "DIN"
)

# Column-name mapping: these are x-axis variables, while the y-axis is SOC (yi)
column_mapping <- list(
  "SOC/TN" = "SOC/TN",
  "Dissolved_organic_carbon" = "Dissolved_organic_carbon",
  "Total_carbon" = "Total_carbon",
  "Total_nitrogen" = "Total_nitrogen",
  "C/N" = "C/N",
  "N/P" = "N/P",
  "Total_phosphorus" = "Total_phosphorus",
  "Microbial_biomass_carbon" = "Microbial_biomass_carbon",
  "Microbial_biomass_nitrogen" = "Microbial_biomass_nitrogen",
  "Microbial_biomass_phosphorus" = "Microbial_biomass_phosphorus",
  "Soil_water_content" = "Soil_water_content",
  "pH" = "pH",
  "Soil_temperature" = "Soil_temperature",
  "Dissolved_organic_nitrogen" = "Dissolved_organic_nitrogen",
  "Nitrate_nitrogen" = "Nitrate_nitrogen",
  "Ammonium_nitrogen" = "Ammonium_nitrogen",
  "Belowground_biomass" = "Belowground_biomass",
  "Available_phosphorus" = "Available_phosphorus",
  "CO2" = "CO2",
  "β-glucosidase" = "β-glucosidase",
  "NAG" = "NAG",
  "PER" = "PER",
  "PPO" = "PPO",
  "ACP" = "ACP",
  "β-xylosidase" = "β-xylosidase",
  "Urease" = "Urease",
  "Invertase" = "Invertase",
  "Active_layer_thickness" = "Active_layer_thickness",
  "DIN" = "DIN"
)

# Validate the mapping
cat("\n=== Column-Name Mapping Validation ===\n")
for(var in variables_list) {
  if(var %in% names(column_mapping)) {
    col_name <- column_mapping[[var]]
    if(col_name %in% names(d2)) {
      cat(sprintf("%-30s -> %s\n", var, col_name))
    } else {
      cat(sprintf("%-30s -> %s (column does not exist)\n", var, col_name))
    }
  } else {
    cat(sprintf("%-30s -> no mapping definition\n", var))
  }
}

# 4. Create the output directory -------------------------------------------------------------
main_output_path <- "C:/Users/Administrator/Desktop/25meta/regression_analysis/SOC_results_2026_03_18"
if (!dir.exists(main_output_path)) {
  dir.create(main_output_path, recursive = TRUE)
  cat("\nCreated main output directory:", main_output_path, "\n")
} else {
  cat("\nMain output directory already exists:", main_output_path, "\n")
}

# Create a data frame for summarized results
all_results_summary <- data.frame()

# ============================
# Simplified column-name lookup function
# ============================
get_column_name <- function(variable_name) {
  if(variable_name %in% names(column_mapping)) {
    col_name <- column_mapping[[variable_name]]
    if(col_name %in% names(d2)) {
      return(col_name)
    }
  }
  return(NULL)
}

# ============================
# Model-fitting function using rma.mv
# ============================
fit_simple_model <- function(data, col_name) {
  if(nrow(data) < 3) {
    cat("  Insufficient data (n < 3)\n")
    return(NULL)
  }

  # Ensure that data are numeric
  data$yi <- as.numeric(as.character(data$yi))
  data$V <- as.numeric(as.character(data$V))
  data$vi <- as.numeric(as.character(data$vi))
  data[[col_name]] <- as.numeric(as.character(data[[col_name]]))

  # Ensure that random-effect columns are factors
  data$Study_id <- as.factor(data$Study_id)
  data$Row_id <- as.factor(data$Row_id)
  data$Sampling_season <- as.factor(data$Sampling_season)
  data$StandAge_year <- as.factor(data$StandAge_year)
  data$Soil_depth <- as.factor(data$Soil_depth)

  # Remove NA values
  data <- data[!is.na(data$yi) & !is.na(data$V) & !is.na(data[[col_name]]), ]

  if(nrow(data) < 3) {
    cat("  Insufficient valid data (n < 3)\n")
    return(NULL)
  }

  cat("  Valid data:", nrow(data), "observations\n")

  tryCatch({
    # Use the rma.mv model with multilevel random effects
    model <- rma.mv(
      yi = data$yi,                    # Effect size of SOC
      V = data$V,                       # Variance of SOC
      mods = ~ data[[col_name]],        # Independent variable, i.e., each indicator
      random = list(~ 1 | Study_id/Row_id,    # Nested random effects: study/row ID
                    ~ 1 | Sampling_season,     # Sampling season
                    ~ 1 | StandAge_year,       # Stand age
                    ~ 1 | Soil_depth),         # Soil depth
      data = data,
      method = "REML",
      test = "t",                        # Use t-test
      control = list(maxiter = 1000)
    )

    cat("  rma.mv model fitted successfully.\n")
    return(model)
  }, error = function(e) {
    cat("  rma.mv model fitting failed:", e$message, "\n")

    # If the full model fails, try a simplified model
    cat("  Trying a simplified model...\n")
    tryCatch({
      # Simplified model: retain only the main random effect
      model_simple <- rma.mv(
        yi = data$yi,
        V = data$V,
        mods = ~ data[[col_name]],
        random = list(~ 1 | Study_id/Row_id),
        data = data,
        method = "REML",
        test = "t",
        control = list(maxiter = 1000)
      )
      cat("  Simplified rma.mv model fitted successfully.\n")
      return(model_simple)
    }, error = function(e2) {
      cat("  Simplified model also failed:", e2$message, "\n")
      return(NULL)
    })
  })
}

# ============================
# Loop through each indicator
# ============================
cat("\n", strrep("=", 80), "\n")
cat("Starting the loop analysis for each indicator using the rma.mv model...\n")
cat("Total number of indicators:", length(variables_list), "\n")
cat(strrep("=", 80), "\n\n")

success_count <- 0

# Create a data frame to store regression coefficients
regression_coefficients <- data.frame()

for (i in seq_along(variables_list)) {
  variable_name <- variables_list[i]

  cat("\n", strrep("-", 60), "\n")
  cat("Analyzing indicator:", i, "/", length(variables_list), "-", variable_name, "\n")
  cat(strrep("-", 60), "\n")

  # Get the column name
  col_name <- get_column_name(variable_name)
  if(is.null(col_name)) {
    cat("Column name could not be found; skipped.\n")
    next
  }

  cat("Column name used:", col_name, "\n")

  # Create the output directory for this indicator
  safe_name <- gsub("[^[:alnum:]_]", "_", variable_name)
  variable_output_path <- file.path(main_output_path, safe_name)
  if (!dir.exists(variable_output_path)) {
    dir.create(variable_output_path, recursive = TRUE)
  }

  # Set the current working directory
  setwd(variable_output_path)

  # Create a log file
  log_file <- file.path(variable_output_path, paste0("analysis_log_", safe_name, ".txt"))
  sink(log_file, split = TRUE)

  cat("Analyzed indicator:", variable_name, "\n")
  cat("Corresponding column name:", col_name, "\n")
  cat("Analysis time:", Sys.time(), "\n\n")

  # ============================
  # Data preprocessing
  # ============================
  cat("=== Data Preprocessing ===\n")

  # Check whether the column exists
  if(!col_name %in% names(d2)) {
    cat("Column", col_name, "does not exist\n")
    sink()
    next
  }

  # Convert the indicator column to numeric
  d2[[col_name]] <- as.numeric(as.character(d2[[col_name]]))

  # Check data validity
  n_valid <- sum(!is.na(d2[[col_name]]) & !is.na(d2$yi) & !is.na(d2$V) & !is.na(d2$vi))
  cat("Total valid data points:", n_valid, "\n")

  if (n_valid < 3) {
    cat("Insufficient valid data\n")
    sink()
    next
  }

  # Split the data
  snow_add <- subset(d2, Treatment_1 == "Increased snowpack thickness" & 
                      !is.na(d2[[col_name]]) & !is.na(yi) & !is.na(V) & !is.na(vi))
  decreased_snowpack_thickness <- subset(d2, Treatment_1 == "Decreased snowpack thickness" & 
                         !is.na(d2[[col_name]]) & !is.na(yi) & !is.na(V) & !is.na(vi))

  cat("Increased snowpack thickness: ", nrow(snow_add), "observations\n")
  cat("Decreased snowpack thickness: ", nrow(decreased_snowpack_thickness), "observations\n")

  # ============================
  # Fit models
  # ============================
  cat("\n=== Model Fitting Using rma.mv ===\n")

  model_add <- NULL
  model_remove <- NULL

  # Fit the Increased snowpack thickness model
  if(nrow(snow_add) >= 3) {
    cat("Fitting the Increased snowpack thickness model...\n")
    model_add <- fit_simple_model(snow_add, col_name)
    if(!is.null(model_add)) {
      cat("Model summary:\n")
      print(summary(model_add))

      # Save the model
      saveRDS(model_add, file.path(variable_output_path, 
                                  paste0("model_snow_add_", safe_name, ".rds")))

      # Save random-effect variance components
      if(length(model_add$sigma2) > 0) {
        random_effects <- data.frame(
          Component = names(model_add$s.names),
          Variance = model_add$sigma2,
          SD = sqrt(model_add$sigma2)
        )
        write.csv(random_effects, 
                  file.path(variable_output_path, 
                           paste0("random_effects_add_", safe_name, ".csv")),
                  row.names = FALSE)
      }
    }
  } else {
    cat("Insufficient data for Increased snowpack thickness (n < 3)\n")
  }

  # Fit the Decreased snowpack thickness model
  if(nrow(decreased_snowpack_thickness) >= 3) {
    cat("\nFitting the Decreased snowpack thickness model...\n")
    model_remove <- fit_simple_model(decreased_snowpack_thickness, col_name)
    if(!is.null(model_remove)) {
      cat("Model summary:\n")
      print(summary(model_remove))

      # Save the model
      saveRDS(model_remove, file.path(variable_output_path, 
                                     paste0("model_decreased_snowpack_thickness_", safe_name, ".rds")))

      # Save random-effect variance components
      if(length(model_remove$sigma2) > 0) {
        random_effects <- data.frame(
          Component = names(model_remove$s.names),
          Variance = model_remove$sigma2,
          SD = sqrt(model_remove$sigma2)
        )
        write.csv(random_effects, 
                  file.path(variable_output_path, 
                           paste0("random_effects_remove_", safe_name, ".csv")),
                  row.names = FALSE)
      }
    }
  } else {
    cat("Insufficient data for Decreased snowpack thickness (n < 3)\n")
  }

  # ============================
  # Result summary
  # ============================
  cat("\n=== Result Summary ===\n")

  variable_results <- data.frame()

  # Increased snowpack thickness results
  if(!is.null(model_add)) {
    # For the rma.mv model, coefficients are stored in the $b matrix and need to be indexed correctly
    beta_value <- model_add$b[2]
    se_value <- model_add$se[2]
    ci_lb <- model_add$ci.lb[2]
    ci_ub <- model_add$ci.ub[2]
    p_value <- model_add$pval[2]

    results_add <- data.frame(
      Variable = variable_name,
      Treatment = "Increased snowpack thickness",
      Beta = beta_value,
      SE = se_value,
      CI_lower = ci_lb,
      CI_upper = ci_ub,
      p_value = p_value,
      QM = model_add$QM,
      QM_p = model_add$QMp,
      n_obs = nrow(snow_add),
      tau2 = ifelse(length(model_add$sigma2) > 0, sum(model_add$sigma2), NA),
      I2 = ifelse(!is.null(model_add$I2), model_add$I2, NA),
      stringsAsFactors = FALSE
    )
    variable_results <- rbind(variable_results, results_add)
    success_count <- success_count + 1

    # Add to the regression-coefficient table
    regression_coefficients <- rbind(regression_coefficients, 
                                     data.frame(
                                       Variable = variable_name,
                                       Treatment = "Increased snowpack thickness",
                                       Beta = round(beta_value, 4),
                                       SE = round(se_value, 4),
                                       CI = paste0("[", round(ci_lb, 4), ", ", round(ci_ub, 4), "]"),
                                       p_value = round(p_value, 4),
                                       n = nrow(snow_add),
                                       stringsAsFactors = FALSE
                                     ))
  }

  # Decreased snowpack thickness results
  if(!is.null(model_remove)) {
    beta_value <- model_remove$b[2]
    se_value <- model_remove$se[2]
    ci_lb <- model_remove$ci.lb[2]
    ci_ub <- model_remove$ci.ub[2]
    p_value <- model_remove$pval[2]

    results_remove <- data.frame(
      Variable = variable_name,
      Treatment = "Decreased snowpack thickness",
      Beta = beta_value,
      SE = se_value,
      CI_lower = ci_lb,
      CI_upper = ci_ub,
      p_value = p_value,
      QM = model_remove$QM,
      QM_p = model_remove$QMp,
      n_obs = nrow(decreased_snowpack_thickness),
      tau2 = ifelse(length(model_remove$sigma2) > 0, sum(model_remove$sigma2), NA),
      I2 = ifelse(!is.null(model_remove$I2), model_remove$I2, NA),
      stringsAsFactors = FALSE
    )
    variable_results <- rbind(variable_results, results_remove)
    success_count <- success_count + 1

    # Add to the regression-coefficient table
    regression_coefficients <- rbind(regression_coefficients, 
                                     data.frame(
                                       Variable = variable_name,
                                       Treatment = "Decreased snowpack thickness",
                                       Beta = round(beta_value, 4),
                                       SE = round(se_value, 4),
                                       CI = paste0("[", round(ci_lb, 4), ", ", round(ci_ub, 4), "]"),
                                       p_value = round(p_value, 4),
                                       n = nrow(decreased_snowpack_thickness),
                                       stringsAsFactors = FALSE
                                     ))
  }

  # Save results
  if(nrow(variable_results) > 0) {
    cat("\nRegression results:\n")
    print(variable_results, row.names = FALSE)
    all_results_summary <- rbind(all_results_summary, variable_results)

    write.csv(variable_results, 
              file.path(variable_output_path, paste0("results_", safe_name, ".csv")), 
              row.names = FALSE)
    cat("Results saved successfully.\n")
  } else {
    cat("No available model results.\n")
  }

  # ============================
  # Generate figures: one 8 cm x 8 cm figure for each indicator, with statistics at the top
  # ============================
  if(!is.null(model_add) || !is.null(model_remove)) {
    cat("\n=== Generating Regression Plots ===\n")

    tryCatch({
      # Define output path
      output_path_tif <- "C:/Users/Administrator/Desktop/25meta/regression_analysis/final_regression_figures/SOC"
      if (!dir.exists(output_path_tif)) {
        dir.create(output_path_tif, recursive = TRUE)
      }

      # Set figure parameters: 8 cm x 8 cm, 300 DPI
      # 8 cm = 8 * 300 / 2.54, approximately 945 pixels
      img_width <- 945
      img_height <- 945

      # Set basic graphical parameters
      par_original <- par(no.readonly = TRUE)

      # ===== Increased snowpack thickness plot =====
      if(!is.null(model_add) && nrow(snow_add) >= 3) {
        # Set TIF output
        tif_file_add <- file.path(output_path_tif, paste0("Increased_snowpack_thickness_", safe_name, ".tif"))
        tiff(tif_file_add, 
             width = img_width, 
             height = img_height, 
             units = "px", 
             res = 300,
             compression = "lzw")

        # Set graphical parameters; axis titles are not shown, so margins can be reduced
        par(family = "serif",
            mar = c(2.5, 2.5, 1.5, 0.5),
            mgp = c(1.5, 0.4, 0),
            tcl = -0.02,
            cex.axis = 0.7,
            font.axis = 2,
            las = 1,
            lwd = 1.2)

        # Prepare data
        x_data <- snow_add[[col_name]]
        y_data <- snow_add$yi
        vi_data <- snow_add$vi

        # Calculate point sizes based on 1/V, ranging from 1 to 4
        if(length(vi_data) > 0 && all(!is.na(vi_data)) && all(vi_data > 0)) {
          inv_vi <- 1/vi_data
          if(length(inv_vi) > 1 && var(inv_vi) > 0) {
            point_sizes <- 1 + 3 * (inv_vi - min(inv_vi)) / (max(inv_vi) - min(inv_vi))
          } else {
            point_sizes <- rep(2.5, length(x_data))
          }
        } else {
          point_sizes <- rep(2.5, length(x_data))
        }

        # Calculate the x-axis range
        min_x <- min(x_data, na.rm = TRUE)
        max_x <- max(x_data, na.rm = TRUE)
        range_x <- max_x - min_x
        x_min <- min_x - range_x * 0.05
        x_max <- max_x + range_x * 0.05

        # Calculate predicted values
        pred_range <- seq(x_min, x_max, length = 100)
        preds <- predict(model_add, newmods = pred_range)

        # Calculate the y-axis range
        y_min <- min(c(y_data, preds$ci.lb), na.rm = TRUE)
        y_max <- max(c(y_data, preds$ci.ub), na.rm = TRUE)
        y_range <- y_max - y_min
        y_lim <- c(y_min - y_range * 0.05, y_max + y_range * 0.15)

        # Create the base plot without axis labels or axis titles
        plot(x_data, y_data,
             type = "n",
             xlab = "",
             ylab = "",
             main = "",
             xlim = c(x_min, x_max),
             ylim = y_lim,
             axes = FALSE,
             lwd = 1.2)

        # Add confidence interval
        confidence_fill_add <- rgb(254/255, 241/255, 231/255)
        polygon(c(pred_range, rev(pred_range)),
                c(preds$ci.lb, rev(preds$ci.ub)),
                col = confidence_fill_add,
                border = NA)

        # Add data points
        point_inner_add <- rgb(255/255, 192/255, 127/255)
        point_border_add <- rgb(50/255, 50/255, 50/255)
        points(x_data, y_data,
               pch = 21,
               bg = point_inner_add,
               col = point_border_add,
               cex = point_sizes,
               lwd = 0.5)

        # Add the regression line and confidence-interval lines
        line_color_add <- rgb(70/255, 70/255, 70/255)
        lines(pred_range, preds$pred, lwd = 1.2, col = line_color_add)
        lines(pred_range, preds$ci.lb, lwd = 0.8, lty = "dashed", col = line_color_add)
        lines(pred_range, preds$ci.ub, lwd = 0.8, lty = "dashed", col = line_color_add)

        # Add the zero line
        abline(h = 0, lty = "dotted", col = "gray40", lwd = 0.8)

        # Add axes
        axis(1, tck = -0.02,
             mgp = c(1.5, 0, 0),
             cex.axis = 0.7,
             lwd = 1.2,
             font = 2)

        axis(2, tck = -0.02,
             mgp = c(1.5, 0.4, 0),
             cex.axis = 0.7,
             las = 1,
             lwd = 1.2,
             font = 2)

        # Add a thicker border
        box(lwd = 1.2)

        # Add statistics at the top of the plot
        beta_val <- round(model_add$b[2], 3)
        p_val <- model_add$pval[2]
        n_val <- nrow(snow_add)

        # Format p-values
        if(p_val < 0.001) {
          p_text <- paste0("p<0.001")
        } else {
          p_text <- paste0("p=", round(p_val, 3))
        }

        text_x <- x_min + (x_max - x_min) * 0.5
        text_y <- y_lim[2] - (y_lim[2] - y_lim[1]) * 0.04

        stat_label <- bquote(italic(β) == .(beta_val) * ", " * italic(p) * .(substr(p_text, 2, nchar(p_text))) * ", n=" * .(n_val))

        text(text_x, text_y, 
             labels = stat_label,
             cex = 0.9,
             family = "serif",
             font = 2)

        dev.off()
        cat("Increased snowpack thickness plot saved successfully:", tif_file_add, "\n")
      }

      # ===== Decreased snowpack thickness plot =====
      if(!is.null(model_remove) && nrow(decreased_snowpack_thickness) >= 3) {
        # Set TIF output
        tif_file_remove <- file.path(output_path_tif, paste0("Decreased_snowpack_thickness_", safe_name, ".tif"))
        tiff(tif_file_remove, 
             width = img_width, 
             height = img_height, 
             units = "px", 
             res = 300,
             compression = "lzw")

        # Set graphical parameters
        par(family = "serif",
            mar = c(2.5, 2.5, 1.5, 0.5),
            mgp = c(1.5, 0.4, 0),
            tcl = -0.02,
            cex.axis = 0.7,
            font.axis = 2,
            las = 1,
            lwd = 1.2)

        # Prepare data
        x_data <- decreased_snowpack_thickness[[col_name]]
        y_data <- decreased_snowpack_thickness$yi
        vi_data <- decreased_snowpack_thickness$vi

        # Calculate point sizes based on 1/V, ranging from 1 to 4
        if(length(vi_data) > 0 && all(!is.na(vi_data)) && all(vi_data > 0)) {
          inv_vi <- 1/vi_data
          if(length(inv_vi) > 1 && var(inv_vi) > 0) {
            point_sizes <- 1 + 3 * (inv_vi - min(inv_vi)) / (max(inv_vi) - min(inv_vi))
          } else {
            point_sizes <- rep(2.5, length(x_data))
          }
        } else {
          point_sizes <- rep(2.5, length(x_data))
        }

        # Calculate the x-axis range
        min_x <- min(x_data, na.rm = TRUE)
        max_x <- max(x_data, na.rm = TRUE)
        range_x <- max_x - min_x
        x_min <- min_x - range_x * 0.05
        x_max <- max_x + range_x * 0.05

        # Calculate predicted values
        pred_range <- seq(x_min, x_max, length = 100)
        preds <- predict(model_remove, newmods = pred_range)

        # Calculate the y-axis range
        y_min <- min(c(y_data, preds$ci.lb), na.rm = TRUE)
        y_max <- max(c(y_data, preds$ci.ub), na.rm = TRUE)
        y_range <- y_max - y_min
        y_lim <- c(y_min - y_range * 0.05, y_max + y_range * 0.15)

        # Create the base plot
        plot(x_data, y_data,
             type = "n",
             xlab = "",
             ylab = "",
             main = "",
             xlim = c(x_min, x_max),
             ylim = y_lim,
             axes = FALSE,
             lwd = 1.2)

        # Add confidence interval
        confidence_fill_remove <- rgb(233/255, 243/255, 249/255)
        polygon(c(pred_range, rev(pred_range)),
                c(preds$ci.lb, rev(preds$ci.ub)),
                col = confidence_fill_remove,
                border = NA)

        # Add data points
        point_inner_remove <- rgb(143/255, 196/255, 222/255)
        point_border_remove <- rgb(50/255, 50/255, 50/255)
        points(x_data, y_data,
               pch = 21,
               bg = point_inner_remove,
               col = point_border_remove,
               cex = point_sizes,
               lwd = 0.5)

        # Add the regression line and confidence-interval lines
        line_color_remove <- rgb(70/255, 70/255, 70/255)
        lines(pred_range, preds$pred, lwd = 1.2, col = line_color_remove)
        lines(pred_range, preds$ci.lb, lwd = 0.8, lty = "dashed", col = line_color_remove)
        lines(pred_range, preds$ci.ub, lwd = 0.8, lty = "dashed", col = line_color_remove)

        # Add the zero line
        abline(h = 0, lty = "dotted", col = "gray40", lwd = 0.8)

        # Add axes
        axis(1, tck = -0.02,
             mgp = c(1.5, 0, 0),
             cex.axis = 0.7,
             lwd = 1.2,
             font = 2)

        axis(2, tck = -0.02,
             mgp = c(1.5, 0.4, 0),
             cex.axis = 0.7,
             las = 1,
             lwd = 1.2,
             font = 2)

        # Add a thicker border
        box(lwd = 1.2)

        # Add statistics at the top of the plot
        beta_val <- round(model_remove$b[2], 3)
        p_val <- model_remove$pval[2]
        n_val <- nrow(decreased_snowpack_thickness)

        # Format p-values
        if(p_val < 0.001) {
          p_text <- paste0("p<0.001")
        } else {
          p_text <- paste0("p=", round(p_val, 3))
        }

        text_x <- x_min + (x_max - x_min) * 0.5
        text_y <- y_lim[2] - (y_lim[2] - y_lim[1]) * 0.04

        stat_label <- bquote(italic(β) == .(beta_val) * ", " * italic(p) * .(substr(p_text, 2, nchar(p_text))) * ", n=" * .(n_val))

        text(text_x, text_y, 
             labels = stat_label,
             cex = 0.9,
             family = "serif",
             font = 2)

        dev.off()
        cat("Decreased snowpack thickness plot saved successfully:", tif_file_remove, "\n")
      }

      # Restore original graphical parameters
      par(par_original)

    }, error = function(e) {
      cat("Plot generation failed:", e$message, "\n")
    })
  }

  # Close the log file
  sink()

  cat("\nIndicator", variable_name, "analysis completed.\n")
}

# ============================
# Save the regression-coefficient table
# ============================
if(nrow(regression_coefficients) > 0) {
  # Sort by variable name
  regression_coefficients <- regression_coefficients[order(regression_coefficients$Variable, regression_coefficients$Treatment), ]

  # Save as a CSV file
  output_path_tif <- "C:/Users/Administrator/Desktop/25meta/regression_analysis/final_regression_figures/SOC"
  write.csv(regression_coefficients, 
            file.path(output_path_tif, "regression_coefficient_table.csv"), 
            row.names = FALSE)

  # Save as an Excel file if openxlsx is installed
  if(require(openxlsx)) {
    # Create a workbook
    wb <- createWorkbook()
    addWorksheet(wb, "Regression_Coefficients")

    # Write data
    writeData(wb, "Regression_Coefficients", regression_coefficients)

    # Set column widths
    setColWidths(wb, "Regression_Coefficients", cols = 1:ncol(regression_coefficients), widths = "auto")

    # Add style
    headerStyle <- createStyle(fontSize = 11, fontColour = "#000000", 
                               textDecoration = "bold", halign = "center")
    addStyle(wb, "Regression_Coefficients", headerStyle, rows = 1, cols = 1:ncol(regression_coefficients))

    # Save the file
    saveWorkbook(wb, file.path(output_path_tif, "regression_coefficient_table.xlsx"), overwrite = TRUE)
    cat("\nRegression-coefficient table saved in Excel format.\n")
  }

  cat("\nRegression-coefficient table saved to:", file.path(output_path_tif, "regression_coefficient_table.csv"), "\n")

  # Display table content
  cat("\n=== Standardized Regression-Coefficient Table ===\n")
  print(regression_coefficients)
}

# ============================
# Save the summary of all results
# ============================
cat("\n", strrep("=", 80), "\n")
cat("All indicator analyses completed!\n")
cat(strrep("=", 80), "\n\n")

cat("Number of successfully fitted models:", success_count, "\n")

# Save summarized results
if(nrow(all_results_summary) > 0) {
  setwd(main_output_path)
  write.csv(all_results_summary, "all_variables_results_summary.csv", row.names = FALSE)
  cat("Summary of all results saved to: all_variables_results_summary.csv\n")

  # Display successfully analyzed indicators
  cat("\nSuccessfully analyzed indicators:\n")
  successful_vars <- unique(all_results_summary$Variable)
  for(var in successful_vars) {
    treatments <- all_results_summary$Treatment[all_results_summary$Variable == var]
    cat("- ", var, ": ", paste(unique(treatments), collapse = ", "), "\n")
  }

  # Create summary statistics
  summary_stats <- all_results_summary %>%
    group_by(Variable) %>%
    summarise(
      n_models = n(),
      n_significant = sum(p_value < 0.05, na.rm = TRUE),
      mean_beta = mean(Beta, na.rm = TRUE),
      mean_tau2 = mean(tau2, na.rm = TRUE),
      .groups = "drop"
    )

  write.csv(summary_stats, "summary_statistics.csv", row.names = FALSE)
  cat("\nSummary statistics saved to: summary_statistics.csv\n")

  # Create summary visualization
  if(require(ggplot2)) {
    # Create a cleaner theme
    my_theme <- theme_minimal(base_family = "serif") +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, size = 9, 
                                   margin = margin(t = -5)),
        axis.text.y = element_text(size = 9, hjust = 1),
        axis.title.x = element_text(size = 10, face = "bold", 
                                    margin = margin(t = 5)),
        axis.title.y = element_text(size = 10, face = "bold", 
                                    margin = margin(r = 5)),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 11,
                                  margin = margin(b = 10)),
        legend.title = element_text(size = 9, face = "bold"),
        legend.text = element_text(size = 8),
        legend.position = "right",
        legend.margin = margin(l = -10),
        panel.grid.major = element_line(color = "gray90", linewidth = 0.3),
        panel.grid.minor = element_blank(),
        panel.spacing = unit(0.5, "lines"),
        plot.margin = margin(10, 10, 10, 10)
      )

    p <- ggplot(all_results_summary, aes(x = Variable, y = Beta, color = Treatment)) +
      geom_point(position = position_dodge(width = 0.5), size = 2) +
      geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper), 
                    width = 0.2, 
                    linewidth = 0.5,
                    position = position_dodge(width = 0.5)) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.5) +
      labs(x = "Indicator", 
           y = expression("Regression coefficient (β): effect on SOC"),
           title = "Regression coefficients of each indicator affecting SOC (rma.mv model)",
           color = "Treatment") +
      scale_color_manual(values = c("Increased snowpack thickness" = "#FFC07F", 
                                    "Decreased snowpack thickness" = "#8FC4DE")) +
      my_theme

    ggsave("all_variables_beta_coefficients.png", p, 
           width = 14, height = 8, dpi = 300, bg = "white")
    cat("Summary visualization saved to: all_variables_beta_coefficients.png\n")

    # Additionally save a plot with random-effect information
    p2 <- ggplot(all_results_summary, aes(x = Variable, y = tau2, fill = Treatment)) +
      geom_bar(stat = "identity", position = position_dodge(width = 0.7), width = 0.6) +
      labs(x = "Indicator", 
           y = expression("Total heterogeneity (τ"^2~")"),
           title = "Total heterogeneity variance of each model",
           fill = "Treatment") +
      scale_fill_manual(values = c("Increased snowpack thickness" = "#FFC07F", 
                                   "Decreased snowpack thickness" = "#8FC4DE")) +
      theme_minimal(base_family = "serif") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))

    ggsave("heterogeneity_tau2.png", p2, width = 14, height = 6, dpi = 300)
    cat("Heterogeneity plot saved to: heterogeneity_tau2.png\n")
  }
} else {
  cat("No models were successfully fitted.\n")
}

# Save the workspace
save.image(file.path(main_output_path, "analysis_workspace.RData"))

cat("\n========== Analysis Completed ==========\n")
cat("Analysis time:", Sys.time(), "\n")
cat("All results have been saved to:", main_output_path, "\n")
cat("\nGenerated file structure:\n")
if(file.exists(main_output_path)) {
  files <- list.files(main_output_path, recursive = TRUE, full.names = FALSE)
  for(file in files) {
    cat("  ", file, "\n")
  }
}
