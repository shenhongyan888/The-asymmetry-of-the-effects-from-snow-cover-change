# ============================
# Overall funnel plots only
# Output PNG files named by indicator
# Treatment labels:
# Increased snowpack thickness
# Decreased snowpack thickness
# Border not bold
# With 3 axis ticks and outward tick marks
# ============================

packages <- c("readxl", "metafor", "ggplot2", "dplyr", "grid")

for (pkg in packages) {
  if (!require(pkg, character.only = TRUE)) install.packages(pkg)
  library(pkg, character.only = TRUE)
}

# ==================== User settings ====================

root_output_dir <- "E:/\u535a\u58eb\u671f\u95f4\u6587\u6863/25meta/\u6f0f\u6597\u56fe7.9"

file_configs <- list(
  list(
    path = "E:/\u535a\u58eb\u671f\u95f4\u6587\u6863/25meta/meta\u6570\u636e/1.xlsx",
    output_subdir = "1",
    file_label = "1"
  ),
  list(
    path = "E:/\u535a\u58eb\u671f\u95f4\u6587\u6863/25meta/\u56de\u5f52\u5206\u6790/SOCTNR\u56de\u5f52\u5206\u6790\u6548\u5e94\u503c.xlsx",
    output_subdir = "SOCTNR",
    file_label = "SOCTNR"
  ),
  list(
    path = "E:/\u535a\u58eb\u671f\u95f4\u6587\u6863/25meta/\u56de\u5f52\u5206\u6790/DIN\u56de\u5f52\u5206\u6790\u6548\u5e94\u503c.xlsx",
    output_subdir = "DIN",
    file_label = "DIN"
  )
)

# Fixed inner panel size, not whole image size
panel_size_cm <- 5

# Output resolution
dpi <- 600

# Font and size
font_family <- "Arial"
font_size <- 11

# Treatment labels
treatment_increased <- "Increased snowpack thickness"
treatment_decreased <- "Decreased snowpack thickness"

# Folder names
treatment_folder_names <- c(
  "Increased snowpack thickness" = "Increased_snowpack_thickness",
  "Decreased snowpack thickness" = "Decreased_snowpack_thickness"
)

# Treatment colors
treatment_colors <- c(
  "Increased snowpack thickness" = rgb(252, 183, 132, maxColorValue = 255),
  "Decreased snowpack thickness" = rgb(145, 171, 210, maxColorValue = 255)
)

# Reserved display labels for possible grouped plots
duration_display_label <- "Experimental duration"
frozen_soil_display_label <- "Frozen soil types"

# ==================== Helper functions ====================

clean_filename <- function(x) {
  x <- as.character(x)
  x[is.na(x) | x == ""] <- "Unknown"
  x <- trimws(x)
  illegal <- "[\\\\/:*?\"<>|]"
  x <- gsub(illegal, "_", x)
  x <- gsub("\\s+", "_", x)
  return(x)
}

to_num <- function(x) {
  suppressWarnings(as.numeric(as.character(x)))
}

find_col <- function(df_names, candidates) {
  df_names_trim <- trimws(df_names)
  df_names_low <- tolower(df_names_trim)
  cand_low <- tolower(candidates)
  
  idx <- match(cand_low, df_names_low)
  idx <- idx[!is.na(idx)]
  
  if (length(idx) == 0) return(NULL)
  return(df_names_trim[idx[1]])
}

normalize_treatment <- function(x) {
  s <- tolower(trimws(as.character(x)))
  
  out <- rep(NA_character_, length(s))
  
  out[grepl(
    "snow addition|addition|add|increase|increased|deepened|increased snowpack thickness|\u589e\u96ea",
    s
  )] <- treatment_increased
  
  out[grepl(
    "snow remove|snow removal|snow reduce|reduce|reduced|remove|removal|decrease|decreased|decreased snowpack thickness|\u51cf\u96ea",
    s
  )] <- treatment_decreased
  
  return(out)
}

make_axis_3_ticks <- function(values) {
  values <- values[is.finite(values)]
  
  if (length(values) == 0) {
    return(list(
      limits = c(0, 2),
      breaks = c(0, 1, 2),
      digits = 0
    ))
  }
  
  min_v <- min(values, na.rm = TRUE)
  max_v <- max(values, na.rm = TRUE)
  
  if (!is.finite(min_v) || !is.finite(max_v)) {
    return(list(
      limits = c(0, 2),
      breaks = c(0, 1, 2),
      digits = 0
    ))
  }
  
  if (min_v == max_v) {
    min_v <- min_v - 0.5
    max_v <- max_v + 0.5
  }
  
  lower_int <- floor(min_v)
  upper_int <- ceiling(max_v)
  int_range <- upper_int - lower_int
  
  if (int_range >= 2) {
    if ((int_range %% 2) != 0) {
      lower_pad <- min_v - lower_int
      upper_pad <- upper_int - max_v
      
      if (upper_pad <= lower_pad) {
        upper_int <- upper_int + 1
      } else {
        lower_int <- lower_int - 1
      }
    }
    
    breaks <- c(
      lower_int,
      (lower_int + upper_int) / 2,
      upper_int
    )
    
    return(list(
      limits = c(lower_int, upper_int),
      breaks = breaks,
      digits = 0
    ))
  }
  
  lower_dec_i <- floor(min_v * 10)
  upper_dec_i <- ceiling(max_v * 10)
  
  if ((upper_dec_i - lower_dec_i) < 2) {
    lower_dec_i <- lower_dec_i - 1
    upper_dec_i <- upper_dec_i + 1
  }
  
  if (((upper_dec_i - lower_dec_i) %% 2) != 0) {
    upper_dec_i <- upper_dec_i + 1
  }
  
  breaks <- c(
    lower_dec_i,
    (lower_dec_i + upper_dec_i) / 2,
    upper_dec_i
  ) / 10
  
  return(list(
    limits = c(lower_dec_i / 10, upper_dec_i / 10),
    breaks = breaks,
    digits = 1
  ))
}

axis_label_fun <- function(digits) {
  force(digits)
  function(x) {
    if (digits == 0) {
      return(as.character(round(x)))
    } else {
      return(sprintf("%.1f", x))
    }
  }
}

translate_indicator_name <- function(x) {
  x <- trimws(as.character(x))
  
  indicator_map <- c(
    "\u571f\u58e4\u6709\u673a\u78b3" = "Soil organic carbon",
    "\u603b\u78b3" = "Total carbon",
    "\u603b\u6c2e" = "Total nitrogen",
    "\u603b\u78f7" = "Total phosphorus",
    "\u53ef\u6eb6\u6027\u6709\u673a\u78b3" = "Dissolved organic carbon",
    "\u53ef\u6eb6\u6027\u6709\u673a\u6c2e" = "Dissolved organic nitrogen",
    "\u5fae\u751f\u7269\u91cf\u78b3" = "Microbial biomass carbon",
    "\u5fae\u751f\u7269\u91cf\u6c2e" = "Microbial biomass nitrogen",
    "\u5fae\u751f\u7269\u91cf\u78f7" = "Microbial biomass phosphorus",
    "\u785d\u6001\u6c2e" = "Nitrate nitrogen",
    "\u94f5\u6001\u6c2e" = "Ammonium nitrogen",
    "\u6c2e\u7d20\u6709\u6548\u6027" = "Nitrogen availability",
    "\u6709\u6548\u78f7" = "Available phosphorus",
    "\u571f\u58e4\u542b\u6c34\u7387" = "Soil water content",
    "\u571f\u58e4\u6e29\u5ea6" = "Soil temperature",
    "\u571f\u58e4\u6e29\u5ea6" = "Soil temperature",
    "\u571f\u58e4pH" = "Soil pH",
    "\u9633\u79bb\u5b50\u4ea4\u6362\u91cf" = "Cation exchange capacity",
    "\u6d3b\u52a8\u5c42\u539a\u5ea6" = "Active layer thickness",
    "\u79ef\u96ea\u878d\u5316\u65e5\u671f" = "Snowmelt date",
    "\u51bb\u7ed3\u5929\u6570" = "Number of frozen days",
    "\u51bb\u878d\u5faa\u73af\u6b21\u6570" = "Number of freeze-thaw cycles",
    "\u51bb\u7ed3\u6df1\u5ea6" = "Frost depth",
    "\u8102\u80aa\u9176" = "Lipase",
    "\u8102\u80aa\u9176" = "Lipase",
    "\u8102\u80aa\u9176\u6d3b\u6027" = "Lipase activity",
    "\u5c3f\u9176" = "Urease",
    "\u8f6c\u5316\u9176" = "Invertase",
    "\u8517\u7cd6\u9176" = "Sucrase",
    "\u78b1\u6027\u78f7\u9178\u9176" = "Alkaline phosphatase",
    "\u9178\u6027\u78f7\u9178\u9176" = "Acid phosphatase",
    "\u7ea4\u7ef4\u7d20\u9176" = "Cellulase",
    "\u8fc7\u6c27\u5316\u7269\u9176" = "Peroxidase",
    "\u591a\u915a\u6c27\u5316\u9176" = "Polyphenol oxidase",
    "\u547c\u5438\u901f\u7387" = "Respiration rate",
    "\u4e8c\u6c27\u5316\u78b3\u901a\u91cf" = "CO2 flux",
    "\u7532\u70f7\u901a\u91cf" = "CH4 flux",
    "\u6c27\u5316\u4e9a\u6c2e\u901a\u91cf" = "N2O flux",
    "\u5730\u4e0b\u751f\u7269\u91cf" = "Belowground biomass",
    "\u5730\u4e0a\u751f\u7269\u91cf" = "Aboveground biomass",
    "\u51cb\u843d\u7269\u5206\u89e3\u901f\u7387" = "Litter decomposition rate",
    "\u51cb\u843d\u7269\u8d28\u91cf\u635f\u5931\u7387" = "Litter mass loss rate",
    "\u51cb\u843d\u7269\u78b3\u6b8b\u7559\u7387" = "Litter carbon residual rate"
  )
  
  matched <- indicator_map[x]
  x[!is.na(matched)] <- matched[!is.na(matched)]
  
  return(x)
}

get_indicator_col <- function(df_raw) {
  candidate_cols <- c(
    "Trait", "trait",
    "Indicator", "indicator", "Indicators", "indicators",
    "Variable", "variable", "Variables", "variables",
    "Y_variable", "y_variable", "Y_variables", "y_variables",
    "Response", "response", "Response_variable", "response_variable",
    "Index", "index", "Item", "item",
    "Parameter", "parameter",
    "\u6307\u6807", "\u6307\u6807\u540d\u79f0",
    "\u54cd\u5e94\u6307\u6807",
    "\u53d8\u91cf", "\u53d8\u91cf\u540d\u79f0"
  )
  
  hit <- intersect(candidate_cols, names(df_raw))
  
  if (length(hit) > 0) {
    return(hit[1])
  } else {
    return(NULL)
  }
}

prepare_effect_size <- function(df_raw) {
  names_raw <- names(df_raw)
  
  yi_col <- find_col(
    names_raw,
    c("yi", "Yi", "effect_size", "Effect_Size", "EffectSize")
  )
  
  vi_col <- find_col(
    names_raw,
    c("vi", "Vi", "variance", "Variance", "var", "Var")
  )
  
  if (!is.null(yi_col) && !is.null(vi_col)) {
    yi <- to_num(df_raw[[yi_col]])
    vi <- to_num(df_raw[[vi_col]])
    
    return(list(
      yi = yi,
      vi = vi,
      source = "existing_yi_vi"
    ))
  }
  
  mean_t_col <- find_col(
    names_raw,
    c("mean_treatment", "Mean_treatment", "treatment_mean", "mean_treat", "m1i")
  )
  sd_t_col <- find_col(
    names_raw,
    c("sd_treatment", "SD_treatment", "treatment_sd", "sd_treat", "sd1i")
  )
  n_t_col <- find_col(
    names_raw,
    c("n_treatment", "N_treatment", "treatment_n", "n_treat", "n1i")
  )
  
  mean_c_col <- find_col(
    names_raw,
    c("mean_control", "Mean_control", "control_mean", "mean_ctrl", "m2i")
  )
  sd_c_col <- find_col(
    names_raw,
    c("sd_control", "SD_control", "control_sd", "sd_ctrl", "sd2i")
  )
  n_c_col <- find_col(
    names_raw,
    c("n_control", "N_control", "control_n", "n_ctrl", "n2i")
  )
  
  need_es <- c(mean_t_col, sd_t_col, n_t_col, mean_c_col, sd_c_col, n_c_col)
  
  if (any(sapply(need_es, is.null))) {
    stop("yi/vi columns were not found, and the mean/sd/n columns required for SMD calculation were not found. Please check the column names.")
  }
  
  esc <- escalc(
    measure = "SMD",
    m1i = to_num(df_raw[[mean_t_col]]),
    sd1i = to_num(df_raw[[sd_t_col]]),
    n1i = to_num(df_raw[[n_t_col]]),
    m2i = to_num(df_raw[[mean_c_col]]),
    sd2i = to_num(df_raw[[sd_c_col]]),
    n2i = to_num(df_raw[[n_c_col]])
  )
  
  return(list(
    yi = esc$yi,
    vi = esc$vi,
    source = "calculated_SMD"
  ))
}

save_fixed_panel_png <- function(p, filename, panel_size_cm = 5, dpi = 600,
                                 image_width_cm = 7.4, image_height_cm = 7.0) {
  g <- ggplotGrob(p)
  
  # Fix every layout column and row so that different tick-label widths
  # cannot move or resize the plotting frame.
  g$widths <- unit(rep(0, length(g$widths)), "cm")
  g$heights <- unit(rep(0, length(g$heights)), "cm")
  
  set_fixed_width <- function(layout_name, width_cm) {
    idx <- which(g$layout$name == layout_name)
    if (length(idx) > 0) {
      cols <- unique(unlist(Map(seq, g$layout$l[idx], g$layout$r[idx])))
      g$widths[cols] <<- unit(width_cm, "cm")
    }
  }
  
  set_fixed_height <- function(layout_name, height_cm) {
    idx <- which(g$layout$name == layout_name)
    if (length(idx) > 0) {
      rows <- unique(unlist(Map(seq, g$layout$t[idx], g$layout$b[idx])))
      g$heights[rows] <<- unit(height_cm, "cm")
    }
  }
  
  # Fixed horizontal structure: outer margin + y title + y ticks + panel + right space
  g$widths[1] <- unit(0.20, "cm")
  g$widths[length(g$widths)] <- unit(0.20, "cm")
  set_fixed_width("ylab-l", 0.65)
  set_fixed_width("axis-l", 0.85)
  set_fixed_width("panel", panel_size_cm)
  set_fixed_width("axis-r", 0.10)
  
  # Fixed vertical structure: outer margin + panel + x ticks + x title
  g$heights[1] <- unit(0.20, "cm")
  g$heights[length(g$heights)] <- unit(0.20, "cm")
  set_fixed_height("axis-t", 0.05)
  set_fixed_height("panel", panel_size_cm)
  set_fixed_height("axis-b", 0.50)
  set_fixed_height("xlab-b", 0.60)
  
  # Fixed full canvas size for every exported image
  png(
    filename = filename,
    width = image_width_cm,
    height = image_height_cm,
    units = "cm",
    res = dpi,
    bg = "white"
  )
  
  grid.newpage()
  grid.draw(g)
  dev.off()
}
# ==================== Funnel plot function ====================

plot_funnel <- function(data, filename, point_color = "gray50", point_border = "gray70",
                        xlab = "Effect Size (yi)", ylab = "Standard Error (SE)") {
  
  if (nrow(data) < 3) {
    warning(paste("Number of studies < 3. Plot skipped:", filename))
    return(list(model = NULL, p_value = NA, n = nrow(data)))
  }
  
  model <- tryCatch(
    rma(yi = yi, vi = vi, data = data, method = "REML"),
    error = function(e) NULL
  )
  
  if (is.null(model)) {
    warning(paste("Model fitting failed:", filename))
    return(list(model = NULL, p_value = NA, n = nrow(data)))
  }
  
  egger <- tryCatch(
    regtest(model, model = "lm"),
    error = function(e) NULL
  )
  
  p_val <- ifelse(is.null(egger), NA, egger$pval)
  
  plot_df <- data.frame(
    yi = data$yi,
    SE = sqrt(data$vi)
  )
  
  pooled_effect <- as.numeric(model$b[1])
  
  y_axis <- make_axis_3_ticks(plot_df$SE)
  
  se_seq <- seq(
    y_axis$limits[1],
    y_axis$limits[2],
    length.out = 100
  )
  
  upper_ci <- pooled_effect + 1.96 * se_seq
  lower_ci <- pooled_effect - 1.96 * se_seq
  
  ci_lines <- data.frame(
    SE = c(se_seq, se_seq),
    yi_bound = c(upper_ci, lower_ci),
    bound = rep(c("upper", "lower"), each = length(se_seq))
  )
  
  x_axis <- make_axis_3_ticks(c(plot_df$yi, upper_ci, lower_ci))
  
  p <- ggplot(plot_df, aes(x = yi, y = SE)) +
    geom_point(
      alpha = 0.8,
      size = 1.5,
      fill = point_color,
      color = point_border,
      shape = 21,
      stroke = 0.5
    ) +
    geom_vline(
      xintercept = pooled_effect,
      linetype = "dashed",
      color = "gray30",
      linewidth = 0.5
    ) +
    geom_line(
      data = ci_lines,
      aes(x = yi_bound, y = SE, group = bound),
      color = "gray30",
      linetype = "dashed",
      linewidth = 0.8
    ) +
    scale_y_reverse(
      limits = rev(y_axis$limits),
      breaks = y_axis$breaks,
      labels = axis_label_fun(y_axis$digits),
      expand = c(0, 0)
    ) +
    scale_x_continuous(
      limits = x_axis$limits,
      breaks = x_axis$breaks,
      labels = axis_label_fun(x_axis$digits),
      expand = c(0, 0)
    ) +
    labs(x = xlab, y = ylab) +
    coord_cartesian(clip = "off") +
    theme_minimal(base_family = font_family, base_size = font_size) +
    theme(
      text = element_text(family = font_family, face = "bold", size = font_size),
      axis.title = element_text(face = "bold", size = font_size),
      axis.text = element_text(face = "bold", size = font_size * 0.9),
      axis.ticks = element_line(color = "black", linewidth = 0.45),
      axis.ticks.length = unit(0.08, "cm"),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.45),
      aspect.ratio = 1,
      plot.margin = margin(0.50, 0.65, 0.50, 0.65, "cm")
    )
  
  save_fixed_panel_png(
    p = p,
    filename = filename,
    panel_size_cm = panel_size_cm,
    dpi = dpi
  )
  
  message("Saved: ", filename)
  
  return(list(
    model = model,
    p_value = p_val,
    n = nrow(data)
  ))
}

# ==================== Process one file ====================

process_file <- function(cfg) {
  message("\n========================================")
  message("Processing file: ", basename(cfg$path))
  
  if (!file.exists(cfg$path)) {
    warning("File does not exist: ", cfg$path)
    return(NULL)
  }
  
  df_raw <- as.data.frame(read_excel(cfg$path))
  names_raw <- names(df_raw)
  
  message("Original column names: ", paste(names_raw, collapse = ", "))
  
  es <- prepare_effect_size(df_raw)
  message("Effect size source: ", es$source)
  
  treatment_col <- find_col(
    names_raw,
    c(
      "Treatment_1", "Treatment", "treatment",
      "Group", "group",
      "Treatment type", "Treatment_type",
      "\u5904\u7406", "\u5904\u7406\u7c7b\u578b"
    )
  )
  
  if (is.null(treatment_col)) {
    warning("Treatment column was not found. File skipped: ", basename(cfg$path))
    return(NULL)
  }
  
  indicator_col <- get_indicator_col(df_raw)
  
  if (!is.null(indicator_col)) {
    message("Detected indicator column: ", indicator_col)
    indicator <- trimws(as.character(df_raw[[indicator_col]]))
    indicator[is.na(indicator) | indicator == ""] <- cfg$file_label
    indicator <- translate_indicator_name(indicator)
  } else {
    message("Indicator column was not detected. File label will be used as the indicator name: ", cfg$file_label)
    indicator <- rep(cfg$file_label, nrow(df_raw))
  }
  
  df <- data.frame(
    yi = to_num(es$yi),
    vi = to_num(es$vi),
    treatment = normalize_treatment(df_raw[[treatment_col]]),
    Indicator = indicator,
    stringsAsFactors = FALSE
  )
  
  df <- df %>%
    filter(
      !is.na(yi),
      !is.na(vi),
      vi > 0,
      !is.na(treatment),
      treatment %in% c(treatment_increased, treatment_decreased),
      !is.na(Indicator),
      Indicator != ""
    )
  
  if (nrow(df) == 0) {
    warning("No valid increased/decreased snowpack records were found: ", basename(cfg$path))
    return(NULL)
  }
  
  out_base <- file.path(root_output_dir, cfg$output_subdir, "Overall")
  
  dir.create(
    file.path(out_base, treatment_folder_names[treatment_increased]),
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  dir.create(
    file.path(out_base, treatment_folder_names[treatment_decreased]),
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  results_table <- data.frame(
    File = character(),
    Indicator = character(),
    Treatment = character(),
    N_studies = integer(),
    Egger_p = numeric(),
    Output_file = character(),
    stringsAsFactors = FALSE
  )
  
  indicators <- unique(df$Indicator)
  indicators <- indicators[!is.na(indicators)]
  
  treatments <- c(treatment_increased, treatment_decreased)
  
  for (ind in indicators) {
    for (trt in treatments) {
      
      sub_df <- df %>%
        filter(
          Indicator == ind,
          treatment == trt
        )
      
      safe_ind <- clean_filename(ind)
      safe_trt_folder <- treatment_folder_names[trt]
      
      sub_dir <- file.path(out_base, safe_trt_folder)
      filename <- file.path(sub_dir, paste0(safe_ind, ".png"))
      
      if (nrow(sub_df) < 3) {
        message(sprintf("Overall: %s - %s n=%d (<3), skipped", ind, trt, nrow(sub_df)))
        
        results_table <- rbind(results_table, data.frame(
          File = cfg$file_label,
          Indicator = ind,
          Treatment = trt,
          N_studies = nrow(sub_df),
          Egger_p = NA,
          Output_file = filename
        ))
        
        next
      }
      
      point_color <- treatment_colors[trt]
      
      res <- plot_funnel(
        data = sub_df,
        filename = filename,
        point_color = point_color,
        point_border = "gray70"
      )
      
      results_table <- rbind(results_table, data.frame(
        File = cfg$file_label,
        Indicator = ind,
        Treatment = trt,
        N_studies = res$n,
        Egger_p = res$p_value,
        Output_file = filename
      ))
    }
  }
  
  result_file <- file.path(
    file.path(root_output_dir, cfg$output_subdir),
    paste0(cfg$output_subdir, "_Overall_Egger_test_results.csv")
  )
  
  write.csv(results_table, result_file, row.names = FALSE, fileEncoding = "UTF-8")
  
  message("Egger test results saved to: ", result_file)
  
  return(results_table)
}

# ==================== Run all files ====================

dir.create(root_output_dir, recursive = TRUE, showWarnings = FALSE)

all_results <- list()

for (cfg in file_configs) {
  res <- process_file(cfg)
  if (!is.null(res)) {
    all_results[[cfg$file_label]] <- res
  }
}

if (length(all_results) > 0) {
  combined <- do.call(rbind, all_results)
  combined_file <- file.path(root_output_dir, "All_Overall_Egger_results.csv")
  
  write.csv(combined, combined_file, row.names = FALSE, fileEncoding = "UTF-8")
  
  message("\nAll overall funnel plots have been completed.")
  message("Main output directory: ", root_output_dir)
  message("Combined Egger test results saved to: ", combined_file)
}