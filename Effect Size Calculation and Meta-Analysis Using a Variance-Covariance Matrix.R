# Effect Size Calculation and Meta-Analysis Using a Variance-Covariance Matrix

# Build dataset
## Overall dataset
### Read file

library(readxl)
library(metafor)
library(writexl)
library(Matrix)

# 1. Read data
d1 <- read_excel("C:/Users/Administrator/Desktop/25meta/Meta date 2025.12.2.xlsx")

# 2. Calculate effect sizes
d2 <- escalc(
  measure = "SMD", 
  data = d1,
  m1i = mean_treatment, 
  sd1i = sd_treatment, 
  n1i = n_treatment,
  m2i = mean_control, 
  sd2i = sd_control, 
  n2i = n_control
)

# 3. Data cleaning: remove problematic data
d2_clean <- d2[
  !is.na(d2$yi) & 
    !is.na(d2$vi) & 
    d2$vi > 0 & 
    d2$n_treatment > 0 & 
    d2$n_control > 0, 
]

# 4. Check required columns
required_cols <- c("Common_id", "Study_id", "Row_id")

for(col in required_cols) {
  if(!col %in% names(d2_clean)) {
    stop(paste0("Required column '", col, "' was not found in the input data."))
  }
}

# 5. Sort data to ensure that the order of yi matches the order of the V matrix
d2_clean <- d2_clean[order(d2_clean$Common_id, d2_clean$Row_id), ]

cat("Number of rows in original data:", nrow(d2), "\n")
cat("Number of rows after cleaning:", nrow(d2_clean), "\n")

# 6. Define the variance-covariance matrix function
cal.v <- function(x) {
  k <- nrow(x)
  
  if(k == 1) {
    # If there is only one row, directly return the variance
    return(matrix(x$vi, nrow = 1, ncol = 1))
  } else {
    # Calculate the base covariance value
    base_cov <- mean(
      x$sd_control^2 / 
        (x$n_control * pmax(x$mean_control^2, 0.01)),
      na.rm = TRUE
    )
    
    # Avoid NA or infinite covariance values
    if(is.na(base_cov) || !is.finite(base_cov)) {
      base_cov <- 0
    }
    
    # Construct the matrix
    v <- matrix(base_cov, nrow = k, ncol = k)
    diag(v) <- x$vi
    
    return(v)
  }
}

# 7. Construct the V matrix
V <- bldiag(lapply(split(d2_clean, d2_clean$Common_id), cal.v))
V <- as.matrix(V)

# 8. Check whether the dimension of V matches the number of rows in d2_clean
if(nrow(V) != nrow(d2_clean) || ncol(V) != nrow(d2_clean)) {
  stop("The dimension of the V matrix does not match the number of rows in d2_clean.")
}

# 9. Check positive definiteness of the matrix
eigen_vals <- eigen(V, symmetric = TRUE, only.values = TRUE)$values
cat("Minimum eigenvalue:", min(eigen_vals), "\n")

# 10. If the matrix is not positive definite, use nearPD to repair it
if(min(eigen_vals) <= 1e-10) {
  cat("The matrix is not positive definite. Repairing...\n")
  
  V_pd <- nearPD(V, keepDiag = TRUE, conv.tol = 1e-7)
  V <- as.matrix(V_pd$mat)
  
  new_eigen_vals <- eigen(V, symmetric = TRUE, only.values = TRUE)$values
  
  cat("Repair completed. New minimum eigenvalue:", min(new_eigen_vals), "\n")
}

# 11. Use V instead of vi for meta-analysis
m0 <- rma.mv(
  yi,
  V = V,
  data = d2_clean,
  random = list(~ 1 | Study_id / Row_id),
  method = "REML"
)

summary(m0)

# 12. Combine yi and the main diagonal of V
d2_clean$V <- diag(V)

# 13. Convert the full variance-covariance matrix to a data frame in long format
V_df <- as.data.frame(as.table(V))
names(V_df) <- c("row", "col", "covariance")

# 14. Export to Excel files
write_xlsx(
  d2_clean,
  "C:/Users/Administrator/Desktop/25meta/Meta_effect_sizes_r_2025.12.2_version_0-1.xlsx"
)

write_xlsx(
  V_df,
  "C:/Users/Administrator/Desktop/25meta/Variance_covariance_matrix_long_format_2025.12.2.xlsx"
)

cat("Effect-size dataset exported successfully.\n")
cat("Variance-covariance matrix exported successfully.\n")