# Effect Size Calculation and Meta-Analysis Using a Variance-Covariance Matrix

# Build dataset
## Overall dataset
### Read file
library(readxl)
require(ggmap)
require(maps)
library(metafor)
library(ggplot2)
library(glmulti)
library(writexl)
library(Matrix)


# 1. Read data
d1 <- read_excel("C:/Users/Administrator/Desktop/25meta/Meta date 2025.12.2.xlsx")

# 2. Calculate effect sizes
d2 <- escalc(measure = "SMD", 
             data = d1,
             m1i = mean_treatment, 
             sd1i = sd_treatment, 
             n1i = n_treatment,
             m2i = mean_control, 
             sd2i = sd_control, 
             n2i = n_control)

# 3. Data cleaning: remove problematic data
d2_clean <- d2[!is.na(d2$yi) & !is.na(d2$vi) & 
                 d2$vi > 0 & 
                 d2$n_treatment > 0 & 
                 d2$n_control > 0, ]

cat("Number of rows in original data:", nrow(d2), "\n")
cat("Number of rows after cleaning:", nrow(d2_clean), "\n")

# 4. Use the original cal.v function with slight improvements
cal.v <- function(x){
  k <- nrow(x)
  
  if (k == 1) {
    # If there is only one row, directly return the variance
    return(matrix(x$vi, nrow = 1, ncol = 1))
  } else {
    # Calculate the base covariance value
    base_cov <- mean(x$sd_control^2 / 
                       (x$n_control * pmax(x$mean_control^2, 0.01)))
    
    # Construct the matrix
    v <- matrix(base_cov, nrow = k, ncol = k)
    diag(v) <- x$vi
    return(v)
  }
}

# 5. Construct the V matrix
V <- bldiag(lapply(split(d2_clean, d2_clean$Common_id), cal.v))

# 6. Check positive definiteness of the matrix
eigen_vals <- eigen(as.matrix(V), symmetric = TRUE, only.values = TRUE)$values
cat("Minimum eigenvalue:", min(eigen_vals), "\n")

# 7. If the matrix is not positive definite, use nearPD to repair it
if (min(eigen_vals) <= 1e-10) {
  cat("The matrix is not positive definite. Repairing...\n")
  
  V_pd <- nearPD(V, keepDiag = TRUE, conv.tol = 1e-7)
  V <- as.matrix(V_pd$mat)
  
  cat("Repair completed. New minimum eigenvalue:", 
      min(eigen(as.matrix(V), symmetric = TRUE, only.values = TRUE)$values), 
      "\n")
}

# Use V instead of vi for meta-analysis
m0 <- rma.mv(yi, V, 
             data = d2_clean, 
             random = list(~1 | Study_id / Row_id), 
             method = "REML")

summary(m0)


# Combine yi and the main diagonal of V
d2_clean$V <- diag(V)

# Convert the full variance-covariance matrix to a data frame in long format
V_df <- as.data.frame(as.table(V))
names(V_df) <- c("row", "col", "covariance")

# Export to an Excel file
write_xlsx(d2_clean, "C:/Users/Administrator/Desktop/25meta/Meta effect size r2025.12.2.xlsx")