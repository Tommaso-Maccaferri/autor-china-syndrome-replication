# ============================================================================
# Extended Table 3: OLS, IV (lagged & non-lagged), FE, FD models, ER SE
# ============================================================================

library(haven)
library(dplyr)
library(AER)
library(sandwich)
library(fixest)
library(gt)
library(plm)

# Load data
df <- read_dta("/Users/macca/Library/Mobile Documents/com~apple~CloudDocs/Uiversita /IV year/Empirical Econometrics/Final Paper/Public-Release-Data/dta/workfile_china.dta")
sic87dd <- read_dta('/Users/macca/Library/Mobile Documents/com~apple~CloudDocs/Uiversita /IV year/Empirical Econometrics/Final Paper/Public-Release-Data/dta/sic87dd_trade_data.dta')

# ============================================================================
# Model specifications
# ============================================================================
reg_vars <- grep("^reg", names(df), value = TRUE)

# OLS baseline
ols_model <- lm(d_sh_empl_mfg ~ d_tradeusch_pw + t2,
                data = df, weights = df$timepwt48)

# IV models (lagged instrument)
model1 <- ivreg(d_sh_empl_mfg ~ d_tradeusch_pw + t2 | d_tradeotch_pw_lag + t2,
                data = df, weights = df$timepwt48)

model2 <- ivreg(d_sh_empl_mfg ~ d_tradeusch_pw + l_shind_manuf_cbp + t2 | 
                  d_tradeotch_pw_lag + l_shind_manuf_cbp + t2,
                data = df, weights = df$timepwt48)

formula3_str <- paste("d_sh_empl_mfg ~ d_tradeusch_pw + l_shind_manuf_cbp +", 
                      paste(reg_vars, collapse = " + "), "+ t2 |",
                      "d_tradeotch_pw_lag + l_shind_manuf_cbp +", 
                      paste(reg_vars, collapse = " + "), "+ t2")
model3 <- ivreg(as.formula(formula3_str), data = df, weights = df$timepwt48)

formula4_str <- paste("d_sh_empl_mfg ~ d_tradeusch_pw + l_shind_manuf_cbp +", 
                      paste(reg_vars, collapse = " + "), 
                      "+ l_sh_popedu_c + l_sh_popfborn + l_sh_empl_f + t2 |",
                      "d_tradeotch_pw_lag + l_shind_manuf_cbp +", 
                      paste(reg_vars, collapse = " + "), 
                      "+ l_sh_popedu_c + l_sh_popfborn + l_sh_empl_f + t2")
model4 <- ivreg(as.formula(formula4_str), data = df, weights = df$timepwt48)

formula5_str <- paste("d_sh_empl_mfg ~ d_tradeusch_pw + l_shind_manuf_cbp +", 
                      paste(reg_vars, collapse = " + "), 
                      "+ l_sh_routine33 + l_task_outsource + t2 |",
                      "d_tradeotch_pw_lag + l_shind_manuf_cbp +", 
                      paste(reg_vars, collapse = " + "), 
                      "+ l_sh_routine33 + l_task_outsource + t2")
model5 <- ivreg(as.formula(formula5_str), data = df, weights = df$timepwt48)

formula6_str <- paste("d_sh_empl_mfg ~ d_tradeusch_pw + l_shind_manuf_cbp +", 
                      paste(reg_vars, collapse = " + "), 
                      "+ l_sh_popedu_c + l_sh_popfborn + l_sh_empl_f",
                      "+ l_sh_routine33 + l_task_outsource + t2 |",
                      "d_tradeotch_pw_lag + l_shind_manuf_cbp +", 
                      paste(reg_vars, collapse = " + "), 
                      "+ l_sh_popedu_c + l_sh_popfborn + l_sh_empl_f",
                      "+ l_sh_routine33 + l_task_outsource + t2")
model6 <- ivreg(as.formula(formula6_str), data = df, weights = df$timepwt48)

# IV model (non-lagged instrument)
formula6_nolag_str <- paste(
  "d_sh_empl_mfg ~ d_tradeusch_pw + l_shind_manuf_cbp +",
  paste(reg_vars, collapse = " + "),
  "+ l_sh_popedu_c + l_sh_popfborn + l_sh_empl_f",
  "+ l_sh_routine33 + l_task_outsource + t2 |",
  "l_tradeotch_pw + l_shind_manuf_cbp +",
  paste(reg_vars, collapse = " + "),
  "+ l_sh_popedu_c + l_sh_popfborn + l_sh_empl_f",
  "+ l_sh_routine33 + l_task_outsource + t2"
)
model6_nolag <- ivreg(as.formula(formula6_nolag_str), data = df, weights = df$timepwt48)

# Fixed-effects model
fe_model <- feols(
  d_sh_empl_mfg ~ d_tradeusch_pw + l_shind_manuf_cbp +
    l_sh_popedu_c + l_sh_popfborn + l_sh_empl_f +
    l_sh_routine33 + l_task_outsource | czone + t2,
  data = df, weights = ~timepwt48, cluster = ~statefip
)

# First-difference model
fd_model <- feols(
  d_sh_empl_mfg ~ d_tradeusch_pw + l_shind_manuf_cbp +
    l_sh_popedu_c + l_sh_popfborn + l_sh_empl_f +
    l_sh_routine33 + l_task_outsource | t2,
  data = df, weights = ~timepwt48, cluster = ~czone
)

# ============================================================================
# ============================================================================
# Exposure-robust SE for model6
# ============================================================================
exposure_var <- df$timepwt48
X <- model.matrix(model6)
residuals <- residuals(model6)
n_units <- length(unique(df$czone))
vcov_exposure_robust <- vcov(model6) * 0
for (i in unique(df$czone)) {
  idx <- which(df$czone == i)
  X_i <- X[idx, , drop = FALSE]
  e_i <- residuals[idx]
  exp_i <- exposure_var[idx]
  vcov_exposure_robust <- vcov_exposure_robust + 
    t(X_i) %*% diag(exp_i^2 * e_i^2) %*% X_i
}
vcov_exposure_robust <- solve(t(X) %*% diag(exposure_var) %*% X) %*% 
  vcov_exposure_robust %*% 
  solve(t(X) %*% diag(exposure_var) %*% X)

# ============================================================================
# Extract results
# ============================================================================
get_result <- function(model, coef_name, cluster_var = NULL, vcov_custom = NULL) {
  if (!coef_name %in% names(coef(model))) return(NULL)
  
  beta <- coef(model)[coef_name]
  
  if (!is.null(vcov_custom)) {
    se <- sqrt(vcov_custom[coef_name, coef_name])
  } else if (inherits(model, "fixest")) {
    se <- sqrt(vcov(model)[coef_name, coef_name])
  } else {
    vcov_cl <- vcovCL(model, cluster = if(is.null(cluster_var)) df$statefip else cluster_var)
    se <- sqrt(vcov_cl[coef_name, coef_name])
  }
  
  t_stat <- beta / se
  df_resid <- nobs(model) - length(coef(model))
  p_val <- 2 * pt(-abs(t_stat), df = df_resid)
  
  stars <- if (p_val < 0.01) "***" else if (p_val < 0.05) "**" else if (p_val < 0.1) "*" else ""
  
  list(coef = beta, se = se, stars = stars, n = nobs(model), 
       r2 = if(inherits(model, "fixest")) r2(model)["r2"] else summary(model)$r.squared)
}

# Get first-stage F-statistics
get_fstat <- function(model) {
  if (inherits(model, "ivreg")) {
    fs <- summary(model, diagnostics = TRUE)
    return(fs$diagnostics["Weak instruments", "statistic"])
  }
  return(NA)
}

# Build results list
models_list <- list(ols_model, model1, model2, model3, model4, 
                    model5, model6, model6_nolag, fe_model, fd_model)
col_labels <- c("OLS", "(1)", "(2)", "(3)", "(4)", "(5)", "(6)", 
                "No Lag", "FE model", "CZ clusters")

all_results <- lapply(models_list, function(m) {
  result <- get_result(m, "d_tradeusch_pw")
  result$fstat <- get_fstat(m)
  result
})

# Add ER SE result
result_er <- get_result(model6, "d_tradeusch_pw", vcov_custom = vcov_exposure_robust)
result_er$n <- nobs(model6)
result_er$r2 <- summary(model6)$r.squared
result_er$fstat <- get_fstat(model6)
all_results <- c(all_results, list(result_er))
col_labels <- c(col_labels, "ER SE")

# ============================================================================
# Build table
# ============================================================================
build_var_row <- function(var_name, var_label) {
  coef_row <- tibble(term = var_label)
  se_row <- tibble(term = "")
  
  for (i in 1:10) {
    result <- get_result(models_list[[i]], var_name)
    if (!is.null(result)) {
      coef_row[[col_labels[i]]] <- sprintf("%.3f%s", result$coef, result$stars)
      se_row[[col_labels[i]]] <- sprintf("(%.3f)", result$se)
    } else {
      coef_row[[col_labels[i]]] <- ""
      se_row[[col_labels[i]]] <- ""
    }
  }
  
  # ER SE column
  result_er_var <- get_result(model6, var_name, vcov_custom = vcov_exposure_robust)
  if (!is.null(result_er_var)) {
    coef_row[["ER SE"]] <- sprintf("%.3f%s", result_er_var$coef, result_er_var$stars)
    se_row[["ER SE"]] <- sprintf("(%.3f)", result_er_var$se)
  } else {
    coef_row[["ER SE"]] <- ""
    se_row[["ER SE"]] <- ""
  }
  
  bind_rows(coef_row, se_row)
}

table_data <- bind_rows(
  build_var_row("d_tradeusch_pw", "Δ Chinese Import Exposure (per worker)"),
  build_var_row("l_shind_manuf_cbp", "Manufacturing Share (lag)"),
  build_var_row("l_sh_popedu_c", "College Educated Share (lag)"),
  build_var_row("l_sh_popfborn", "Foreign Born Share (lag)"),
  build_var_row("l_sh_empl_f", "Female Employment Share (lag)"),
  build_var_row("l_sh_routine33", "Routine Occupation Share (lag)"),
  build_var_row("l_task_outsource", "Offshorability Index (lag)")
)

# Add goodness-of-fit rows
n_row <- tibble(term = "N")
r2_row <- tibble(term = "R²")
fstat_row <- tibble(term = "First-stage F-statistic")

for (i in 1:length(all_results)) {
  n_row[[col_labels[i]]] <- as.character(all_results[[i]]$n)
  r2_row[[col_labels[i]]] <- sprintf("%.3f", all_results[[i]]$r2)
  fstat_row[[col_labels[i]]] <- if (is.na(all_results[[i]]$fstat)) "/" else sprintf("%.1f", all_results[[i]]$fstat)
}

table_data <- bind_rows(table_data, n_row, r2_row, fstat_row)

# ============================================================================
# Create GT table
# ============================================================================
table3 <- table_data %>%
  gt() %>%
  tab_header(
    title = "Extended Table 3: Change in Manufacturing Employment Share (Δ Mfg Empl/Pop)",
    subtitle = "Dependent Variable: Ten-year Equivalent Change in Manufacturing Employment Share"
  ) %>%
  cols_align(align = "left") %>%
  cols_align(align = "center") %>%
  tab_options(
    table.font.size = px(11),
    heading.title.font.size = px(14),
    heading.subtitle.font.size = px(12),
    table_body.hlines.color = "transparent",
    table_body.border.bottom.color = "black",
    table_body.border.top.color = "black",
    column_labels.border.bottom.color = "black"
  ) %>%
  tab_footnote(
    footnote = "Standard errors in parentheses, clustered at state level (except CZ clusters column, clustered at CZ level; and ER SE column using exposure-robust formula)."
  ) %>%
  tab_footnote(
    footnote = "* p<0.1, ** p<0.05, *** p<0.01"
  )

print(table3)

gt::gtsave(table3,
           filename = paste0(normalizePath("~/Desktop"), "/Table 3.png"),
           vwidth   = 6000,
           vheight  = 1600)