# ============================================================================
# Table 10: Alternative Exposure Measures
# ============================================================================

# Load required libraries
library(haven)        # For reading STATA files
library(dplyr)        # For data manipulation
library(AER)          # For IV regression (ivreg function)
library(lmtest)       # For coeftest with clustering
library(sandwich)     # For clustered standard errors
library(modelsummary) # For creating regression tables

# Load the data
df <- read_dta('/Users/macca/Library/Mobile Documents/com~apple~CloudDocs/Uiversita /IV year/Empirical Econometrics/Final Paper/Public-Release-Data/dta/workfile_china.dta')

# ============================================================================
# Create list of region variables
# ============================================================================
reg_vars <- grep("^reg", names(df), value = TRUE)

# ============================================================================
# Section 1: Domestic plus International Exposure
# ============================================================================
# Panel A: Models 1-6 using d_tradex_usch_pw (imports + exports)

# Model 1: Mfg employment share (with first stage)
formula1_str <- paste(
  "d_sh_empl_mfg ~ d_tradex_usch_pw + l_shind_manuf_cbp +",
  paste(reg_vars, collapse = " + "),
  "+ l_sh_popedu_c + l_sh_popfborn + l_sh_empl_f",
  "+ l_sh_routine33 + l_task_outsource + t2 |",
  "d_tradex_otch_pw_lag + l_shind_manuf_cbp +",
  paste(reg_vars, collapse = " + "),
  "+ l_sh_popedu_c + l_sh_popfborn + l_sh_empl_f",
  "+ l_sh_routine33 + l_task_outsource + t2"
)

model1 <- ivreg(as.formula(formula1_str),
                data = df,
                weights = df$timepwt48)

# Models 2-6: Same instrument, different outcomes
outcomes_panel_a <- c(
  "d_sh_empl_nmfg",              # (2) Non-mfg employment share
  "d_avg_lnwkwage_mfg",          # (3) Mfg wages
  "d_avg_lnwkwage_nmfg",         # (4) Non-mfg wages
  "lnchg_trans_totindiv_pc",     # (5) Total transfers per capita
  "relchg_avg_hhincwage_pc_pw"   # (6) Household wage income
)

models_panel_a <- list(model1)

for (outcome in outcomes_panel_a) {
  formula_str <- paste(
    outcome, "~ d_tradex_usch_pw + l_shind_manuf_cbp +",
    paste(reg_vars, collapse = " + "),
    "+ l_sh_popedu_c + l_sh_popfborn + l_sh_empl_f",
    "+ l_sh_routine33 + l_task_outsource + t2 |",
    "d_tradex_otch_pw_lag + l_shind_manuf_cbp +",
    paste(reg_vars, collapse = " + "),
    "+ l_sh_popedu_c + l_sh_popfborn + l_sh_empl_f",
    "+ l_sh_routine33 + l_task_outsource + t2"
  )
  
  model <- ivreg(as.formula(formula_str),
                 data = df,
                 weights = df$timepwt48)
  
  models_panel_a[[length(models_panel_a) + 1]] <- model
}

# ============================================================================
# Section 2: Final Goods and Intermediate Imports
# ============================================================================
# Panel B: Models 7-12 using d_tradeusch_netinput_pw (net of intermediates)

outcomes_panel_b <- c(
  "d_sh_empl_mfg",               # (7) Mfg employment share
  "d_sh_empl_nmfg",              # (8) Non-mfg employment share
  "d_avg_lnwkwage_mfg",          # (9) Mfg wages
  "d_avg_lnwkwage_nmfg",         # (10) Non-mfg wages
  "lnchg_trans_totindiv_pc",     # (11) Total transfers per capita
  "relchg_avg_hhincwage_pc_pw"   # (12) Household wage income
)

models_panel_b <- list()

for (outcome in outcomes_panel_b) {
  formula_str <- paste(
    outcome, "~ d_tradeusch_netinput_pw + l_shind_manuf_cbp +",
    paste(reg_vars, collapse = " + "),
    "+ l_sh_popedu_c + l_sh_popfborn + l_sh_empl_f",
    "+ l_sh_routine33 + l_task_outsource + t2 |",
    "d_tradeotch_pw_lag + d_inputotch_pw_lag + l_shind_manuf_cbp +",
    paste(reg_vars, collapse = " + "),
    "+ l_sh_popedu_c + l_sh_popfborn + l_sh_empl_f",
    "+ l_sh_routine33 + l_task_outsource + t2"
  )
  
  model <- ivreg(as.formula(formula_str),
                 data = df,
                 weights = df$timepwt48)
  
  models_panel_b[[length(models_panel_b) + 1]] <- model
}

# ============================================================================
# Section 3: Net Imports
# ============================================================================
# Panel C: Models 13-18 using d_netimpusch_pw (imports - exports)

outcomes_panel_c <- c(
  "d_sh_empl_mfg",               # (13) Mfg employment share
  "d_sh_empl_nmfg",              # (14) Non-mfg employment share
  "d_avg_lnwkwage_mfg",          # (15) Mfg wages
  "d_avg_lnwkwage_nmfg",         # (16) Non-mfg wages
  "lnchg_trans_totindiv_pc",     # (17) Total transfers per capita
  "relchg_avg_hhincwage_pc_pw"   # (18) Household wage income
)

models_panel_c <- list()

for (outcome in outcomes_panel_c) {
  formula_str <- paste(
    outcome, "~ d_netimpusch_pw + l_shind_manuf_cbp +",
    paste(reg_vars, collapse = " + "),
    "+ l_sh_popedu_c + l_sh_popfborn + l_sh_empl_f",
    "+ l_sh_routine33 + l_task_outsource + t2 |",
    "d_tradeotch_pw_lag + d_expotch_pw_lag + l_shind_manuf_cbp +",
    paste(reg_vars, collapse = " + "),
    "+ l_sh_popedu_c + l_sh_popfborn + l_sh_empl_f",
    "+ l_sh_routine33 + l_task_outsource + t2"
  )
  
  model <- ivreg(as.formula(formula_str),
                 data = df,
                 weights = df$timepwt48)
  
  models_panel_c[[length(models_panel_c) + 1]] <- model
}

# ============================================================================
# Section 4: Gravity Residual (Reduced Form OLS)
# ============================================================================
# Panel D: Models 19-24 using d_traderes_pw_lag (OLS, not IV)

outcomes_panel_d <- c(
  "d_sh_empl_mfg",               # (19) Mfg employment share
  "d_sh_empl_nmfg",              # (20) Non-mfg employment share
  "d_avg_lnwkwage_mfg",          # (21) Mfg wages
  "d_avg_lnwkwage_nmfg",         # (22) Non-mfg wages
  "lnchg_trans_totindiv_pc",     # (23) Total transfers per capita
  "relchg_avg_hhincwage_pc_pw"   # (24) Household wage income
)

models_panel_d <- list()

for (outcome in outcomes_panel_d) {
  formula_str <- paste(
    outcome, "~ d_traderes_pw_lag + l_shind_manuf_cbp +",
    paste(reg_vars, collapse = " + "),
    "+ l_sh_popedu_c + l_sh_popfborn + l_sh_empl_f",
    "+ l_sh_routine33 + l_task_outsource + t2"
  )
  
  model <- lm(as.formula(formula_str),
              data = df,
              weights = df$timepwt48)
  
  models_panel_d[[length(models_panel_d) + 1]] <- model
}

# ============================================================================
# Section 5: Factor Content of Net Imports
# ============================================================================
# Panel E: Models 25-30 using d_nettradefactor_usch_io

outcomes_panel_e <- c(
  "d_sh_empl_mfg",               # (25) Mfg employment share
  "d_sh_empl_nmfg",              # (26) Non-mfg employment share
  "d_avg_lnwkwage_mfg",          # (27) Mfg wages
  "d_avg_lnwkwage_nmfg",         # (28) Non-mfg wages
  "lnchg_trans_totindiv_pc",     # (29) Total transfers per capita
  "relchg_avg_hhincwage_pc_pw"   # (30) Household wage income
)

models_panel_e <- list()

for (outcome in outcomes_panel_e) {
  formula_str <- paste(
    outcome, "~ d_nettradefactor_usch_io + l_shind_manuf_cbp +",
    paste(reg_vars, collapse = " + "),
    "+ l_sh_popedu_c + l_sh_popfborn + l_sh_empl_f",
    "+ l_sh_routine33 + l_task_outsource + t2 |",
    "d_tradefactor_otch_lag_io + d_expfactor_otch_lag_io + l_shind_manuf_cbp +",
    paste(reg_vars, collapse = " + "),
    "+ l_sh_popedu_c + l_sh_popfborn + l_sh_empl_f",
    "+ l_sh_routine33 + l_task_outsource + t2"
  )
  
  model <- ivreg(as.formula(formula_str),
                 data = df,
                 weights = df$timepwt48)
  
  models_panel_e[[length(models_panel_e) + 1]] <- model
}

# ============================================================================
# Combine all models
# ============================================================================
all_models <- c(models_panel_a, models_panel_b, models_panel_c, 
                models_panel_d, models_panel_e)

# ============================================================================
# Create regression table with descriptive names
# ============================================================================

# Custom coefficient names (row names) - different for each panel
coef_names_iv <- c(
  "d_tradex_usch_pw" = "Δ Import + Export Exposure (per worker)",
  "d_tradeusch_netinput_pw" = "Δ Import Exposure, Net of Intermediates (per worker)",
  "d_netimpusch_pw" = "Δ Net Import Exposure (per worker)",
  "d_traderes_pw_lag" = "Gravity Residual Import Exposure (per worker)",
  "d_nettradefactor_usch_io" = "Δ Factor Content of Net Imports (per worker)",
  "l_shind_manuf_cbp" = "Manufacturing Share (lag)",
  "l_sh_popedu_c" = "College Educated Share (lag)",
  "l_sh_popfborn" = "Foreign Born Share (lag)",
  "l_sh_empl_f" = "Female Employment Share (lag)",
  "l_sh_routine33" = "Routine Occupation Share (lag)",
  "l_task_outsource" = "Offshorability Index (lag)"
)

# Descriptive column names
col_names <- c(
  # Panel A: Domestic + International (1-6)
  "(1)\nΔ Mfg Empl/Pop\n[Import+Export]",
  "(2)\nΔ Non-Mfg Empl/Pop\n[Import+Export]",
  "(3)\nΔ Log Wage Mfg\n[Import+Export]",
  "(4)\nΔ Log Wage Non-Mfg\n[Import+Export]",
  "(5)\nLog Δ Transfers p.c.\n[Import+Export]",
  "(6)\n% Δ HH Wage Inc\n[Import+Export]",
  
  # Panel B: Net of Intermediates (7-12)
  "(7)\nΔ Mfg Empl/Pop\n[Net Interm]",
  "(8)\nΔ Non-Mfg Empl/Pop\n[Net Interm]",
  "(9)\nΔ Log Wage Mfg\n[Net Interm]",
  "(10)\nΔ Log Wage Non-Mfg\n[Net Interm]",
  "(11)\nLog Δ Transfers p.c.\n[Net Interm]",
  "(12)\n% Δ HH Wage Inc\n[Net Interm]",
  
  # Panel C: Net Imports (13-18)
  "(13)\nΔ Mfg Empl/Pop\n[Net Imports]",
  "(14)\nΔ Non-Mfg Empl/Pop\n[Net Imports]",
  "(15)\nΔ Log Wage Mfg\n[Net Imports]",
  "(16)\nΔ Log Wage Non-Mfg\n[Net Imports]",
  "(17)\nLog Δ Transfers p.c.\n[Net Imports]",
  "(18)\n% Δ HH Wage Inc\n[Net Imports]",
  
  # Panel D: Gravity Residual (19-24)
  "(19)\nΔ Mfg Empl/Pop\n[Gravity OLS]",
  "(20)\nΔ Non-Mfg Empl/Pop\n[Gravity OLS]",
  "(21)\nΔ Log Wage Mfg\n[Gravity OLS]",
  "(22)\nΔ Log Wage Non-Mfg\n[Gravity OLS]",
  "(23)\nLog Δ Transfers p.c.\n[Gravity OLS]",
  "(24)\n% Δ HH Wage Inc\n[Gravity OLS]",
  
  # Panel E: Factor Content (25-30)
  "(25)\nΔ Mfg Empl/Pop\n[Factor Content]",
  "(26)\nΔ Non-Mfg Empl/Pop\n[Factor Content]",
  "(27)\nΔ Log Wage Mfg\n[Factor Content]",
  "(28)\nΔ Log Wage Non-Mfg\n[Factor Content]",
  "(29)\nLog Δ Transfers p.c.\n[Factor Content]",
  "(30)\n% Δ HH Wage Inc\n[Factor Content]"
)

# Create table with clustered standard errors
table10a <- modelsummary(
  setNames(models_panel_a, col_names[1:6]),
  estimate = "{estimate}",
  statistic = "({std.error})",
  vcov = lapply(1:6, function(x) ~statefip),
  coef_omit = "^reg|^t2|^t1|Intercept",
  coef_map = coef_names_iv,
  gof_map = c("nobs", "r.squared"),
  output = "gt",
  fmt = 3,
  title = "Table 10A: Import + Export Exposure"
)
gt::gtsave(table10a,
           filename = paste0(normalizePath("~/Desktop"), "/Table 10a.png"),
           vwidth   = 6000,
           vheight  = 1600)

table10b <- modelsummary(
  setNames(models_panel_b, col_names[7:12]),
  estimate = "{estimate}",
  statistic = "({std.error})",
  vcov = lapply(1:6, function(x) ~statefip),
  coef_omit = "^reg|^t2|^t1|Intercept",
  coef_map = coef_names_iv,
  gof_map = c("nobs", "r.squared"),
  output = "gt",
  fmt = 3,
  title = "Table 10B: Net of Intermediate Imports"
)
gt::gtsave(table10b,
           filename = paste0(normalizePath("~/Desktop"), "/Table 10b.png"),
           vwidth   = 6000,
           vheight  = 1600)

table10c <- modelsummary(
  setNames(models_panel_c, col_names[13:18]),
  estimate = "{estimate}",
  statistic = "({std.error})",
  vcov = lapply(1:6, function(x) ~statefip),
  coef_omit = "^reg|^t2|^t1|Intercept",
  coef_map = coef_names_iv,
  gof_map = c("nobs", "r.squared"),
  output = "gt",
  fmt = 3,
  title = "Table 10C: Net Imports (Imports - Exports)"
)
gt::gtsave(table10c,
           filename = paste0(normalizePath("~/Desktop"), "/Table 10c.png"),
           vwidth   = 6000,
           vheight  = 1600)

table10d <- modelsummary(
  setNames(models_panel_d, col_names[19:24]),
  estimate = "{estimate}",
  statistic = "({std.error})",
  vcov = lapply(1:6, function(x) ~statefip),
  coef_omit = "^reg|^t2|^t1|Intercept",
  coef_map = coef_names_iv,
  gof_map = c("nobs", "r.squared"),
  output = "gt",
  fmt = 3,
  title = "Table 10D: Gravity Residual (OLS, not IV)"
)
gt::gtsave(table10d,
           filename = paste0(normalizePath("~/Desktop"), "/Table 10d.png"),
           vwidth   = 6000,
           vheight  = 1600)

table10e <- modelsummary(
  setNames(models_panel_e, col_names[25:30]),
  estimate = "{estimate}",
  statistic = "({std.error})",
  vcov = lapply(1:6, function(x) ~statefip),
  coef_omit = "^reg|^t2|^t1|Intercept",
  coef_map = coef_names_iv,
  gof_map = c("nobs", "r.squared"),
  output = "gt",
  fmt = 3,
  title = "Table 10E: Factor Content of Net Imports"
)
gt::gtsave(table10e,
           filename = paste0(normalizePath("~/Desktop"), "/Table 10e.png"),
           vwidth   = 6000,
           vheight  = 1600)
# ============================================================================
# Create summary by panel
# ============================================================================
create_panel_summary <- function(models, panel_name, exposure_var) {
  data.frame(
    panel = panel_name,
    outcome = c("Mfg Empl/Pop", "Non-Mfg Empl/Pop", "Mfg Wages", 
                "Non-Mfg Wages", "Transfers", "HH Wage Inc"),
    coef = sapply(models, function(m) coef(m)[exposure_var]),
    se = sapply(models, function(m) sqrt(diag(vcovCL(m, cluster = df$statefip)))[exposure_var]),
    n = sapply(models, nobs),
    r2 = sapply(models, function(m) summary(m)$r.squared)
  )
}

summary_a <- create_panel_summary(models_panel_a, "A: Import+Export", "d_tradex_usch_pw")
summary_b <- create_panel_summary(models_panel_b, "B: Net Intermediates", "d_tradeusch_netinput_pw")
summary_c <- create_panel_summary(models_panel_c, "C: Net Imports", "d_netimpusch_pw")
summary_d <- create_panel_summary(models_panel_d, "D: Gravity (OLS)", "d_traderes_pw_lag")
summary_e <- create_panel_summary(models_panel_e, "E: Factor Content", "d_nettradefactor_usch_io")

all_summaries <- rbind(summary_a, summary_b, summary_c, summary_d, summary_e)
all_summaries$t_stat <- all_summaries$coef / all_summaries$se
all_summaries$p_value <- 2 * pt(-abs(all_summaries$t_stat), df = all_summaries$n - 1)

print(all_summaries)

# ============================================================================
# Notes:
# ============================================================================
# Table 10 tests robustness using 5 alternative trade exposure measures:
# 
# Panel A (Models 1-6): Import + Export exposure
#   - Exposure: d_tradex_usch_pw (imports + exports per worker)
#   - Instrument: d_tradex_otch_pw_lag (other countries' imp+exp from China)
#
# Panel B (Models 7-12): Net of intermediate inputs
#   - Exposure: d_tradeusch_netinput_pw (imports minus intermediates)
#   - Instruments: d_tradeotch_pw_lag + d_inputotch_pw_lag
#
# Panel C (Models 13-18): Net imports (imports - exports)
#   - Exposure: d_netimpusch_pw (imports minus exports per worker)
#   - Instruments: d_tradeotch_pw_lag + d_expotch_pw_lag
#
# Panel D (Models 19-24): Gravity residual (OLS, not IV)
#   - Exposure: d_traderes_pw_lag (gravity model residuals)
#   - Method: OLS regression (reduced form)
#
# Panel E (Models 25-30): Factor content of trade
#   - Exposure: d_nettradefactor_usch_io (labor content of net imports)
#   - Instruments: d_tradefactor_otch_lag_io + d_expfactor_otch_lag_io
#
# All models include full controls and cluster SE at state level