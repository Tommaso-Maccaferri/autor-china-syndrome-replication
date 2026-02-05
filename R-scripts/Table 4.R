# ============================================================================
# Table 4: Population Change
# ============================================================================

library(haven)
library(dplyr)
library(AER)
library(sandwich)
library(gt)
library(tibble)

# Load data
df <- read_dta('/Users/macca/Library/Mobile Documents/com~apple~CloudDocs/Uiversita /IV year/Empirical Econometrics/Final Paper/Public-Release-Data/dta/workfile_china.dta')

reg_vars <- grep("^reg", names(df), value = TRUE)

outcomes <- c(
  "lnchg_popworkage",
  "lnchg_popworkage_edu_c",
  "lnchg_popworkage_edu_nc",
  "lnchg_popworkage_age1634",
  "lnchg_popworkage_age3549",
  "lnchg_popworkage_age5064"
)

outcome_labels <- c("All", "College", "Non-College", "Age 16-34", "Age 35-49", "Age 50-64")

# ============================================================================
# Run models and extract results
# ============================================================================
calculate_centered_r2 <- function(model) {
  y <- model$y
  y_hat <- fitted(model)
  y_mean <- mean(y, na.rm = TRUE)
  1 - sum((y - y_hat)^2, na.rm = TRUE) / sum((y - y_mean)^2, na.rm = TRUE)
}

run_model_set <- function(formula_template) {
  lapply(1:6, function(i) {
    formula_str <- gsub("OUTCOME", outcomes[i], formula_template)
    model <- ivreg(as.formula(formula_str), data = df, weights = df$timepwt48)
    
    coef_val <- coef(model)["d_tradeusch_pw"]
    vcov_cluster <- vcovCL(model, cluster = df$statefip)
    se_val <- sqrt(vcov_cluster["d_tradeusch_pw", "d_tradeusch_pw"])
    
    t_stat <- coef_val / se_val
    p_val <- 2 * pt(-abs(t_stat), df = nobs(model) - length(coef(model)))
    stars <- if (p_val < 0.01) "***" else if (p_val < 0.05) "**" else if (p_val < 0.1) "*" else ""
    
    list(
      coef = coef_val,
      se = se_val,
      stars = stars,
      n = nobs(model),
      r2 = calculate_centered_r2(model),
      fstat = summary(model, diagnostics = TRUE)$diagnostics["Weak instruments", "statistic"]
    )
  })
}

# Panel A: No Controls
results_no_controls <- run_model_set("OUTCOME ~ d_tradeusch_pw + t2 | d_tradeotch_pw_lag + t2")

# Panel B: Census Region Dummies
results_region <- run_model_set(paste0(
  "OUTCOME ~ d_tradeusch_pw + ",
  paste(reg_vars, collapse = " + "),
  " + t2 | d_tradeotch_pw_lag + ",
  paste(reg_vars, collapse = " + "),
  " + t2"
))

# Panel C: Full Controls
results_full <- run_model_set(paste0(
  "OUTCOME ~ d_tradeusch_pw + l_shind_manuf_cbp + ",
  paste(reg_vars, collapse = " + "),
  " + l_sh_popedu_c + l_sh_popfborn + l_sh_empl_f",
  " + l_sh_routine33 + l_task_outsource + t2 | ",
  "d_tradeotch_pw_lag + l_shind_manuf_cbp + ",
  paste(reg_vars, collapse = " + "),
  " + l_sh_popedu_c + l_sh_popfborn + l_sh_empl_f",
  " + l_sh_routine33 + l_task_outsource + t2"
))

# ============================================================================
# Build table
# ============================================================================
build_panel <- function(results, panel_name) {
  bind_rows(
    tibble(
      Specification = panel_name,
      All = sprintf("%.3f%s", results[[1]]$coef, results[[1]]$stars),
      College = sprintf("%.3f%s", results[[2]]$coef, results[[2]]$stars),
      `Non-College` = sprintf("%.3f%s", results[[3]]$coef, results[[3]]$stars),
      `Age 16-34` = sprintf("%.3f%s", results[[4]]$coef, results[[4]]$stars),
      `Age 35-49` = sprintf("%.3f%s", results[[5]]$coef, results[[5]]$stars),
      `Age 50-64` = sprintf("%.3f%s", results[[6]]$coef, results[[6]]$stars)
    ),
    tibble(
      Specification = "",
      All = sprintf("(%.3f)", results[[1]]$se),
      College = sprintf("(%.3f)", results[[2]]$se),
      `Non-College` = sprintf("(%.3f)", results[[3]]$se),
      `Age 16-34` = sprintf("(%.3f)", results[[4]]$se),
      `Age 35-49` = sprintf("(%.3f)", results[[5]]$se),
      `Age 50-64` = sprintf("(%.3f)", results[[6]]$se)
    ),
    tibble(
      Specification = "N",
      All = as.character(results[[1]]$n),
      College = as.character(results[[2]]$n),
      `Non-College` = as.character(results[[3]]$n),
      `Age 16-34` = as.character(results[[4]]$n),
      `Age 35-49` = as.character(results[[5]]$n),
      `Age 50-64` = as.character(results[[6]]$n)
    ),
    tibble(
      Specification = "R²",
      All = sprintf("%.3f", results[[1]]$r2),
      College = sprintf("%.3f", results[[2]]$r2),
      `Non-College` = sprintf("%.3f", results[[3]]$r2),
      `Age 16-34` = sprintf("%.3f", results[[4]]$r2),
      `Age 35-49` = sprintf("%.3f", results[[5]]$r2),
      `Age 50-64` = sprintf("%.3f", results[[6]]$r2)
    ),
    tibble(
      Specification = "First-stage F",
      All = sprintf("%.2f", results[[1]]$fstat),
      College = sprintf("%.2f", results[[2]]$fstat),
      `Non-College` = sprintf("%.2f", results[[3]]$fstat),
      `Age 16-34` = sprintf("%.2f", results[[4]]$fstat),
      `Age 35-49` = sprintf("%.2f", results[[5]]$fstat),
      `Age 50-64` = sprintf("%.2f", results[[6]]$fstat)
    )
  )
}

spacer <- tibble(Specification = "", All = "", College = "", `Non-College` = "", 
                 `Age 16-34` = "", `Age 35-49` = "", `Age 50-64` = "")

table_data <- bind_rows(
  build_panel(results_no_controls, "Panel A: No Controls"),
  spacer,
  build_panel(results_region, "Panel B: Census Region Dummies"),
  spacer,
  build_panel(results_full, "Panel C: Full Controls")
)

# ============================================================================
# Create Table
# ============================================================================
table4 <- table_data %>%
  gt() %>%
  tab_header(
    title = "Table 4: Population Change",
    subtitle = "Dependent Variable: 10-year Equivalent Log Change in Working-Age Population"
  ) %>%
  tab_spanner(label = "By Education Level", columns = c(All, College, `Non-College`)) %>%
  tab_spanner(label = "By Age Group", columns = c(`Age 16-34`, `Age 35-49`, `Age 50-64`)) %>%
  cols_label(
    Specification = "",
    All = "All",
    College = "College",
    `Non-College` = "Non-College",
    `Age 16-34` = "Age 16-34",
    `Age 35-49` = "Age 35-49",
    `Age 50-64` = "Age 50-64"
  ) %>%
  cols_align(align = "center", columns = c(All, College, `Non-College`, `Age 16-34`, `Age 35-49`, `Age 50-64`)) %>%
  cols_align(align = "left", columns = Specification) %>%
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_body(
      columns = Specification,
      rows = Specification %in% c("Panel A: No Controls", "Panel B: Census Region Dummies", "Panel C: Full Controls")
    )
  ) %>%
  tab_style(
    style = cell_fill(color = "#f0f0f0"),
    locations = cells_body(
      rows = Specification %in% c("Panel A: No Controls", "Panel B: Census Region Dummies", "Panel C: Full Controls")
    )
  ) %>%
  tab_options(
    table.font.size = px(12),
    heading.title.font.size = px(14),
    heading.subtitle.font.size = px(12),
    table_body.hlines.color = "transparent",
    table_body.border.bottom.color = "transparent",
    table_body.border.top.color = "black"
  ) %>%
  tab_footnote(footnote = "Robust standard errors in parentheses, clustered at state level.") %>%
  tab_footnote(footnote = "* p<0.1, ** p<0.05, *** p<0.01") %>%
  tab_footnote(footnote = "Full controls (Panel C) include: manufacturing share, census regions, education share, foreign-born share, female employment share, routine occupation share, and offshorability index.")

print(table4)

gt::gtsave(table4,
           filename = paste0(normalizePath("~/Desktop"), "/Table 4.png"),
           vwidth   = 6000,
           vheight  = 1600)