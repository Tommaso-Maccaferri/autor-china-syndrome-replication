# ============================================================================
# Table 7: Manufacturing vs. Non-Manufacturing Employment and Wages
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

# ============================================================================
# Define outcomes
# ============================================================================
outcomes_panel_a <- c("lnchg_no_empl_mfg", "lnchg_no_empl_mfg_edu_c", "lnchg_no_empl_mfg_edu_nc")
outcomes_panel_b <- c("lnchg_no_empl_nmfg", "lnchg_no_empl_nmfg_edu_c", "lnchg_no_empl_nmfg_edu_nc")
outcomes_panel_c <- c("d_avg_lnwkwage_mfg", "d_avg_lnwkwage_mfg_c", "d_avg_lnwkwage_mfg_nc")
outcomes_panel_d <- c("d_avg_lnwkwage_nmfg", "d_avg_lnwkwage_nmfg_c", "d_avg_lnwkwage_nmfg_nc")

# ============================================================================
# Run models and extract results
# ============================================================================
build_iv_formula <- function(outcome_var) {
  as.formula(paste(
    outcome_var, "~ d_tradeusch_pw + l_shind_manuf_cbp +",
    paste(reg_vars, collapse = " + "),
    "+ l_sh_popedu_c + l_sh_popfborn + l_sh_empl_f",
    "+ l_sh_routine33 + l_task_outsource + t2 |",
    "d_tradeotch_pw_lag + l_shind_manuf_cbp +",
    paste(reg_vars, collapse = " + "),
    "+ l_sh_popedu_c + l_sh_popfborn + l_sh_empl_f",
    "+ l_sh_routine33 + l_task_outsource + t2"
  ))
}

run_panel_models <- function(outcomes_list) {
  lapply(outcomes_list, function(outcome) {
    model <- ivreg(build_iv_formula(outcome), data = df, weights = df$timepwt48)
    
    beta <- coef(model)["d_tradeusch_pw"]
    vcov_cl <- vcovCL(model, cluster = df$statefip)
    se <- sqrt(vcov_cl["d_tradeusch_pw", "d_tradeusch_pw"])
    
    df_resid <- nobs(model) - length(coef(model))
    t_stat <- beta / se
    p_val <- 2 * pt(-abs(t_stat), df = df_resid)
    stars <- if (p_val < 0.01) "***" else if (p_val < 0.05) "**" else if (p_val < 0.1) "*" else ""
    
    list(
      coef = beta,
      se = se,
      stars = stars,
      n = nobs(model),
      r2 = summary(model)$r.squared
    )
  })
}

results_a <- run_panel_models(outcomes_panel_a)
results_b <- run_panel_models(outcomes_panel_b)
results_c <- run_panel_models(outcomes_panel_c)
results_d <- run_panel_models(outcomes_panel_d)

# ============================================================================
# Build table
# ============================================================================
build_panel_rows <- function(results, panel_name) {
  bind_rows(
    tibble(
      Specification = panel_name,
      `All Workers` = sprintf("%.3f%s", results[[1]]$coef, results[[1]]$stars),
      College = sprintf("%.3f%s", results[[2]]$coef, results[[2]]$stars),
      `Non-College` = sprintf("%.3f%s", results[[3]]$coef, results[[3]]$stars)
    ),
    tibble(
      Specification = "",
      `All Workers` = sprintf("(%.3f)", results[[1]]$se),
      College = sprintf("(%.3f)", results[[2]]$se),
      `Non-College` = sprintf("(%.3f)", results[[3]]$se)
    ),
    tibble(
      Specification = "N",
      `All Workers` = as.character(results[[1]]$n),
      College = as.character(results[[2]]$n),
      `Non-College` = as.character(results[[3]]$n)
    ),
    tibble(
      Specification = "R²",
      `All Workers` = sprintf("%.3f", results[[1]]$r2),
      College = sprintf("%.3f", results[[2]]$r2),
      `Non-College` = sprintf("%.3f", results[[3]]$r2)
    )
  )
}

spacer <- tibble(Specification = "", `All Workers` = "", College = "", `Non-College` = "")

table_data <- bind_rows(
  build_panel_rows(results_a, "Panel A: Manufacturing Sector"),
  spacer,
  build_panel_rows(results_b, "Panel B: Nonmanufacturing Sectors"),
  spacer,
  build_panel_rows(results_c, "Panel C: Mfg Wages"),
  spacer,
  build_panel_rows(results_d, "Panel D: Non-Mfg Wages")
)

# ============================================================================
# Create GT table
# ============================================================================
table7 <- table_data %>%
  gt() %>%
  tab_header(
    title = "Table 7: Comparing Manufacturing and Nonmanufacturing Employment and Wages, 1990-2007",
    subtitle = "Dependent Variables: Ten-year Equivalent Changes (2SLS Estimates)"
  ) %>%
  tab_spanner(label = "Log Change in Number of Workers", columns = c(`All Workers`, College, `Non-College`), id = "employment") %>%
  tab_spanner(label = "Change in Average Log Weekly Wage", columns = c(`All Workers`, College, `Non-College`), id = "wages", level = 2) %>%
  cols_label(
    Specification = "",
    `All Workers` = "All Workers",
    College = "College",
    `Non-College` = "Non-College"
  ) %>%
  cols_align(align = "center", columns = c(`All Workers`, College, `Non-College`)) %>%
  cols_align(align = "left", columns = Specification) %>%
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_body(
      columns = Specification,
      rows = Specification %in% c("Panel A: Manufacturing Sector", "Panel B: Nonmanufacturing Sectors",
                                  "Panel C: Mfg Wages", "Panel D: Non-Mfg Wages")
    )
  ) %>%
  tab_style(
    style = cell_fill(color = "#f0f0f0"),
    locations = cells_body(
      rows = Specification %in% c("Panel A: Manufacturing Sector", "Panel B: Nonmanufacturing Sectors",
                                  "Panel C: Mfg Wages", "Panel D: Non-Mfg Wages")
    )
  ) %>%
  tab_options(
    table.font.size = px(11),
    heading.title.font.size = px(14),
    heading.subtitle.font.size = px(12),
    table_body.hlines.color = "transparent",
    table_body.border.bottom.color = "transparent",
    table_body.border.top.color = "black",
    row_group.border.bottom.color = "black",
    column_labels.border.bottom.color = "black"
  ) %>%
  tab_footnote(footnote = "Δ Chinese Import Exposure (per worker). Robust standard errors in parentheses, clustered at state level.") %>%
  tab_footnote(footnote = "* p<0.1, ** p<0.05, *** p<0.01") %>%
  tab_footnote(footnote = "All models include full controls: manufacturing share, census regions, education share, foreign-born share, female employment share, routine occupation share, and offshorability index.")

print(table7)

gt::gtsave(table7,
           filename = paste0(normalizePath("~/Desktop"), "/Table 7.png"),
           vwidth   = 6000,
           vheight  = 1600)