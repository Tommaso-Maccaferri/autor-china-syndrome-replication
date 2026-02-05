# ============================================================================
# Table 5: Change in Employment, Unemployment and Non-Employment
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
outcomes_panel_a <- c("lnchg_no_empl_mfg", "lnchg_no_empl_nmfg", "lnchg_no_unempl", 
                      "lnchg_no_nilf", "lnchg_no_ssadiswkrs")
outcomes_panel_b <- c("d_sh_empl_mfg", "d_sh_empl_nmfg", "d_sh_unempl", 
                      "d_sh_nilf", "d_sh_ssadiswkrs")
outcomes_panel_c <- c("d_sh_empl_mfg_edu_c", "d_sh_empl_nmfg_edu_c", 
                      "d_sh_unempl_edu_c", "d_sh_nilf_edu_c")
outcomes_panel_d <- c("d_sh_empl_mfg_edu_nc", "d_sh_empl_nmfg_edu_nc", 
                      "d_sh_unempl_edu_nc", "d_sh_nilf_edu_nc")

col_labels <- c("Mfg Emp", "Non-Mfg Emp", "Unemployed", "Not in LF", "SSDI Receipt")

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
  n <- length(results)
  
  get_vals <- function(field, formatter = identity) {
    sapply(1:5, function(i) {
      if (i <= n && !is.null(results[[i]])) formatter(results[[i]][[field]]) else if (field == "coef") "—" else ""
    })
  }
  
  bind_rows(
    tibble(
      Specification = panel_name,
      `Mfg Emp` = if (n >= 1) sprintf("%.3f%s", results[[1]]$coef, results[[1]]$stars) else "—",
      `Non-Mfg Emp` = if (n >= 2) sprintf("%.3f%s", results[[2]]$coef, results[[2]]$stars) else "—",
      Unemployed = if (n >= 3) sprintf("%.3f%s", results[[3]]$coef, results[[3]]$stars) else "—",
      `Not in LF` = if (n >= 4) sprintf("%.3f%s", results[[4]]$coef, results[[4]]$stars) else "—",
      `SSDI Receipt` = if (n >= 5) sprintf("%.3f%s", results[[5]]$coef, results[[5]]$stars) else "—"
    ),
    tibble(
      Specification = "",
      `Mfg Emp` = get_vals("se", function(x) sprintf("(%.3f)", x))[1],
      `Non-Mfg Emp` = get_vals("se", function(x) sprintf("(%.3f)", x))[2],
      Unemployed = get_vals("se", function(x) sprintf("(%.3f)", x))[3],
      `Not in LF` = get_vals("se", function(x) sprintf("(%.3f)", x))[4],
      `SSDI Receipt` = if (n >= 5) sprintf("(%.3f)", results[[5]]$se) else ""
    ),
    tibble(
      Specification = "N",
      `Mfg Emp` = get_vals("n", as.character)[1],
      `Non-Mfg Emp` = get_vals("n", as.character)[2],
      Unemployed = get_vals("n", as.character)[3],
      `Not in LF` = get_vals("n", as.character)[4],
      `SSDI Receipt` = if (n >= 5) as.character(results[[5]]$n) else ""
    ),
    tibble(
      Specification = "R²",
      `Mfg Emp` = get_vals("r2", function(x) sprintf("%.3f", x))[1],
      `Non-Mfg Emp` = get_vals("r2", function(x) sprintf("%.3f", x))[2],
      Unemployed = get_vals("r2", function(x) sprintf("%.3f", x))[3],
      `Not in LF` = get_vals("r2", function(x) sprintf("%.3f", x))[4],
      `SSDI Receipt` = if (n >= 5) sprintf("%.3f", results[[5]]$r2) else ""
    )
  )
}

spacer <- tibble(Specification = "", `Mfg Emp` = "", `Non-Mfg Emp` = "", 
                 Unemployed = "", `Not in LF` = "", `SSDI Receipt` = "")

table_data <- bind_rows(
  build_panel_rows(results_a, "Panel A: 100 × Log Change in Population Counts"),
  spacer,
  build_panel_rows(results_b, "Panel B: All Education Levels"),
  spacer,
  build_panel_rows(results_c, "Panel C: College Education"),
  spacer,
  build_panel_rows(results_d, "Panel D: No College Education")
)

# ============================================================================
# Create GT table
# ============================================================================
table5 <- table_data %>%
  gt() %>%
  tab_header(
    title = "Table 5: Employment Status of Working-Age Population within CZs, 1990-2007",
    subtitle = "Dependent Variables: Ten-year Equivalent Changes (2SLS Estimates)"
  ) %>%
  cols_label(
    Specification = "",
    `Mfg Emp` = "Mfg Emp",
    `Non-Mfg Emp` = "Non-Mfg Emp",
    Unemployed = "Unemployed",
    `Not in LF` = "Not in LF",
    `SSDI Receipt` = "SSDI Receipt"
  ) %>%
  cols_align(align = "center", columns = c(`Mfg Emp`, `Non-Mfg Emp`, Unemployed, `Not in LF`, `SSDI Receipt`)) %>%
  cols_align(align = "left", columns = Specification) %>%
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_body(
      columns = Specification,
      rows = Specification %in% c("Panel A: 100 × Log Change in Population Counts", "Panel B: All Education Levels",
                                  "Panel C: College Education", "Panel D: No College Education")
    )
  ) %>%
  tab_style(
    style = cell_fill(color = "#f0f0f0"),
    locations = cells_body(
      rows = Specification %in% c("Panel A: 100 × Log Change in Population Counts", "Panel B: All Education Levels",
                                  "Panel C: College Education", "Panel D: No College Education")
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

print(table5)

gt::gtsave(table5,
           filename = paste0(normalizePath("~/Desktop"), "/Table 5.png"),
           vwidth   = 6000,
           vheight  = 1600)