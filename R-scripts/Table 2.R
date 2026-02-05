# ============================================================================
# Table 2: Pre-Period Placebo Test
# ============================================================================

library(haven)
library(dplyr)
library(AER)
library(sandwich)
library(gt)
library(tibble)

# Load data
df <- read_dta('workfile_china_preperiod.dta')

# ============================================================================
# Extract results function
# ============================================================================
run_and_extract <- function(model, cluster_var, coef_name) {
  beta <- coef(model)[coef_name]
  vcov_cl <- vcovCL(model, cluster = cluster_var)
  se <- sqrt(vcov_cl[coef_name, coef_name])
  
  t_stat <- beta / se
  df_resid <- nobs(model) - length(coef(model))
  p_val <- 2 * pt(-abs(t_stat), df = df_resid)
  
  stars <- if (p_val < 0.01) "***" else if (p_val < 0.05) "**" else if (p_val < 0.1) "*" else ""
  
  beta_scaled <- beta * 0.48
  t_scaled <- beta_scaled / se
  p_scaled <- 2 * pt(-abs(t_scaled), df = df_resid)
  stars_scaled <- if (p_scaled < 0.01) "***" else if (p_scaled < 0.05) "**" else if (p_scaled < 0.1) "*" else ""
  
  list(
    coef = beta,
    se = se,
    stars = stars,
    scaled_coef = beta_scaled,
    scaled_stars = stars_scaled,
    n = nobs(model),
    r2 = summary(model)$r.squared
  )
}

# ============================================================================
# Panel A: Main Period (1990-2007)
# ============================================================================
model1 <- ivreg(d_sh_empl_mfg ~ d_tradeusch_pw | d_tradeotch_pw_lag,
                data = df[df$yr == 1990, ],
                weights = df$timepwt48[df$yr == 1990])
result_1 <- run_and_extract(model1, df$statefip[df$yr == 1990], "d_tradeusch_pw")

model2 <- ivreg(d_sh_empl_mfg ~ d_tradeusch_pw | d_tradeotch_pw_lag,
                data = df[df$yr == 2000, ],
                weights = df$timepwt48[df$yr == 2000])
result_2 <- run_and_extract(model2, df$statefip[df$yr == 2000], "d_tradeusch_pw")

model3 <- ivreg(d_sh_empl_mfg ~ d_tradeusch_pw + t2000 | d_tradeotch_pw_lag + t2000,
                data = df[df$yr >= 1990, ],
                weights = df$timepwt48[df$yr >= 1990])
result_3 <- run_and_extract(model3, df$statefip[df$yr >= 1990], "d_tradeusch_pw")

# ============================================================================
# Panel B: Pre-Period Placebo (1970-1990)
# ============================================================================
model4 <- ivreg(d_sh_empl_mfg ~ d_tradeusch_pw_future | d_tradeotch_pw_lag_future,
                data = df[df$yr == 1970, ],
                weights = df$timepwt48[df$yr == 1970])
result_4 <- run_and_extract(model4, df$statefip[df$yr == 1970], "d_tradeusch_pw_future")

model5 <- ivreg(d_sh_empl_mfg ~ d_tradeusch_pw_future | d_tradeotch_pw_lag_future,
                data = df[df$yr == 1980, ],
                weights = df$timepwt48[df$yr == 1980])
result_5 <- run_and_extract(model5, df$statefip[df$yr == 1980], "d_tradeusch_pw_future")

model6 <- ivreg(d_sh_empl_mfg ~ d_tradeusch_pw_future + t1980 | d_tradeotch_pw_lag_future + t1980,
                data = df[df$yr >= 1970 & df$yr < 1990, ],
                weights = df$timepwt48[df$yr >= 1970 & df$yr < 1990])
result_6 <- run_and_extract(model6, df$statefip[df$yr >= 1970 & df$yr < 1990], "d_tradeusch_pw_future")

# ============================================================================
# Build table
# ============================================================================
build_panel_rows <- function(results, panel_name) {
  bind_rows(
    tibble(
      Specification = panel_name,
      Col1 = sprintf("%.3f%s", results[[1]]$coef, results[[1]]$stars),
      Col2 = sprintf("%.3f%s", results[[2]]$coef, results[[2]]$stars),
      Col3 = sprintf("%.3f%s", results[[3]]$coef, results[[3]]$stars)
    ),
    tibble(
      Specification = "[Adjusted coefficient]",
      Col1 = sprintf("[%.3f%s]", results[[1]]$scaled_coef, results[[1]]$scaled_stars),
      Col2 = sprintf("[%.3f%s]", results[[2]]$scaled_coef, results[[2]]$scaled_stars),
      Col3 = sprintf("[%.3f%s]", results[[3]]$scaled_coef, results[[3]]$scaled_stars)
    ),
    tibble(
      Specification = "",
      Col1 = sprintf("(%.3f)", results[[1]]$se),
      Col2 = sprintf("(%.3f)", results[[2]]$se),
      Col3 = sprintf("(%.3f)", results[[3]]$se)
    ),
    tibble(
      Specification = "N",
      Col1 = as.character(results[[1]]$n),
      Col2 = as.character(results[[2]]$n),
      Col3 = as.character(results[[3]]$n)
    ),
    tibble(
      Specification = "R²",
      Col1 = sprintf("%.3f", results[[1]]$r2),
      Col2 = sprintf("%.3f", results[[2]]$r2),
      Col3 = sprintf("%.3f", results[[3]]$r2)
    )
  )
}

panel_a_data <- build_panel_rows(list(result_1, result_2, result_3), "Panel A: Main Period (1990-2007)")
panel_b_data <- build_panel_rows(list(result_4, result_5, result_6), "Panel B: Pre-Period Placebo (1970-1990)")
spacer <- tibble(Specification = "", Col1 = "", Col2 = "", Col3 = "")

table_data <- bind_rows(panel_a_data, spacer, panel_b_data)
colnames(table_data) <- c("Specification", "Period 1", "Period 2", "Pooled")

# ============================================================================
# Create Table
# ============================================================================
table2 <- table_data %>%
  gt() %>%
  tab_header(
    title = "Table 2: Pre-Period Placebo Test - Manufacturing Employment Changes",
    subtitle = "Dependent Variable: Change in Manufacturing Employment Share (Δ Mfg Empl/Pop)"
  ) %>%
  cols_label(
    Specification = "",
    `Period 1` = "1990-2000",
    `Period 2` = "2000-2007",
    Pooled = "Pooled"
  ) %>%
  cols_align(align = "left", columns = Specification) %>%
  cols_align(align = "center", columns = c(`Period 1`, `Period 2`, Pooled)) %>%
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_body(
      columns = Specification,
      rows = Specification %in% c("Panel A: Main Period (1990-2007)", "Panel B: Pre-Period Placebo (1970-1990)")
    )
  ) %>%
  tab_style(
    style = cell_fill(color = "#f0f0f0"),
    locations = cells_body(
      rows = Specification %in% c("Panel A: Main Period (1990-2007)", "Panel B: Pre-Period Placebo (1970-1990)")
    )
  ) %>%
  tab_style(
    style = cell_text(style = "italic", size = px(10)),
    locations = cells_body(columns = Specification, rows = Specification == "[Adjusted coefficient]")
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
  tab_footnote(
    footnote = "Panel A shows the main period results (should be negative and significant). Panel B shows placebo tests using future import exposure on past employment changes (should be insignificant)."
  ) %>%
  tab_footnote(
    footnote = "* p<0.1, ** p<0.05, *** p<0.01. Standard errors in parentheses, clustered at state level."
  )

print(table2)

#gt::gtsave(table2,
#           filename = paste0(normalizePath("~/Desktop"), "/Table 2.png"),
#          vwidth   = 6000,
#         vheight  = 1600)