# ============================================================================
# Trade-Exposure Preparation
# ============================================================================

library(tidyverse)
library(haven)
library(ShiftShareSE)
library(estimatr)

# ============================================================================
# CZ-Industry Employment Share Matrix (W)
# ============================================================================
cz_emp_data <- read_dta('/Users/macca/Library/Mobile Documents/com~apple~CloudDocs/Uiversita /IV year/Empirical Econometrics/Final Paper/Public-Release-Data/dta/cbp_czone_merged.dta')

cz_shares_long <- cz_emp_data %>%
  filter(year == 1991, sic87dd >= 2000, sic87dd < 4000) %>%
  group_by(czone) %>%
  mutate(L_cz_total = sum(emp, na.rm = TRUE)) %>%
  ungroup() %>%
  filter(L_cz_total > 0) %>%
  mutate(share_cz_ind = emp / L_cz_total) %>%
  select(czone, sic_industry = sic87dd, share_cz_ind)

W_shares_df <- cz_shares_long %>%
  pivot_wider(
    names_from = sic_industry,
    values_from = share_cz_ind,
    values_fill = 0
  )

# ============================================================================
# Industry Shock Vector (Delta Z_j)
# ============================================================================
trade_data <- read_dta('/Users/macca/Library/Mobile Documents/com~apple~CloudDocs/Uiversita /IV year/Empirical Econometrics/Final Paper/Public-Release-Data/dta/sic87dd_trade_data.dta')

industry_shocks <- trade_data %>%
  filter(importer == "USA", exporter %in% c("ROW", "OTH")) %>% 
  group_by(sic87dd, year) %>%
  summarise(total_imports = sum(imports, na.rm = TRUE), .groups = 'drop') %>%
  pivot_wider(names_from = year, values_from = total_imports, values_fill = 0) %>%
  mutate(Z_shock_other = (`2007` - `1991`)) %>%
  select(sic_industry = sic87dd, Z_shock_other) %>%
  filter(sic_industry >= 2000, sic_industry < 4000)

# ============================================================================
# CZ-Level Exposure (X = W × Z)
# ============================================================================
industry_cols_numeric <- as.numeric(setdiff(names(W_shares_df), "czone"))
matched_industries <- intersect(industry_cols_numeric, industry_shocks$sic_industry)

W_subset <- W_shares_df %>%
  select(czone, all_of(as.character(matched_industries)))

Z_subset <- industry_shocks %>%
  filter(sic_industry %in% matched_industries) %>%
  arrange(sic_industry)

W_matrix_only <- as.matrix(W_subset[, -1])
Z_vector <- Z_subset$Z_shock_other
cz_exposure <- W_matrix_only %*% Z_vector

cz_exposure_df <- data.frame(
  czone = W_subset$czone,
  exposure_other = as.vector(cz_exposure)
)