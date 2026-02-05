# Replication of "The China Syndrome" (Autor et al., 2013)

This repository contains a replication study of the seminal paper *The China Syndrome: Local Labor Market Effects of Import Competition in the United States*, originally published by Autor, Dorn, and Hanson.

## The Paper
Autor et al. "The China Syndrome" replication 2025.pdf

## The Code
The analysis is performed using R. The replication covers the main tables presented in the original paper.

### Repository Structure:
* `Autor et al. "The China Syndrome" replication 2025.pdf`: The final output describing the methodology and results.
* `R-scripts`: Contains the R code for each step of the analysis.

### Script Guide:
To replicate specific tables, refer to the corresponding scripts:
* **Table 1:** `Trade-Exposure Prep.R` - *Preparation Code used to implement the Shift-Share Model*
* **Table 2:** `Table 2.R` - *Pre-Period Placebo Test - Manufacturing Employment Changes*
* **Table 3:** `Table 3.R` - *Change in Manufacturing Employment Share (Mfg Empl/Pop)*
* **Table 4:** `Table 4.R` - *Population Change*
* **Table 5:** `Table 5.R` - *Employment Status of Working-Age Population within CZs, 1990-2007*
* **Table 6:** `Table 7.R` - *Comparing Manufacturing and Nonmanufacturing Employment and Wages, 1990-2007*
* **Table 7**  `Table 10.R`- *Comprehensive Results*

### Data Gathering
* For the data used to write this replication paper please refer to the freely available and downloable data provided in the replication package in the AER page of the paper: https://www.aeaweb.org/articles?id=10.1257/aer.103.6.2121

**Tools used:** R, Tidyverse, haven, dplyr, AER, lmtest, modelsummary, sandwich

---
*This project was conducted as part of the Research Master's curriculum in Economics and Financial Research.*
