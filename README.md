# Visualizing measures of departure from MH, CoME, and EMH for Collapsed Square Contingency Tables

This repository provides R implementations of the measures proposed in the paper

> *“Decomposition of measure of departure from marginal homogeneity for collapsed ordinal square contingency tables”* 

The code allows reproduction of numerical results and figures reported in the paper.

---

## Contents

This repository contains two R scripts:

- `gh01_functions.R` 
Implements all mathematical definitions, including 
- collapsed tables \(G^{(s,t)}\) 
- measures \(\Psi_{\mathrm{MH}}, \Psi_{\mathrm{CoME}}, \Psi_{\mathrm{EMH}}\) 
- asymptotic standard errors and confidence intervals 
- visualization functions for Figures 2 and 3 

- `gh02_reproduce_paper.R` 
Reproduces tables and figures in the paper using 
- artificial data (Figure 2) 
- Japan data (Table 3a, Figure 3) 
- Britain data (Table 3b, Figure 3)

---

## Requirements

The following R packages are required:

- `ggplot2`
- `tidyr`

Install them in R if necessary:

```r
install.packages(c("ggplot2", "tidyr"))
