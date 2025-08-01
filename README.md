# Association of Vitamins and Micronutrients with Type 1 Diabetes Risk: A Mendelian Randomization Study 
This repository contains code for a paper on Mendelian Randomization (MR) analysis investigating the causal effects of micronutrients on Type 1 Diabetes (T1D) risk.

# Abstract:

Background/Objectives: Previous studies suggest that nutrient deficiencies can alter immune responses in animals. However, the impact of micronutrients on autoimmune diseases like type 1 diabetes (T1D) in humans remains unclear since the described associations are based on observational data and as such, they cannot establish causality. This study aims to examine the causal relationship between various micronutrients and T1D using Mendelian randomization (MR). 

Methods: We performed a two-sample MR analysis using genetic variants from genome-wide association studies (GWAS) of 17 micronutrients as instrumental variables (IVs). We analyzed European (18,942 cases/520,580 controls), multi-ancestry (25,717 cases/583,311 controls), Latin American/Hispanic (2,295 cases/55,134 controls), African American/Afro-Carribean (6,451 cases/109,410 controls), and East Asian (1,219 cases/132,032 controls) T1D GWAS datasets. We applied the inverse variance weighted (IVW) analysis, additional MR estimators (MR-Egger, weighted median, weighted mode, MR-PRESSO) to address pleiotropy, and the Steiger test for directionality. 

Results: The MR analysis showed a positive association between potassium levels and T1D risk, sustained after Bonferroni correction (p < 0.0029, OR = 1.098, 95% CI [1.075, 1.122] in the multi-ancestry T1D GWAS). Zinc, vitamin B12, retinol, and alpha-tocopherol showed nominal associations. Vitamin C, D, K1, B6, beta- and gamma-tocopherol, magnesium, iron, copper, selenium, carotene, and folate showed no significant effect on T1D risk. For the multi-ancestry GWAS, ORs required for 80% power for the micronutrients ranged from 1.15 to 3.35, indicating limited power to detect small effects. 

Conclusions: Higher serum potassium levels were associated with increased T1D risk in our MR study, though supporting observational evidence is currently limited. Other micronutrients are unlikely to have large effects on T1D. 

# GITHUB content
Main_MR/: Core MR analysis scripts
- TwoSampleMR.R: Implements main MR methods (e.g., inverse-variance weighted, MR-Egger)
- Verma_proxy_search.R: Proxy SNP search for the Verma et al. T1D GWAS
- proxy_search.R: General proxy SNP search for MR analysis

Supplementary_MR/: Supplementary MR methods and visualizations
- MRlap_method.r: MRLap analysis for integrative MR
- MR_PRESSO.r: MR-PRESSO for pleiotropy correction
- Steiger_test.r: Steiger test for causal direction
- MR_plot.r: Generates supplementary plots

GWAS_processing/: GWAS data processing scripts
- rsid_extraction.r: RSID extraction for the Sakaue et al. T1D GWAS
