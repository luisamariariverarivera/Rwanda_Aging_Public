# Rwanda_Aging_Public

**Description:**
This repository contains code to reproduce Uwizeye et al (2024) Prenatal exposure to genocide accelerates epigenetic aging as measured in second-generation clocks among young adults.

Preprint available: https://www.medrxiv.org/content/10.1101/2024.10.01.24314372v1

**Data Availability:**
Raw data are not available due to restrictions in our consent form. Data for the purposes of reproduction only is available by request to luisa.maria.rivera@gmail.com.

Our DNA methylation preprocessing pipeline and cell-type deconvolution is not presented here but is available upon request. 

Epigenetic age estimates were generated using the Biolearn tool from the Biomarkers of Aging consortium available at: https://bio-learn.github.io/.

**Files:**

The "Processing.R" file is the basic data processing file.
  - Joins our various phenotypic and age estimate files, and contains the cell-type PCA used to control from cell-type heterogeneity
  - Calcuates aging residuals where appropriate (e.g. not for GrimAgeAccel and DunedinPACE, which calculate residuals automatically)

The "Analysis.R" file contains the analysis in the main text.
- We begin by running a series of regression models investigating the effect of exposure group on aging estimates with and without adjustment for postnatal ACES.
- We extract effect sizes and plot these to demonstrate the increased sensitivity of second generation clocks compared with first generation clocks to the prenatal environment.

The "revision_analyses.R" script contains analysis requested by reviewers after peer review. 
- We compare exposure group BMI and SES
- We use the epismoker package to estimate smoking prevalence in participants from DNA methylation data.
- We conduct sensititivity analyses with interaction by sex and with additional immune cell type principle components.
  


