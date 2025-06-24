## smoking ##

library(minfi)

# Load IDAT files and preprocess
rgSet <- read.metharray.exp("/Users/luisarivera/Documents/aging_backup/data_raw/idats")
mSet <- preprocessNoob(rgSet)  # Noob background correction
beta_values <- getBeta(mSet)   # Get beta values
beta_df <- as.data.frame(getBeta(mSet))

library(devtools)
install_github("sailalithabollepalli/EpiSmokEr")

samplesheet <- read.csv("SampleSheet_DBSCard_2023.csv")
rownames(samplesheet) <- samplesheet$Filename

result <- epismoker(dataset=beta_values, samplesheet = samplesheet, method = "SSt")
 ## two likely current smokers, one in single-exposed and one in control###


##SES, group, and bmi####

anova_result <- aov(bmi ~ group_factor, data = d)
summary(anova_result)


## Epigenetic aging by group linear models with addtional cell type PCs#######################

sHorvath1 <- lm(residuals_Horvath ~ group_factor + bio_sex +PC1 + PC2 + PC3 + PC4, data = d)
sHorvath2 <- lm(residuals_Horvath ~ ace_total + group_factor + bio_sex+ PC1+ PC2 + PC3 + PC4, data = d)

sHannum1 <- lm(residuals_Hannum ~ group_factor + bio_sex+ PC1+ PC2 + PC3 + PC4, data = d)
sHannum2 <- lm(residuals_Hannum ~ ace_total + group_factor + bio_sex+PC1+ PC2 + PC3 + PC4, data = d)


sPheno1 <- lm(residuals_PhenoAge ~ group_factor + bio_sex+PC1+ PC2 + PC3 + PC4, data = d)
sPheno2 <- lm(residuals_PhenoAge ~ ace_total + group_factor + bio_sex+PC1+ PC2 + PC3 + PC4, data = d)


sYingAdaptAge1 <- lm(residuals_YingAdaptAge~ group_factor + bio_sex+PC1+ PC2 + PC3 + PC4, data = d)
sYingAdaptAge2 <- lm(residuals_YingAdaptAge ~ ace_total + group_factor + bio_sex+PC1+ PC2 + PC3 + PC4, data = d)

sYingDamAge1 <- lm(residuals_YingDamAge~ group_factor + bio_sex+PC1+ PC2 + PC3 + PC4, data = d)
sYingDamAge2 <- lm(residuals_YingDamAge ~ ace_total + group_factor + bio_sex+PC1+ PC2 + PC3 + PC4, data = d)

sDunedin1 <- lm(residuals_dunedin~ group_factor + bio_sex+PC1+ PC2 + PC3 + PC4, data = d)
sDunedin2 <- lm(residuals_dunedin ~ ace_total + group_factor + bio_sex+PC1+ PC2 + PC3 + PC4, data = d)

sGrim1 <- lm(grim~ group_factor + bio_sex+PC1+ PC2 + PC3 + PC4, data = d)
sGrim2 <- lm(grim ~ ace_total + group_factor + bio_sex+PC1+ PC2 + PC3 + PC4, data = d)


## Save as table

sft1 <- as_flextable(sHorvath1)
sft2 <- as_flextable(sHorvath2)
sft3 <- as_flextable(sHannum1)
sft4 <- as_flextable(sHannum2)
sft5 <- as_flextable(sPheno1)
sft6 <- as_flextable(sPheno2)

save_as_docx(
  `Horvath Model 1` = sft1, `Horvath Model 2` = sft2,`Hannum Model 1` = sft3, `Hannum Model 2` = sft4,`PhenoAge Model 1` = sft5, `PhenoAge Model 2` = sft6,
  path ="figures_tables/First_gen_clock_models_sens.docx")

sft1 <- as_flextable(sYingDamAge1)
sft2 <- as_flextable(sYingDamAge2)
sft3 <- as_flextable(sYingAdaptAge1)
sft4 <- as_flextable(sYingAdaptAge2)
sft5 <- as_flextable(sGrim1)
sft6 <- as_flextable(sGrim1)
sft7 <- as_flextable(sDunedin1)  
sft8 <- as_flextable(sDunedin2)


save_as_docx(
  `YingDamAge  Model 1` = sft1, `YingDamAge Model 2` = sft2,`YingAdaptAge Model 1` = sft3, `YingAdaptAge Model 2` = sft4,`GrimAgeAccel Model 1` = sft5, `GrimAgeAccel Model 2` = sft6,
  `DunedinPACE Model 1` = sft7, `DunedinPACE Model 2` = sft8,
  path ="figures_tables/second_gen_models_sens.docx")


#######Interaction by sex ############


isHorvath1 <- lm(residuals_Horvath ~ group_factor * bio_sex +PC1, data = d)
isHorvath2 <- lm(residuals_Horvath ~ ace_total + group_factor * bio_sex+ PC1, data = d)

isHannum1 <- lm(residuals_Hannum ~ group_factor * bio_sex+ PC1, data = d)
isHannum2 <- lm(residuals_Hannum ~ ace_total + group_factor * bio_sex+PC1, data = d)


isPheno1 <- lm(residuals_PhenoAge ~ group_factor * bio_sex+PC1, data = d)
isPheno2 <- lm(residuals_PhenoAge ~ ace_total + group_factor *bio_sex+PC1, data = d)


isYingAdaptAge1 <- lm(residuals_YingAdaptAge~ group_factor * bio_sex+PC1, data = d)
isYingAdaptAge2 <- lm(residuals_YingAdaptAge ~ ace_total + group_factor * bio_sex+PC1, data = d)

isYingDamAge1 <- lm(residuals_YingDamAge~ group_factor * bio_sex+PC1, data = d)
isYingDamAge2 <- lm(residuals_YingDamAge ~ ace_total + group_factor * bio_sex+PC1, data = d)

isDunedin1 <- lm(residuals_dunedin~ group_factor * bio_sex+PC1, data = d)
isDunedin2 <- lm(residuals_dunedin ~ ace_total + group_factor * bio_sex+PC1, data = d)

isGrim1 <- lm(grim~ group_factor * bio_sex+PC1, data = d)
isGrim2 <- lm(grim ~ ace_total + group_factor * bio_sex+PC1, data = d)


## Save as table

isft1 <- as_flextable(isHorvath1)
isft2 <- as_flextable(isHorvath2)
isft3 <- as_flextable(isHannum1)
isft4 <- as_flextable(isHannum2)
isft5 <- as_flextable(isPheno1)
isft6 <- as_flextable(isPheno2)
isft7 <- as_flextable(isYingDamAge1)
isft8 <- as_flextable(isYingDamAge2)
isft9 <- as_flextable(isYingAdaptAge1)
isft10 <- as_flextable(isYingAdaptAge2)
isft11 <- as_flextable(isGrim1)
isft12 <- as_flextable(isGrim2)
isft13 <- as_flextable(isDunedin1)  
isft14 <- as_flextable(isDunedin2)

save_as_docx(
  `Horvath Model 1` = isft1, `Horvath Model 2` = isft2,`Hannum Model 1` = isft3, `Hannum Model 2` = isft4,`PhenoAge Model 1` = isft5, `PhenoAge Model 2` = isft6,  `YingDamAge  Model 1` = isft7, `YingDamAge Model 2` = isft8,`YingAdaptAge Model 1 ` = isft9, `YingAdaptAge Model 2` = isft10,`GrimAgeAccel Model 1` = isft11, `GrimAgeAccel Model 2` = isft12,
  `DunedinPACE Model 1` = isft13, `DunedinPACE Model 2` = isft14,
  path ="figures_tables/is_models.docx")



#######Sex Specific Analyses ##########
# Split dataset by bio sex
d_male <- subset(d, bio_sex == "Male")
d_female <- subset(d, bio_sex == "Female")

# Function to fit models
fit_models <- function(data) {
  list(
    Horvath1 = lm(residuals_Horvath ~ group_factor + PC1, data = data),
    Horvath2 = lm(residuals_Horvath ~ ace_total + group_factor + PC1, data = data),
    Hannum1 = lm(residuals_Hannum ~ group_factor + PC1, data = data),
    Hannum2 = lm(residuals_Hannum ~ ace_total + group_factor + PC1, data = data),
    Pheno1 = lm(residuals_PhenoAge ~ group_factor + PC1, data = data),
    Pheno2 = lm(residuals_PhenoAge ~ ace_total + group_factor + PC1, data = data),
    YingAdaptAge1 = lm(residuals_YingAdaptAge ~ group_factor + PC1, data = data),
    YingAdaptAge2 = lm(residuals_YingAdaptAge ~ ace_total + group_factor + PC1, data = data),
    YingDamAge1 = lm(residuals_YingDamAge ~ group_factor + PC1, data = data),
    YingDamAge2 = lm(residuals_YingDamAge ~ ace_total + group_factor + PC1, data = data),
    Dunedin1 = lm(residuals_dunedin ~ group_factor + PC1, data = data),
    Dunedin2 = lm(residuals_dunedin ~ ace_total + group_factor + PC1, data = data),
    Grim1 = lm(grim ~ group_factor + PC1, data = data),
    Grim2 = lm(grim ~ ace_total + group_factor + PC1, data = data)
  )
}

# Fit models for each sex
models_male <- fit_models(d_male)
models_female <- fit_models(d_female)


# Function to convert models to flextables
convert_to_flextable <- function(models) {
  lapply(models, as_flextable)
}

# Convert models to tables
ft_male <- convert_to_flextable(models_male)
ft_female <- convert_to_flextable(models_female)

# Save results as Word documents
save_as_docx(
  `Male - Horvath Model 1` = ft_male$Horvath1, `Male - Horvath Model 2` = ft_male$Horvath2,
  `Male - Hannum Model 1` = ft_male$Hannum1, `Male - Hannum Model 2` = ft_male$Hannum2,
  `Male - PhenoAge Model 1` = ft_male$Pheno1, `Male - PhenoAge Model 2` = ft_male$Pheno2,
  path ="figures_tables/male_first_gen_clock_models.docx"
)

save_as_docx(
  `Male - YingDamAge Model 1` = ft_male$YingDamAge1, `Male - YingDamAge Model 2` = ft_male$YingDamAge2,
  `Male - YingAdaptAge Model 1` = ft_male$YingAdaptAge1, `Male - YingAdaptAge Model 2` = ft_male$YingAdaptAge2,
  `Male - GrimAgeAccel Model 1` = ft_male$Grim1, `Male - GrimAgeAccel Model 2` = ft_male$Grim2,
  `Male - DunedinPACE Model 1` = ft_male$Dunedin1, `Male - DunedinPACE Model 2` = ft_male$Dunedin2,
  path ="figures_tables/male_second_gen_models.docx"
)

save_as_docx(
  `Female - Horvath Model 1` = ft_female$Horvath1, `Female - Horvath Model 2` = ft_female$Horvath2,
  `Female - Hannum Model 1` = ft_female$Hannum1, `Female - Hannum Model 2` = ft_female$Hannum2,
  `Female - PhenoAge Model 1` = ft_female$Pheno1, `Female - PhenoAge Model 2` = ft_female$Pheno2,
  path ="figures_tables/female_first_gen_clock_models.docx"
)

save_as_docx(
  `Female - YingDamAge Model 1` = ft_female$YingDamAge1, `Female - YingDamAge Model 2` = ft_female$YingDamAge2,
  `Female - YingAdaptAge Model 1` = ft_female$YingAdaptAge1, `Female - YingAdaptAge Model 2` = ft_female$YingAdaptAge2,
  `Female - GrimAgeAccel Model 1` = ft_female$Grim1, `Female - GrimAgeAccel Model 2` = ft_female$Grim2,
  `Female - DunedinPACE Model 1` = ft_female$Dunedin1, `Female - DunedinPACE Model 2` = ft_female$Dunedin2,
  path ="figures_tables/female_second_gen_models.docx"
)


