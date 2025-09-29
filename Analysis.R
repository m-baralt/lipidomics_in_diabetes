####### Packages loading #######
library(ggplot2)
library(magrittr)
library(plyr)
library(dplyr)
library(knitr)
library(viridis)
doParallel::registerDoParallel(6)

dir.create("Results")

dir_path <- ""

####### Lipidomics data preparation #######

# Lipidomics filtered data - Negative ionization
Neg_table <- readxl::read_excel(path = paste0(dir_path, "LIPIDS_NEG_RTCOR_pmp_preprocessing.xlsx"), 
                                sheet = "pqn_normalised_data")
Neg_table <- as.data.frame(Neg_table)
rownames(Neg_table) <- Neg_table$...1
Neg_table <- Neg_table[,-1]
Neg_table <- as.data.frame(t(Neg_table))

for (i in 1:ncol(Neg_table))
  Neg_table[,i] <- log10(as.numeric(as.character(Neg_table[,i])))

# Fix rownames
rownames(Neg_table) <- gsub(pattern = "OYO", replacement = "OY0", rownames(Neg_table))
rownames(Neg_table) <- gsub(pattern = "[.]", replacement = "-", rownames(Neg_table))

# Lipidomics filtered data - Positive ionization
Pos_table <- readxl::read_excel(path = paste0(dir_path, "LIPIDS_POS_RTCOR_pmp_preprocessing.xlsx"), 
                                sheet = "pqn_normalised_data")
Pos_table <- as.data.frame(Pos_table)
rownames(Pos_table) <- Pos_table$...1
Pos_table <- Pos_table[,-1]
Pos_table <- as.data.frame(t(Pos_table))

for (i in 1:ncol(Pos_table))
  Pos_table[,i] <- log10(as.numeric(as.character(Pos_table[,i])))

# Fix rownames
rownames(Pos_table) <- gsub(pattern = "OYO", replacement = "OY0", rownames(Pos_table))
rownames(Pos_table) <- gsub(pattern = "[.]", replacement = "-", rownames(Pos_table))

####### preparation of identification data #######

source("Required_functions.R")

lipidsearch_neg <- readxl::read_excel(path = paste0(dir_path, "LIPIDS_NEG_RTCOR_pmp_preprocessing.xlsx"), 
                                      sheet = "Lipids_search_matched")
lipidsearch_neg <- as.data.frame(lipidsearch_neg)
rownames(lipidsearch_neg) <- lipidsearch_neg$...1
lipidsearch_neg <- lipidsearch_neg[,-1]

lipidsearch_neg <- lipidsearch_prep(df = lipidsearch_neg)

lipidsearch_pos <- readxl::read_excel(path = paste0(dir_path, "LIPIDS_POS_RTCOR_pmp_preprocessing.xlsx"), 
                                      sheet = "Lipids_search_matched")
lipidsearch_pos <- as.data.frame(lipidsearch_pos)
rownames(lipidsearch_pos) <- lipidsearch_pos$...1
lipidsearch_pos <- lipidsearch_pos[,-1]

lipidsearch_pos <- lipidsearch_prep(df = lipidsearch_pos)

lipidsearch_pos$NCarbons[lipidsearch_pos$uniqueID %in% "M385T416_ZyE(0:0)+H"] <- 0
lipidsearch_pos$NBonds[lipidsearch_pos$uniqueID %in% "M385T416_ZyE(0:0)+H"] <- 0

##### clinical variables #######
clinical_vars <- read.csv(file = "data/clinical_variables.csv", 
                          row.names = 1, stringsAsFactors = TRUE)

# Observations with missing values in the variables used as confounding factors are removed
clinical_vars <- clinical_vars[!(is.na(clinical_vars$BMI) | 
                                   is.na(clinical_vars$Glucose) |
                                   is.na(clinical_vars$systolic_BP) |
                                   is.na(clinical_vars$aMED) |
                                   is.na(clinical_vars$MDRD_4_IDMS)),]

# Merge data
merge_pos <- merge(clinical_vars, Pos_table, by = 0)
rownames(merge_pos) <- merge_pos$Row.names
merge_pos <- merge_pos[,-1]

merge_neg <- merge(clinical_vars, Neg_table, by = 0)
rownames(merge_neg) <- merge_neg$Row.names
merge_neg <- merge_neg[,-1]

table(clinical_vars$DM, clinical_vars$Origen)

# Only the samples from Lleida and Mollerusa are considered
merge_pos <- merge_pos[merge_pos$Origen %in% c("Lleida", "Mollerussa"),]
merge_neg <- merge_neg[merge_neg$Origen %in% c("Lleida", "Mollerussa"),]

## prediabetes separation
merge_pos$DM <- as.character(merge_pos$DM)
merge_pos$DM[(merge_pos$DM %in% "Control") & ((merge_pos$HbA1c>=5.7) | (merge_pos$Glucose>=100))] <- "Prediabetes"
merge_pos$DM <- as.factor(merge_pos$DM)
table(merge_pos$DM)
nonScaled_Pos <- merge_pos

merge_pos <- merge_pos[!(merge_pos$DM %in% "Prediabetes"),]
merge_pos$DM <- as.factor(as.character(merge_pos$DM))

merge_neg$DM <- as.character(merge_neg$DM)
merge_neg$DM[(merge_neg$DM %in% "Control") & ((merge_neg$HbA1c>=5.7) | (merge_neg$Glucose>=100))] <- "Prediabetes"
merge_neg$DM <- as.factor(merge_neg$DM)
table(merge_neg$DM)
nonScaled_Neg <- merge_neg
merge_neg <- merge_neg[!(merge_neg$DM %in% "Prediabetes"),]
merge_neg$DM <- as.factor(as.character(merge_neg$DM))

nonScaled_Pos <- merge_pos
nonScaled_Neg <- merge_neg

# Statistical Analysis
feat_cols_pos <- colnames(merge_pos)[-c(1:26)]
feat_cols_neg <- colnames(merge_neg)[-c(1:26)]

numeric_cols <- c("AGE", "BMI", "Waist", "DM_duration", "HbA1c", "Glucose", 
                  "systolic_BP", "diastolic_BP", "Triglycerides", "HDL_cholesterol", 
                  "LDL_cholesterol", "aMED", "MDRD_4_IDMS", "Total_number_plaques", 
                  "Total_cholesterol", "GPT_ALT")

merge_pos[,c(numeric_cols, feat_cols_pos)] <- scale(
  x = merge_pos[,c(numeric_cols, feat_cols_pos)], center = TRUE, scale = TRUE)

merge_neg[,c(numeric_cols, feat_cols_neg)] <- scale(
  x = merge_neg[,c(numeric_cols, feat_cols_neg)], center = TRUE, scale = TRUE)

######################################################
##################### T1D vs T2D #####################
######################################################

###################### Positive ######################

confounders_str <- paste(
  "Sex", "AGE", "Hypertension", "Dyslipidemia", "DM", "BMI", "DM_duration",
  "HbA1c", "Glucose", "Smoking", "systolic_BP", "Triglycerides",
  "HDL_cholesterol", "LDL_cholesterol", "aMED", "MDRD_4_IDMS", sep = "+")

pvals_DM1_DM2 <- linear_analysis(data = merge_pos, feat_cols = feat_cols_pos, 
                                 confounders_str = confounders_str, 
                                 group_subjects = c("T1D", "T2D"), 
                                 doPar = TRUE)

merge_pvals_DM1_DM2_pos <- merge(pvals_DM1_DM2, 
                                 lipidsearch_pos[,c("Feature", "LipidIon", "Class", 
                                                    "NCarbons", "NBonds","Grade", 
                                                    "lipidsearch_mz_ppm", "XCMS_rt", 
                                                    "lipidsearch_rt", "diff_rt", 
                                                    "confirmed")], 
                                 all.x = TRUE)

write.csv(merge_pvals_DM1_DM2_pos, "Results/pvalues/DM1_DM2_pos.csv")

pvals_DM1_DM2_sex <- linear_analysis_interactions(
  data = merge_pos, feat_cols = feat_cols_pos, 
  confounders_str = confounders_str, 
  group_subjects = c("T1D", "T2D"), doPar = TRUE)

merge_pvals_DM1_DM2_pos_sex <- merge(
  pvals_DM1_DM2_sex, 
  lipidsearch_pos[,c("Feature", "LipidIon", "Class", "NCarbons", "NBonds","Grade", 
                     "lipidsearch_mz_ppm", "XCMS_rt", "lipidsearch_rt", 
                     "diff_rt", "confirmed")], 
  all.x = TRUE)

write.csv(merge_pvals_DM1_DM2_pos_sex, "Results/pvalues/DM1_DM2_pos_sex.csv")

###################### Negative ######################

pvals_DM1_DM2 <- linear_analysis(data = merge_neg, feat_cols = feat_cols_neg, 
                                 confounders_str = confounders_str, 
                                 group_subjects = c("T1D", "T2D"), 
                                 doPar = TRUE)

merge_pvals_DM1_DM2_neg <- merge(pvals_DM1_DM2, 
                                 lipidsearch_neg[,c("Feature", "LipidIon", "Class", 
                                                    "NCarbons", "NBonds","Grade", 
                                                    "lipidsearch_mz_ppm", "XCMS_rt", 
                                                    "lipidsearch_rt", "diff_rt", 
                                                    "confirmed")], 
                                 all.x = TRUE)

write.csv(merge_pvals_DM1_DM2_neg, "Results/pvalues/DM1_DM2_neg.csv")

pvals_DM1_DM2_sex <- linear_analysis_interactions(
  data = merge_neg, feat_cols = feat_cols_neg, 
  confounders_str = confounders_str, 
  group_subjects = c("T1D", "T2D"), doPar = TRUE)

merge_pvals_DM1_DM2_neg_sex <- merge(
  pvals_DM1_DM2_sex, 
  lipidsearch_neg[,c("Feature", "LipidIon", "Class", "NCarbons", "NBonds","Grade", 
                     "lipidsearch_mz_ppm", "XCMS_rt", "lipidsearch_rt", 
                     "diff_rt", "confirmed")], 
  all.x = TRUE)

write.csv(merge_pvals_DM1_DM2_neg_sex, "Results/pvalues/DM1_DM2_neg_sex.csv")

######################################################
##################### CT vs T1D ######################
######################################################

###################### Positive ######################

confounders_str <- paste(
  "Sex", "AGE", "Hypertension", "Dyslipidemia", "DM", "BMI", 
  "Glucose", "Smoking", "systolic_BP", "Triglycerides",
  "HDL_cholesterol", "LDL_cholesterol", "aMED", "MDRD_4_IDMS", sep = "+")

pvals_CT_DM1 <- linear_analysis(data = merge_pos, feat_cols = feat_cols_pos, 
                                 confounders_str = confounders_str, 
                                 group_subjects = c("Control", "T1D"), 
                                 doPar = TRUE)

merge_pvals_CT_DM1_pos <- merge(pvals_CT_DM1, 
                                 lipidsearch_pos[,c("Feature", "LipidIon", "Class", 
                                                    "NCarbons", "NBonds","Grade", 
                                                    "lipidsearch_mz_ppm", "XCMS_rt", 
                                                    "lipidsearch_rt", "diff_rt", 
                                                    "confirmed")], 
                                 all.x = TRUE)

write.csv(merge_pvals_CT_DM1_pos, "Results/pvalues/CT_DM1_pos.csv")

pvals_CT_DM1_sex <- linear_analysis_interactions(
  data = merge_pos, feat_cols = feat_cols_pos, 
  confounders_str = confounders_str, 
  group_subjects = c("Control", "T1D"), doPar = TRUE)

merge_pvals_CT_DM1_pos_sex <- merge(
  pvals_CT_DM1_sex, 
  lipidsearch_pos[,c("Feature", "LipidIon", "Class", "NCarbons", "NBonds","Grade", 
                     "lipidsearch_mz_ppm", "XCMS_rt", "lipidsearch_rt", 
                     "diff_rt", "confirmed")], 
  all.x = TRUE)

write.csv(merge_pvals_CT_DM1_pos_sex, "Results/pvalues/CT_DM1_pos_sex.csv")

###################### Negative ######################

pvals_CT_DM1 <- linear_analysis(data = merge_neg, feat_cols = feat_cols_neg, 
                                 confounders_str = confounders_str, 
                                 group_subjects = c("Control", "T1D"), 
                                 doPar = TRUE)

merge_pvals_CT_DM1_neg <- merge(pvals_CT_DM1, 
                                 lipidsearch_neg[,c("Feature", "LipidIon", "Class", 
                                                    "NCarbons", "NBonds","Grade", 
                                                    "lipidsearch_mz_ppm", "XCMS_rt", 
                                                    "lipidsearch_rt", "diff_rt", 
                                                    "confirmed")], 
                                 all.x = TRUE)

write.csv(merge_pvals_CT_DM1_neg, "Results/pvalues/CT_DM1_neg.csv")

pvals_CT_DM1_sex <- linear_analysis_interactions(
  data = merge_neg, feat_cols = feat_cols_neg, 
  confounders_str = confounders_str, 
  group_subjects = c("Control", "T1D"), doPar = TRUE)

merge_pvals_CT_DM1_neg_sex <- merge(
  pvals_CT_DM1_sex, 
  lipidsearch_neg[,c("Feature", "LipidIon", "Class", "NCarbons", "NBonds","Grade", 
                     "lipidsearch_mz_ppm", "XCMS_rt", "lipidsearch_rt", 
                     "diff_rt", "confirmed")], 
  all.x = TRUE)

write.csv(merge_pvals_CT_DM1_neg_sex, "Results/pvalues/CT_DM1_neg_sex.csv")

######################################################
##################### CT vs T2D ######################
######################################################

###################### Positive ######################

confounders_str <- paste(
  "Sex", "AGE", "Hypertension", "Dyslipidemia", "DM", "BMI", 
  "Glucose", "Smoking", "systolic_BP", "Triglycerides",
  "HDL_cholesterol", "LDL_cholesterol", "aMED", "MDRD_4_IDMS", sep = "+")

pvals_CT_DM2 <- linear_analysis(data = merge_pos, feat_cols = feat_cols_pos, 
                                confounders_str = confounders_str, 
                                group_subjects = c("Control", "T2D"), 
                                doPar = TRUE)

merge_pvals_CT_DM2_pos <- merge(pvals_CT_DM2, 
                                lipidsearch_pos[,c("Feature", "LipidIon", "Class", 
                                                   "NCarbons", "NBonds","Grade", 
                                                   "lipidsearch_mz_ppm", "XCMS_rt", 
                                                   "lipidsearch_rt", "diff_rt", 
                                                   "confirmed")], 
                                all.x = TRUE)

write.csv(merge_pvals_CT_DM2_pos, "Results/pvalues/CT_DM2_pos.csv")

pvals_CT_DM2_sex <- linear_analysis_interactions(
  data = merge_pos, feat_cols = feat_cols_pos, 
  confounders_str = confounders_str, 
  group_subjects = c("Control", "T2D"), doPar = TRUE)

merge_pvals_CT_DM2_pos_sex <- merge(
  pvals_CT_DM2_sex, 
  lipidsearch_pos[,c("Feature", "LipidIon", "Class", "NCarbons", "NBonds","Grade", 
                     "lipidsearch_mz_ppm", "XCMS_rt", "lipidsearch_rt", 
                     "diff_rt", "confirmed")], 
  all.x = TRUE)

write.csv(merge_pvals_CT_DM2_pos_sex, "Results/pvalues/CT_DM2_pos_sex.csv")

###################### Negative ######################

pvals_CT_DM2 <- linear_analysis(data = merge_neg, feat_cols = feat_cols_neg, 
                                confounders_str = confounders_str, 
                                group_subjects = c("Control", "T2D"), 
                                doPar = TRUE)

merge_pvals_CT_DM2_neg <- merge(pvals_CT_DM2, 
                                lipidsearch_neg[,c("Feature", "LipidIon", "Class", 
                                                   "NCarbons", "NBonds","Grade", 
                                                   "lipidsearch_mz_ppm", "XCMS_rt", 
                                                   "lipidsearch_rt", "diff_rt", 
                                                   "confirmed")], 
                                all.x = TRUE)

write.csv(merge_pvals_CT_DM2_neg, "Results/pvalues/CT_DM2_neg.csv")

pvals_CT_DM2_sex <- linear_analysis_interactions(
  data = merge_neg, feat_cols = feat_cols_neg, 
  confounders_str = confounders_str, 
  group_subjects = c("Control", "T2D"), doPar = TRUE)

merge_pvals_CT_DM2_neg_sex <- merge(
  pvals_CT_DM2_sex, 
  lipidsearch_neg[,c("Feature", "LipidIon", "Class", "NCarbons", "NBonds","Grade", 
                     "lipidsearch_mz_ppm", "XCMS_rt", "lipidsearch_rt", 
                     "diff_rt", "confirmed")], 
  all.x = TRUE)

write.csv(merge_pvals_CT_DM2_neg_sex, "Results/pvalues/CT_DM2_neg_sex.csv")







