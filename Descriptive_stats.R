####### Packages loading #######
library(ggplot2)
library(magrittr)
library(plyr)
library(dplyr)
library(knitr)
library(viridis)
doParallel::registerDoParallel(6)

dir_path <- "/home/maria/paper_placa/"

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

####### Clinical variables #######
clinical_vars <- read.csv(file = "data/clinical_variables.csv", 
                          row.names = 1, stringsAsFactors = TRUE)

numeric_data <- clinical_vars

for (j in c("Sex", "Hypertension", "Dyslipidemia", "Plaque", 
            "Insulin", "Metformin", "Antiplatelet")){
  levels(numeric_data[,j]) <- c("0", "1")
}

levels(numeric_data$DM) <- c("0", "1", "2")
levels(numeric_data$Origen) <- c("0", "1", "2", "3")
levels(numeric_data$Smoking) <- c("1", "0")

for (j in c("Sex", "Hypertension", "Dyslipidemia", "Plaque", 
            "Insulin", "Metformin", "Antiplatelet",
            "DM", "Origen", "Smoking")){
  numeric_data[,j] <- as.numeric(as.character(numeric_data[,j]))
}

corrplot::corrplot(cor(numeric_data, use = "pairwise.complete.obs"), tl.col = 'black', tl.cex = 0.5)

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

### SI Tables
res1 <- 
  compareGroups::compareGroups(DM ~ Sex+AGE+Hypertension+Dyslipidemia+BMI+
                                 Waist+DM_duration+HbA1c+Glucose+Smoking+
                                 systolic_BP+diastolic_BP+Triglycerides+HDL_cholesterol+
                                 LDL_cholesterol+aMED+MDRD_4_IDMS+
                                 Insulin+Metformin,
                               data = nonScaled_Pos[nonScaled_Pos$DM %in% c("T1D", "T2D"),], 
                               simplify = TRUE,
                               ref = 1)
compTable <- compareGroups::createTable(x = res1)

res2 <- 
  compareGroups::compareGroups(DM ~ Sex+AGE+Hypertension+Dyslipidemia+BMI+
                                 Waist+DM_duration+HbA1c+Glucose+Smoking+
                                 systolic_BP+diastolic_BP+Triglycerides+HDL_cholesterol+
                                 LDL_cholesterol+aMED+MDRD_4_IDMS+
                                 Insulin+Metformin,
                               data = nonScaled_Pos[nonScaled_Pos$DM %in% c("T1D", "Control"),], 
                               simplify = TRUE,
                               ref = 1)
compTable2 <- compareGroups::createTable(x = res2)

res3 <- 
  compareGroups::compareGroups(DM ~ Sex+AGE+Hypertension+Dyslipidemia+BMI+
                                 Waist+DM_duration+HbA1c+Glucose+Smoking+
                                 systolic_BP+diastolic_BP+Triglycerides+HDL_cholesterol+
                                 LDL_cholesterol+aMED+MDRD_4_IDMS+
                                 Insulin+Metformin,
                               data = nonScaled_Pos[nonScaled_Pos$DM %in% c("T2D", "Control"),], 
                               simplify = TRUE,
                               ref = 1)
compTable3 <- compareGroups::createTable(x = res3)

compareGroups::export2csv(compTable, 
                          file="Results/SI_main_analysis_t1d_t2d.csv")

compareGroups::export2csv(compTable2, 
                          file="Results/SI_main_analysis_t1d_ct.csv")

compareGroups::export2csv(compTable3, 
                          file="Results/SI_main_analysis_t2d_ct.csv")

### Manuscript tables - data description
res_total <- 
  compareGroups::compareGroups(DM ~ Sex+AGE+Hypertension+Dyslipidemia+BMI+
                                 Waist+DM_duration+HbA1c+Glucose+Smoking+
                                 systolic_BP+diastolic_BP+Triglycerides+HDL_cholesterol+
                                 LDL_cholesterol+
                                 Insulin+Metformin+aMED+MDRD_4_IDMS,
                               data = nonScaled_Pos, 
                               simplify = TRUE,
                               ref = 1)

compTable_total <- compareGroups::createTable(x = res_total)
compareGroups::export2csv(compTable_total, 
                          file="Results/general_data_descr.csv")

nonScaled_Pos <- merge_pos
nonScaled_Neg <- merge_neg


kableExtra::kable(data.frame(Pos = nrow(merge_pos), Neg = nrow(merge_neg)), 
                  col.names = c("Positive ionization mode", 
                                "Negative ionization mode"), 
                  caption = "Number of patients in each acquisition mode") %>%
  kableExtra::kable_styling(latex_options = "HOLD_position")
n.feats.before.pos <- ncol(Pos_table)
n.feats.before.neg <- ncol(Neg_table)

# correlation between diabetes and glucose/HbA1c
diabetesT1D <- merge_pos$DM[merge_pos$DM %in% c("Control", "T1D")]
diabetesT1D <- as.factor(as.character(diabetesT1D))
levels(diabetesT1D) <- c("0", "1")
cor(as.numeric(as.character(diabetesT1D)), 
    merge_pos$HbA1c[merge_pos$DM %in% c("Control", "T1D")], 
    use = "pairwise.complete.obs")
cor(as.numeric(as.character(diabetesT1D)), 
    merge_pos$Glucose[merge_pos$DM %in% c("Control", "T1D")], 
    use = "pairwise.complete.obs")

diabetesT2D <- merge_pos$DM[merge_pos$DM %in% c("Control", "T2D")]
diabetesT2D <- as.factor(as.character(diabetesT2D))
levels(diabetesT2D) <- c("0", "1")
cor(as.numeric(as.character(diabetesT2D)), 
    merge_pos$HbA1c[merge_pos$DM %in% c("Control", "T2D")], 
    use = "pairwise.complete.obs")
cor(as.numeric(as.character(diabetesT2D)), 
    merge_pos$Glucose[merge_pos$DM %in% c("Control", "T2D")], 
    use = "pairwise.complete.obs")

library(ggfortify)
pca_data <- merge_pos
pca_data[is.na(pca_data)] <- 0
p1 <- autoplot(prcomp(pca_data[pca_data$DM %in% c("T1D", "T2D"),-c(1:24)], center = TRUE, scale = TRUE), 
               colour = "DM", data = pca_data[pca_data$DM %in% c("T1D", "T2D"),],
               frame = TRUE, frame.type = 'norm', scale = FALSE)

p2 <- autoplot(prcomp(pca_data[pca_data$DM %in% c("Control", "T2D"),-c(1:24)], center = TRUE, scale = TRUE), 
               colour = "DM", data = pca_data[pca_data$DM %in% c("Control", "T2D"),],
               frame = TRUE, frame.type = 'norm', scale = FALSE)

p3 <- autoplot(prcomp(pca_data[pca_data$DM %in% c("Control", "T1D"),-c(1:24)], center = TRUE, scale = TRUE), 
               colour = "DM", data = pca_data[pca_data$DM %in% c("Control", "T1D"),],
               frame = TRUE, frame.type = 'norm', scale = FALSE)

pca_data <- merge_neg
pca_data[is.na(pca_data)] <- 0
p1 <- autoplot(prcomp(pca_data[pca_data$DM %in% c("T1D", "T2D"),-c(1:24)], center = TRUE, scale = TRUE), 
               colour = "DM", data = pca_data[pca_data$DM %in% c("T1D", "T2D"),],
               frame = TRUE, frame.type = 'norm', scale = FALSE)

p2 <- autoplot(prcomp(pca_data[pca_data$DM %in% c("Control", "T2D"),-c(1:24)], center = TRUE, scale = TRUE), 
               colour = "DM", data = pca_data[pca_data$DM %in% c("Control", "T2D"),],
               frame = TRUE, frame.type = 'norm', scale = FALSE)

p3 <- autoplot(prcomp(pca_data[pca_data$DM %in% c("Control", "T1D"),-c(1:24)], center = TRUE, scale = TRUE), 
               colour = "DM", data = pca_data[pca_data$DM %in% c("Control", "T1D"),],
               frame = TRUE, frame.type = 'norm', scale = FALSE)

res1 <- 
  compareGroups::compareGroups(DM ~ Sex+AGE+Hypertension+Dyslipidemia+BMI+
                                 Waist+DM_duration+HbA1c+Glucose+Smoking+
                                 systolic_BP+diastolic_BP+Triglycerides+HDL_cholesterol+
                                 LDL_cholesterol+aMED+Insulin+Metformin+MDRD_4_IDMS,
                               data = nonScaled_Pos[nonScaled_Pos$DM %in% c("T1D", "T2D"),], 
                               simplify = TRUE,
                               ref = 1)

compTable <- compareGroups::createTable(x = res1)

dda_t1d_ct <- nonScaled_Pos[nonScaled_Pos$DM %in% c("T1D", "Control"),]
dda_t1d_ct$DM <- as.factor(as.character(dda_t1d_ct$DM))
res2 <- 
  compareGroups::compareGroups(DM ~ Sex+AGE+Hypertension+Dyslipidemia+BMI+
                                 Waist+DM_duration+HbA1c+Glucose+Smoking+
                                 systolic_BP+diastolic_BP+Triglycerides+HDL_cholesterol+
                                 LDL_cholesterol+aMED+Insulin+Metformin+MDRD_4_IDMS,
                               data = dda_t1d_ct, 
                               simplify = FALSE,
                               ref = 1)

compTable2 <- compareGroups::createTable(x = res2)

res3 <- 
  compareGroups::compareGroups(DM ~ Sex+AGE+Hypertension+Dyslipidemia+BMI+
                                 Waist+DM_duration+HbA1c+Glucose+Smoking+
                                 systolic_BP+diastolic_BP+Triglycerides+HDL_cholesterol+
                                 LDL_cholesterol+aMED+Insulin+Metformin+MDRD_4_IDMS,
                               data = nonScaled_Pos[nonScaled_Pos$DM %in% c("T2D", "Control"),], 
                               simplify = TRUE,
                               ref = 1)

compTable3 <- compareGroups::createTable(x = res3)

ctable <- cbind("T1D vs T2D"=compTable,"T1D vs CT"=compTable2,"T2D vs CT"=compTable3)







