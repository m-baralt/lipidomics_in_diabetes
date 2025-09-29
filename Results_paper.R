library(plyr)
library(ggplot2)
library(ggpubr)
library(viridis)
library(ggnewscale)
library(ggpattern)
library(ggtext)
doParallel::registerDoParallel(10)
source("Required_functions.R")

names_files <- list.files("Results/pvalues/")
for (n in names_files){
  assign(paste("all_",gsub(n, pattern = ".csv", replacement = ""), sep = ""), 
         read.csv(paste("Results/pvalues/", n, sep = ""), row.names = 1))
}

dir_path <- "/home/maria/paper_placa/"
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

### T1D vs T2D ###

R11_pos <- sign_features(pvals_data = all_DM1_DM2_pos, variable_pval = "p_adj", 
                         analysis_name = "R11", ionization = "Positive", sign_pval = 0.05)

R11_neg <- sign_features(pvals_data = all_DM1_DM2_neg, variable_pval = "p_adj", 
                         analysis_name = "R11", ionization = "Negative", sign_pval = 0.05)

R12_pos <- sign_features(pvals_data = all_DM1_DM2_pos_sex, variable_pval = "p_adj_women", 
                         analysis_name = "R12", ionization = "Positive", sign_pval = 0.05)

R12_neg <- sign_features(pvals_data = all_DM1_DM2_neg_sex, variable_pval = "p_adj_women", 
                         analysis_name = "R12", ionization = "Negative", sign_pval = 0.05)

R13_pos <- sign_features(pvals_data = all_DM1_DM2_pos_sex, variable_pval = "p_adj_men", 
                         analysis_name = "R13", ionization = "Positive", sign_pval = 0.05)

R13_neg <- sign_features(pvals_data = all_DM1_DM2_neg_sex, variable_pval = "p_adj_men", 
                         analysis_name = "R13", ionization = "Negative", sign_pval = 0.05)

### CT vs T1D ###

R21_pos <- sign_features(pvals_data = all_CT_DM1_pos, variable_pval = "p_adj", 
                         analysis_name = "R21", ionization = "Positive", sign_pval = 0.05)

R21_neg <- sign_features(pvals_data = all_CT_DM1_neg, variable_pval = "p_adj", 
                         analysis_name = "R21", ionization = "Negative", sign_pval = 0.05)

R22_pos <- sign_features(pvals_data = all_CT_DM1_pos_sex, variable_pval = "p_adj_women", 
                         analysis_name = "R22", ionization = "Positive", sign_pval = 0.05)

R22_neg <- sign_features(pvals_data = all_CT_DM1_neg_sex, variable_pval = "p_adj_women", 
                         analysis_name = "R22", ionization = "Negative", sign_pval = 0.05)

R23_pos <- sign_features(pvals_data = all_CT_DM1_pos_sex, variable_pval = "p_adj_men", 
                         analysis_name = "R23", ionization = "Positive", sign_pval = 0.05)

R23_neg <- sign_features(pvals_data = all_CT_DM1_neg_sex, variable_pval = "p_adj_men", 
                         analysis_name = "R23", ionization = "Negative", sign_pval = 0.05)

### CT vs T2D ###

R31_pos <- sign_features(pvals_data = all_CT_DM2_pos, variable_pval = "p_adj", 
                         analysis_name = "R31", ionization = "Positive", sign_pval = 0.05)

R31_neg <- sign_features(pvals_data = all_CT_DM2_neg, variable_pval = "p_adj", 
                         analysis_name = "R31", ionization = "Negative", sign_pval = 0.05)

R32_pos <- sign_features(pvals_data = all_CT_DM2_pos_sex, variable_pval = "p_adj_women", 
                         analysis_name = "R32", ionization = "Positive", sign_pval = 0.05)

R32_neg <- sign_features(pvals_data = all_CT_DM2_neg_sex, variable_pval = "p_adj_women", 
                         analysis_name = "R32", ionization = "Negative", sign_pval = 0.05)

R33_pos <- sign_features(pvals_data = all_CT_DM2_pos_sex, variable_pval = "p_adj_men", 
                         analysis_name = "R33", ionization = "Positive", sign_pval = 0.05)

R33_neg <- all_CT_DM2_neg_sex[(all_CT_DM2_neg_sex$p_adj_men<0.05) & (!(is.na(all_CT_DM2_neg_sex$LipidIon))),]
R33_neg <- R33_neg[R33_neg$confirmed==2,]


sum(nrow(R11_pos), nrow(R11_neg), nrow(R12_pos), nrow(R12_neg), nrow(R13_pos), nrow(R13_neg),
    nrow(R21_pos), nrow(R21_neg), nrow(R22_pos), nrow(R22_neg), nrow(R23_pos), nrow(R23_neg),
    nrow(R31_pos), nrow(R31_neg), nrow(R32_pos), nrow(R32_neg), nrow(R33_pos), nrow(R33_neg))

Rtotal <- rbind(R11_pos, R11_neg, R12_pos, R12_neg, R13_pos, R13_neg,
                R21_pos, R21_neg, R22_pos, R22_neg, R23_pos, R23_neg,
                R31_pos, R31_neg, R32_pos, R32_neg, R33_pos, R33_neg)

lipidsearch_pos_unique <- lipidsearch_pos[!duplicated(lipidsearch_pos$Feature),]
RtotalPos <- ddply(Rtotal[Rtotal$Ionization == "Positive",], "Analysis", function(dx){
  dy <- merge(dx, 
              lipidsearch_pos_unique[,c("Feature", "XCMS_mz","XCMS_rt")], 
              by = "Feature", all.x = TRUE)
  dy
})

lipidsearch_neg_unique <- lipidsearch_neg[!duplicated(lipidsearch_neg$Feature),]
RtotalNeg <- ddply(Rtotal[Rtotal$Ionization == "Negative",], "Analysis", function(dx){
  dy <- merge(dx, 
              lipidsearch_neg_unique[,c("Feature", "XCMS_mz","XCMS_rt")], 
              by = "Feature", all.x = TRUE)
  dy
})

Rtotal <- rbind(RtotalPos, RtotalNeg)

Rtotal.final <- ddply(Rtotal, ~Feature+Ionization, function(dy){
  if (nrow(dy)>1){
    dx <- dy[1,]
    dx$Analysis <- paste(dy$Analysis, collapse = ",")
    dx$`q-value` <- paste("(",formatC(min(dy$`q-value`), format = "e", digits = 2),
                          ",",formatC(max(dy$`q-value`), format = "e", digits = 2),")", sep = "")
    return(dx)
  } else {
    dy$`q-value` <- formatC(min(dy$`q-value`), format = "e", digits = 2)
    return(dy)
  }
})

colnames(Rtotal.final)[c(5,7,8)] <- c("corrected p-value","mz", "rt")
Rtotal.final <- Rtotal.final[,c(1,3,2,7,8,5,6,4)]
Rtotal.final$mz <- round(Rtotal.final$mz, digits = 4)
Rtotal.final$rt <- round(Rtotal.final$rt, digits = 2)

openxlsx::write.xlsx(Rtotal.final[,-1], "Results/TableS3.xlsx")

##### Number of features significant, number of features sign and id, number of features sign and id with confirmed = 2

sign_pval <- 0.05

## Positive
# DM1 vs DM2
all_DM1_DM2_pos_sign <- all_DM1_DM2_pos[all_DM1_DM2_pos$p_adj<sign_pval,]
sign.feats.pos <- all_DM1_DM2_pos_sign$Feature
length(unique(all_DM1_DM2_pos_sign$Feature))
all_DM1_DM2_pos_id <- all_DM1_DM2_pos_sign[!is.na(all_DM1_DM2_pos_sign$LipidIon),]
length(unique(all_DM1_DM2_pos_id$Feature))
table(all_DM1_DM2_pos_id$confirmed[!duplicated(all_DM1_DM2_pos_id$Feature)])

#negative
all_DM1_DM2_neg_sign <- all_DM1_DM2_neg[all_DM1_DM2_neg$p_adj<sign_pval,]
sign.feats.neg <- all_DM1_DM2_neg_sign$Feature
length(unique(all_DM1_DM2_neg_sign$Feature))
all_DM1_DM2_neg_id <- all_DM1_DM2_neg_sign[!is.na(all_DM1_DM2_neg_sign$LipidIon),]
length(unique(all_DM1_DM2_neg_id$Feature))
table(all_DM1_DM2_neg_id$confirmed[!duplicated(all_DM1_DM2_neg_id$Feature)])

# DM1 vs DM2 women
all_DM1_DM2_pos_sign_women <- all_DM1_DM2_pos_sex[all_DM1_DM2_pos_sex$p_adj_women<sign_pval,]
sign.feats.pos <- c(sign.feats.pos,all_DM1_DM2_pos_sign_women$Feature)
length(unique(all_DM1_DM2_pos_sign_women$Feature))
all_DM1_DM2_pos_id_women <- all_DM1_DM2_pos_sign_women[!is.na(all_DM1_DM2_pos_sign_women$LipidIon),]
length(unique(all_DM1_DM2_pos_id_women$Feature))
table(all_DM1_DM2_pos_id_women$confirmed[!duplicated(all_DM1_DM2_pos_id_women$Feature)])

#negative
all_DM1_DM2_neg_sign_women <- all_DM1_DM2_neg_sex[all_DM1_DM2_neg_sex$p_adj_women<sign_pval,]
sign.feats.neg <- c(sign.feats.neg,all_DM1_DM2_neg_sign_women$Feature)
length(unique(all_DM1_DM2_neg_sign_women$Feature))
all_DM1_DM2_neg_id_women <- all_DM1_DM2_neg_sign_women[!is.na(all_DM1_DM2_neg_sign_women$LipidIon),]
length(unique(all_DM1_DM2_neg_id_women$Feature))
table(all_DM1_DM2_neg_id_women$confirmed[!duplicated(all_DM1_DM2_neg_id_women$Feature)])

# DM1 vs DM2 men
all_DM1_DM2_pos_sign_men <- all_DM1_DM2_pos_sex[all_DM1_DM2_pos_sex$p_adj_men<sign_pval,]
sign.feats.pos <- c(sign.feats.pos,all_DM1_DM2_pos_sign_men$Feature)
length(unique(all_DM1_DM2_pos_sign_men$Feature))
all_DM1_DM2_pos_id_men <- all_DM1_DM2_pos_sign_men[!is.na(all_DM1_DM2_pos_sign_men$LipidIon),]
length(unique(all_DM1_DM2_pos_id_men$Feature))
table(all_DM1_DM2_pos_id_men$confirmed[!duplicated(all_DM1_DM2_pos_id_men$Feature)])

# negative
all_DM1_DM2_neg_sign_men <- all_DM1_DM2_neg_sex[all_DM1_DM2_neg_sex$p_adj_men<sign_pval,]
sign.feats.neg <- c(sign.feats.neg,all_DM1_DM2_neg_sign_men$Feature)
length(unique(all_DM1_DM2_neg_sign_men$Feature))
all_DM1_DM2_neg_id_men <- all_DM1_DM2_neg_sign_men[!is.na(all_DM1_DM2_neg_sign_men$LipidIon),]
length(unique(all_DM1_DM2_neg_id_men$Feature))
table(all_DM1_DM2_neg_id_men$confirmed[!duplicated(all_DM1_DM2_neg_id_men$Feature)])

## Positive
# DM1 vs CT
all_CT_DM1_pos_sign <- all_CT_DM1_pos[all_CT_DM1_pos$p_adj<sign_pval,]
sign.feats.pos <- c(sign.feats.pos,all_CT_DM1_pos_sign$Feature)
length(unique(all_CT_DM1_pos_sign$Feature))
all_CT_DM1_pos_id <- all_CT_DM1_pos_sign[!is.na(all_CT_DM1_pos_sign$LipidIon),]
length(unique(all_CT_DM1_pos_id$Feature))
table(all_CT_DM1_pos_id$confirmed[!duplicated(all_CT_DM1_pos_id$Feature)])

#negative
all_CT_DM1_neg_sign <- all_CT_DM1_neg[all_CT_DM1_neg$p_adj<sign_pval,]
sign.feats.neg <- c(sign.feats.neg,all_CT_DM1_neg_sign$Feature)
length(unique(all_CT_DM1_neg_sign$Feature))
all_CT_DM1_neg_id <- all_CT_DM1_neg_sign[!is.na(all_CT_DM1_neg_sign$LipidIon),]
length(unique(all_CT_DM1_neg_id$Feature))
table(all_CT_DM1_neg_id$confirmed[!duplicated(all_CT_DM1_neg_id$Feature)])

# DM1 vs CT women
all_CT_DM1_pos_sign_women <- all_CT_DM1_pos_sex[all_CT_DM1_pos_sex$p_adj_women<sign_pval,]
sign.feats.pos <- c(sign.feats.pos,all_CT_DM1_pos_sign_women$Feature)
length(unique(all_CT_DM1_pos_sign_women$Feature))
all_CT_DM1_pos_id_women <- all_CT_DM1_pos_sign_women[!is.na(all_CT_DM1_pos_sign_women$LipidIon),]
length(unique(all_CT_DM1_pos_id_women$Feature))
table(all_CT_DM1_pos_id_women$confirmed[!duplicated(all_CT_DM1_pos_id_women$Feature)])

# negative
all_CT_DM1_neg_sign_women <- all_CT_DM1_neg_sex[all_CT_DM1_neg_sex$p_adj_women<sign_pval,]
sign.feats.neg <- c(sign.feats.neg,all_CT_DM1_neg_sign_women$Feature)
length(unique(all_CT_DM1_neg_sign_women$Feature))
all_CT_DM1_neg_id_women <- all_CT_DM1_neg_sign_women[!is.na(all_CT_DM1_neg_sign_women$LipidIon),]
length(unique(all_CT_DM1_neg_id_women$Feature))
table(all_CT_DM1_neg_id_women$confirmed[!duplicated(all_CT_DM1_neg_id_women$Feature)])

# DM1 vs CT men
all_CT_DM1_pos_sign_men <- all_CT_DM1_pos_sex[all_CT_DM1_pos_sex$p_adj_men<sign_pval,]
sign.feats.pos <- c(sign.feats.pos,all_CT_DM1_pos_sign_men$Feature)
length(unique(all_CT_DM1_pos_sign_men$Feature))
all_CT_DM1_pos_id_men <- all_CT_DM1_pos_sign_men[!is.na(all_CT_DM1_pos_sign_men$LipidIon),]
length(unique(all_CT_DM1_pos_id_men$Feature))
table(all_CT_DM1_pos_id_men$confirmed[!duplicated(all_CT_DM1_pos_id_men$Feature)])

# negative
all_CT_DM1_neg_sign_men <- all_CT_DM1_neg_sex[all_CT_DM1_neg_sex$p_adj_men<sign_pval,]
sign.feats.neg <- c(sign.feats.neg,all_CT_DM1_neg_sign_men$Feature)
length(unique(all_CT_DM1_neg_sign_men$Feature))
all_CT_DM1_neg_id_men <- all_CT_DM1_neg_sign_men[!is.na(all_CT_DM1_neg_sign_men$LipidIon),]
length(unique(all_CT_DM1_neg_id_men$Feature))
table(all_CT_DM1_neg_id_men$confirmed[!duplicated(all_CT_DM1_neg_id_men$Feature)])

## Positive
# DM2 vs CT
all_CT_DM2_pos_sign <- all_CT_DM2_pos[all_CT_DM2_pos$p_adj<sign_pval,]
sign.feats.pos <- c(sign.feats.pos,all_CT_DM2_pos_sign$Feature)
length(unique(all_CT_DM2_pos_sign$Feature))
all_CT_DM2_pos_id <- all_CT_DM2_pos_sign[!is.na(all_CT_DM2_pos_sign$LipidIon),]
length(unique(all_CT_DM2_pos_id$Feature))
table(all_CT_DM2_pos_id$confirmed[!duplicated(all_CT_DM2_pos_id$Feature)])

# negative
all_CT_DM2_neg_sign <- all_CT_DM2_neg[all_CT_DM2_neg$p_adj<sign_pval,]
sign.feats.neg <- c(sign.feats.neg,all_CT_DM2_neg_sign$Feature)
length(unique(all_CT_DM2_neg_sign$Feature))
all_CT_DM2_neg_id <- all_CT_DM2_neg_sign[!is.na(all_CT_DM2_neg_sign$LipidIon),]
length(unique(all_CT_DM2_neg_id$Feature))
table(all_CT_DM2_neg_id$confirmed[!duplicated(all_CT_DM2_neg_id$Feature)])

# DM2 vs CT women
all_CT_DM2_pos_sign_women <- all_CT_DM2_pos_sex[all_CT_DM2_pos_sex$p_adj_women<sign_pval,]
sign.feats.pos <- c(sign.feats.pos,all_CT_DM2_pos_sign_women$Feature)
length(unique(all_CT_DM2_pos_sign_women$Feature))
all_CT_DM2_pos_id_women <- all_CT_DM2_pos_sign_women[!is.na(all_CT_DM2_pos_sign_women$LipidIon),]
length(unique(all_CT_DM2_pos_id_women$Feature))
table(all_CT_DM2_pos_id_women$confirmed[!duplicated(all_CT_DM2_pos_id_women$Feature)])

# negative
all_CT_DM2_neg_sign_women <- all_CT_DM2_neg_sex[all_CT_DM2_neg_sex$p_adj_women<sign_pval,]
sign.feats.neg <- c(sign.feats.neg,all_CT_DM2_neg_sign_women$Feature)
length(unique(all_CT_DM2_neg_sign_women$Feature))
all_CT_DM2_neg_id_women <- all_CT_DM2_neg_sign_women[!is.na(all_CT_DM2_neg_sign_women$LipidIon),]
length(unique(all_CT_DM2_neg_id_women$Feature))
table(all_CT_DM2_neg_id_women$confirmed[!duplicated(all_CT_DM2_neg_id_women$Feature)])

# DM2 vs CT men
all_CT_DM2_pos_sign_men <- all_CT_DM2_pos_sex[all_CT_DM2_pos_sex$p_adj_men<sign_pval,]
sign.feats.pos <- c(sign.feats.pos,all_CT_DM2_pos_sign_men$Feature)
length(unique(all_CT_DM2_pos_sign_men$Feature))
all_CT_DM2_pos_id_men <- all_CT_DM2_pos_sign_men[!is.na(all_CT_DM2_pos_sign_men$LipidIon),]
length(unique(all_CT_DM2_pos_id_men$Feature))
table(all_CT_DM2_pos_id_men$confirmed[!duplicated(all_CT_DM2_pos_id_men$Feature)])

# negative
all_CT_DM2_neg_sign_men <- all_CT_DM2_neg_sex[all_CT_DM2_neg_sex$p_adj_men<sign_pval,]
sign.feats.neg <- c(sign.feats.neg,all_CT_DM2_neg_sign_men$Feature)
length(unique(all_CT_DM2_neg_sign_men$Feature))
all_CT_DM2_neg_id_men <- all_CT_DM2_neg_sign_men[!is.na(all_CT_DM2_neg_sign_men$LipidIon),]
length(unique(all_CT_DM2_neg_id_men$Feature))
table(all_CT_DM2_neg_id_men$confirmed[!duplicated(all_CT_DM2_neg_id_men$Feature)])

##### Upset plots of significant features
df <- data.frame(Features = unique(c(as.character(all_DM1_DM2_pos$Feature[all_DM1_DM2_pos$p_adj<sign_pval]),
                                     as.character(all_DM1_DM2_pos_sex$Feature[all_DM1_DM2_pos_sex$p_adj_women<sign_pval]),
                                     as.character(all_DM1_DM2_pos_sex$Feature[all_DM1_DM2_pos_sex$p_adj_men<sign_pval]),
                                     as.character(all_CT_DM1_pos$Feature[all_CT_DM1_pos$p_adj<sign_pval]),
                                     as.character(all_CT_DM1_pos_sex$Feature[all_CT_DM1_pos_sex$p_adj_women<sign_pval]),
                                     as.character(all_CT_DM1_pos_sex$Feature[all_CT_DM1_pos_sex$p_adj_men<sign_pval]),
                                     as.character(all_CT_DM2_pos$Feature[all_CT_DM2_pos$p_adj<sign_pval]),
                                     as.character(all_CT_DM2_pos_sex$Feature[all_CT_DM2_pos_sex$p_adj_women<sign_pval]),
                                     as.character(all_CT_DM2_pos_sex$Feature[all_CT_DM2_pos_sex$p_adj_men<sign_pval]))))

df$R11 <- as.numeric(as.character(df$Features) %in% as.character(all_DM1_DM2_pos$Feature[all_DM1_DM2_pos$p_adj<sign_pval]))
df$R12 <- as.numeric(as.character(df$Features) %in% as.character(all_DM1_DM2_pos_sex$Feature[all_DM1_DM2_pos_sex$p_adj_women<sign_pval]))
df$R13 <- as.numeric(as.character(df$Features) %in% as.character(all_DM1_DM2_pos_sex$Feature[all_DM1_DM2_pos_sex$p_adj_men<sign_pval]))
df$R21 <- as.numeric(as.character(df$Features) %in% as.character(all_CT_DM1_pos$Feature[all_CT_DM1_pos$p_adj<sign_pval]))
df$R22 <- as.numeric(as.character(df$Features) %in% as.character(all_CT_DM1_pos_sex$Feature[all_CT_DM1_pos_sex$p_adj_women<sign_pval]))
df$R23 <- as.numeric(as.character(df$Features) %in% as.character(all_CT_DM1_pos_sex$Feature[all_CT_DM1_pos_sex$p_adj_men<sign_pval]))
df$R31 <- as.numeric(as.character(df$Features) %in% as.character(all_CT_DM2_pos$Feature[all_CT_DM2_pos$p_adj<sign_pval]))
df$R32 <- as.numeric(as.character(df$Features) %in% as.character(all_CT_DM2_pos_sex$Feature[all_CT_DM2_pos_sex$p_adj_women<sign_pval]))
df$R33 <- as.numeric(as.character(df$Features) %in% as.character(all_CT_DM2_pos_sex$Feature[all_CT_DM2_pos_sex$p_adj_men<sign_pval]))

p <- UpSetR::upset(df, sets = c(
  "R33", "R32", "R31", 
  "R23", "R22", "R21", 
  "R13", "R12", "R11"), 
  keep.order = TRUE, order.by = "freq", 
  mb.ratio = c(0.6,0.4), nintersects = NA,
  text.scale = c(1.2, 1.2,1.2,1.4,1.4,1.4),
  mainbar.y.label = "Coincident features \n",
  sets.x.label = "Significant features",
  point.size = 2, #scale.sets = "log10",
  sets.bar.color = c(
    "darkorchid1", "darkorchid1", "darkorchid1",
    "goldenrod1", "goldenrod1", "goldenrod1",
    "turquoise3", "turquoise3", "turquoise3"),set_size.numbers_size = 7,
  set_size.show = TRUE, set_size.scale_max = 900)

newVars <- c("Common variables","DM duration & HbA1c", "DM:Sex")
vars.mat <- data.frame(Code = c("R11", "R12", "R13", 
                                "R21", "R22", "R23", 
                                "R31", "R32", "R33"), 
                       V1 = c(1,1,1,1,1,1,1,1,1), 
                       V2 = c(1,1,1,0,0,0,0,0,0), 
                       V3 = c(0,1,1,0,1,1,0,1,1))
colnames(vars.mat)[-1] <- newVars
k <- reshape2::melt(vars.mat)
k$value <- as.factor(as.character(k$value))
p2 <- ggplot(aes(y = Code, x = variable), data = k) + 
  geom_point(size = 2, aes(colour = value)) + 
  theme_bw() + 
  theme(panel.grid.major.x = element_blank(), 
        panel.grid.major.y = element_line(size=10, color=c("gray96", "white")), 
        axis.title.x = element_blank(), axis.title.y = element_blank(), 
        axis.text.x = element_text(angle = 70, size = 12), 
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x.top = element_text(vjust = 0.15, hjust = 0)) + 
  scale_x_discrete(position = "top") + 
  scale_y_discrete(limits = rev(levels(k$Code))) +
  scale_color_manual(name = "value", values = c("0" = "gray75", "1" = "black")) + 
  guides(colour = FALSE)

upset.pos <- ggarrange(ggarrange(NULL, p$Sizes, nrow = 2, heights = c(0.59,0.41)), NULL, 
                       ggarrange(ggarrange(NULL, p$Main_bar, widths = c(0.005,0.995)), 
                                 NULL, p$Matrix, NULL, nrow = 4, heights = c(0.57, 0.015, 0.37, 0.035)), 
                       ggarrange(NULL, ggarrange(NULL,p2, widths = c(0.02,0.99)), NULL, nrow = 3, 
                                 heights = c(0.355,0.57,0.04)), 
                       NULL, 
                       ncol = 5, widths = c(0.12,0.007,0.90,0.07,0.01), align = "h")

ggsave("Results/Upset_pos.png", plot = upset.pos, width = 20, height = 8, dpi = 300)

# Negative
df <- data.frame(Features = unique(c(as.character(all_DM1_DM2_neg$Feature[all_DM1_DM2_neg$p_adj<sign_pval]),
                                     as.character(all_DM1_DM2_neg_sex$Feature[all_DM1_DM2_neg_sex$p_adj_women<sign_pval]),
                                     as.character(all_DM1_DM2_neg_sex$Feature[all_DM1_DM2_neg_sex$p_adj_men<sign_pval]),
                                     as.character(all_CT_DM1_neg$Feature[all_CT_DM1_neg$p_adj<sign_pval]),
                                     as.character(all_CT_DM1_neg_sex$Feature[all_CT_DM1_neg_sex$p_adj_women<sign_pval]),
                                     as.character(all_CT_DM1_neg_sex$Feature[all_CT_DM1_neg_sex$p_adj_men<sign_pval]),
                                     as.character(all_CT_DM2_neg$Feature[all_CT_DM2_neg$p_adj<sign_pval]),
                                     as.character(all_CT_DM2_neg_sex$Feature[all_CT_DM2_neg_sex$p_adj_women<sign_pval]),
                                     as.character(all_CT_DM2_neg_sex$Feature[all_CT_DM2_neg_sex$p_adj_men<sign_pval]))))

df$R11 <- as.numeric(as.character(df$Features) %in% as.character(all_DM1_DM2_neg$Feature[all_DM1_DM2_neg$p_adj<sign_pval]))
df$R12 <- as.numeric(as.character(df$Features) %in% as.character(all_DM1_DM2_neg_sex$Feature[all_DM1_DM2_neg_sex$p_adj_women<sign_pval]))
df$R13 <- as.numeric(as.character(df$Features) %in% as.character(all_DM1_DM2_neg_sex$Feature[all_DM1_DM2_neg_sex$p_adj_men<sign_pval]))
df$R21 <- as.numeric(as.character(df$Features) %in% as.character(all_CT_DM1_neg$Feature[all_CT_DM1_neg$p_adj<sign_pval]))
df$R22 <- as.numeric(as.character(df$Features) %in% as.character(all_CT_DM1_neg_sex$Feature[all_CT_DM1_neg_sex$p_adj_women<sign_pval]))
df$R23 <- as.numeric(as.character(df$Features) %in% as.character(all_CT_DM1_neg_sex$Feature[all_CT_DM1_neg_sex$p_adj_men<sign_pval]))
df$R31 <- as.numeric(as.character(df$Features) %in% as.character(all_CT_DM2_neg$Feature[all_CT_DM2_neg$p_adj<sign_pval]))
df$R32 <- as.numeric(as.character(df$Features) %in% as.character(all_CT_DM2_neg_sex$Feature[all_CT_DM2_neg_sex$p_adj_women<sign_pval]))
df$R33 <- as.numeric(as.character(df$Features) %in% as.character(all_CT_DM2_neg_sex$Feature[all_CT_DM2_neg_sex$p_adj_men<sign_pval]))

p <- UpSetR::upset(df, sets = c(
  "R33", "R32", "R31", 
  "R23", "R22", "R21", 
  "R13", "R12", "R11"), 
  keep.order = TRUE, order.by = "freq", 
  mb.ratio = c(0.6,0.4), nintersects = NA,
  text.scale = c(1.2, 1.2,1.2,1.4,1.4,1.7),
  mainbar.y.label = "Coincident features \n",
  sets.x.label = "Significant features",
  point.size = 2, #scale.sets = "log10",
  sets.bar.color = c(
    "darkorchid1", "darkorchid1", "darkorchid1",
    "goldenrod1", "goldenrod1", "goldenrod1",
    "turquoise3", "turquoise3", "turquoise3"),
  set_size.numbers_size = 7,
  set_size.show = TRUE, set_size.scale_max = 400)

newVars <- c("Common variables","DM duration & HbA1c", "DM:Sex")
vars.mat <- data.frame(Code = c("R11", "R12", "R13", 
                                "R21", "R22", "R23", 
                                "R31", "R32", "R33"), 
                       V1 = c(1,1,1,1,1,1,1,1,1), 
                       V2 = c(1,1,1,0,0,0,0,0,0), 
                       V3 = c(0,1,1,0,1,1,0,1,1))
colnames(vars.mat)[-1] <- newVars

k <- reshape2::melt(vars.mat)
k$value <- as.factor(as.character(k$value))
p2 <- ggplot(aes(y = Code, x = variable), data = k) + 
  geom_point(size = 2, aes(colour = value)) + 
  theme_bw() + 
  theme(panel.grid.major.x = element_blank(), 
        panel.grid.major.y = element_line(size=10, color=c("gray96", "white")), 
        axis.title.x = element_blank(), axis.title.y = element_blank(), 
        axis.text.x = element_text(angle = 70, size = 12), 
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x.top = element_text(vjust = 0.15, hjust = 0)) + 
  scale_x_discrete(position = "top") + 
  scale_y_discrete(limits = rev(levels(k$Code))) +
  scale_color_manual(name = "value", values = c("0" = "gray75", "1" = "black")) + 
  guides(colour = FALSE)

upset.neg <- ggarrange(ggarrange(NULL, p$Sizes, nrow = 2, heights = c(0.59,0.41)), NULL, 
                       ggarrange(ggarrange(NULL, p$Main_bar, widths = c(0.005,0.995)), 
                                 NULL, p$Matrix, NULL, nrow = 4, heights = c(0.57, 0.015, 0.37, 0.035)), 
                       ggarrange(NULL, ggarrange(NULL,p2, widths = c(0.02,0.99)), NULL, nrow = 3, 
                                 heights = c(0.355,0.56,0.04)), 
                       NULL, 
                       ncol = 5, widths = c(0.12,0.007,0.90,0.07,0.01), align = "h")

ggsave("Results/Upset_neg.png", plot = upset.neg, width = 20, height = 8)

##### Manhattan plots
data.dm1.dm2 <- rbind(all_DM1_DM2_pos, all_DM1_DM2_neg)
data.dm1.dm2 <- data.dm1.dm2[!is.na(data.dm1.dm2$LipidIon),]
data.dm1.dm2 <- data.dm1.dm2[data.dm1.dm2$confirmed==2,]

data.dm1.dm2$Class <- as.factor(as.character(data.dm1.dm2$Class))
p.dm1.dm2 <- ggplot(data.dm1.dm2, aes(x = Class, y = -log10(p_adj), color = Class)) + 
  geom_point(size=2, position = position_dodge2(w = 0.25),alpha = 0.7) +
  geom_hline(yintercept = -log10(sign_pval), color = "grey40", linetype = "dashed") +
  scale_size_continuous(range = c(0.5,3)) + 
  ylab("-log10 (Corrected p-value)") + ylim(0, 5) + 
  labs(title = "A: T1D vs T2D") +
  theme_bw()+
  theme(axis.title.x = element_blank(), 
        legend.position = "none",
        axis.text.x = element_text(hjust=0.95,vjust=0.2, size = 12, angle = 90), 
        axis.text.y = element_text(size = 12), axis.title.y = element_text(size = 12),
        plot.title = element_text(size = 15))+
  scale_color_manual(values = c("lavenderblush2", "peachpuff", "lavenderblush2", 
                                "wheat2", "skyblue", "wheat2", "wheat2", "wheat2",
                                "wheat2","wheat2", "lavenderblush2", "plum2", "sandybrown",
                                "pink1", "deepskyblue2", "slateblue1", "gray75", 
                                "wheat2", "lavenderblush2","palegreen3", "lavenderblush2"))

data.dm1.ct <- rbind(all_CT_DM1_pos, all_CT_DM1_neg)
data.dm1.ct <- data.dm1.ct[!is.na(data.dm1.ct$LipidIon),]
data.dm1.ct <- data.dm1.ct[data.dm1.ct$confirmed==2,]

data.dm1.ct$Class <- as.factor(as.character(data.dm1.ct$Class))
p.dm1.ct <- ggplot(data.dm1.ct, aes(x = Class, y = -log10(p_adj), color = Class)) + 
  geom_point(size=2, position = position_dodge2(w = 0.25), alpha = 0.7) +
  geom_hline(yintercept = -log10(sign_pval), color = "grey40", linetype = "dashed") +
  scale_size_continuous(range = c(0.5,3)) + 
  ylab("-log10 (Corrected p-value)") + ylim(0, 5) + 
  labs(title = "B: T1D vs CT") +
  theme_bw()+
  theme(axis.title.x = element_blank(), 
        legend.position = "none",
        axis.text.x = element_text(hjust=0.95,vjust=0.2, size = 12, angle = 90), 
        axis.text.y = element_text(size = 12), axis.title.y = element_text(size = 12),
        plot.title = element_text(size = 15))+
  scale_color_manual(values = c("lavenderblush2", "peachpuff", "lavenderblush2", "wheat2", "skyblue", "wheat2", 
                                "wheat2", "wheat2","wheat2","wheat2", "lavenderblush2", "plum2", "sandybrown",
                                "pink1", "deepskyblue2", "slateblue1", "gray75", "wheat2", "lavenderblush2",
                                "palegreen3", "lavenderblush2"))

data.dm2.ct <- rbind(all_CT_DM2_pos, all_CT_DM2_neg)
data.dm2.ct <- data.dm2.ct[!is.na(data.dm2.ct$LipidIon),]
data.dm2.ct <- data.dm2.ct[data.dm2.ct$confirmed==2,]

data.dm2.ct$Class <- as.factor(as.character(data.dm2.ct$Class))
p.dm2.ct <- ggplot(data.dm2.ct, aes(x = Class, y = -log10(p_adj), color = Class)) + 
  geom_point(size=2, position = position_dodge2(w = 0.25), alpha = 0.7) +
  geom_hline(yintercept = -log10(sign_pval), color = "grey40", linetype = "dashed") +
  scale_size_continuous(range = c(0.5,3)) + 
  ylab("-log10 (Corrected p-value)") + ylim(0, 5) + 
  labs(title = "C: T2D vs CT") +
  theme_bw()+
  theme(axis.title.x = element_blank(), 
        legend.position = "none",
        axis.text.x = element_text(hjust=0.95,vjust=0.2, size = 12, angle = 90), 
        axis.text.y = element_text(size = 12), axis.title.y = element_text(size = 12),
        plot.title = element_text(size = 15))+
  scale_color_manual(values = c("lavenderblush2", "peachpuff", "lavenderblush2", "wheat2", "skyblue", "wheat2", 
                                "wheat2", "wheat2","wheat2","wheat2", "lavenderblush2", "plum2", "sandybrown",
                                "pink1", "deepskyblue2", "slateblue1", "gray75", "wheat2", "lavenderblush2",
                                "palegreen3", "lavenderblush2"))

comb.plot <- ggarrange(p.dm1.dm2, p.dm1.ct, p.dm2.ct, ncol = 3, align = "h")

ggsave("Results/manhattan_plot.png", plot = comb.plot, width = 12, height = 4.5, dpi = 500)

##### Barplot with fold-changes

# formatting
# T1D vs T2D
all_DM1_DM2_pos_women <- all_DM1_DM2_pos_sex[,c("Feature","pval_women","Beta_women","p_adj_women",
                                                "LipidIon","Class","NCarbons","NBonds","Grade",
                                                "lipidsearch_mz_ppm","XCMS_rt","lipidsearch_rt",
                                                "diff_rt","confirmed")]
colnames(all_DM1_DM2_pos_women)[2:4] <- c("pval", "Beta", "p_adj")

all_DM1_DM2_pos_men <- all_DM1_DM2_pos_sex[,c("Feature","pval_men","Beta_men","p_adj_men",
                                              "LipidIon","Class","NCarbons","NBonds","Grade",
                                              "lipidsearch_mz_ppm","XCMS_rt","lipidsearch_rt",
                                              "diff_rt","confirmed")]
colnames(all_DM1_DM2_pos_men)[2:4] <- c("pval", "Beta", "p_adj")


all_DM1_DM2_neg_women <- all_DM1_DM2_neg_sex[,c("Feature","pval_women","Beta_women","p_adj_women",
                                                "LipidIon","Class","NCarbons","NBonds","Grade",
                                                "lipidsearch_mz_ppm","XCMS_rt","lipidsearch_rt",
                                                "diff_rt","confirmed")]
colnames(all_DM1_DM2_neg_women)[2:4] <- c("pval", "Beta", "p_adj")

all_DM1_DM2_neg_men <- all_DM1_DM2_neg_sex[,c("Feature","pval_men","Beta_men","p_adj_men",
                                              "LipidIon","Class","NCarbons","NBonds","Grade",
                                              "lipidsearch_mz_ppm","XCMS_rt","lipidsearch_rt",
                                              "diff_rt","confirmed")]
colnames(all_DM1_DM2_neg_men)[2:4] <- c("pval", "Beta", "p_adj")

# T1D vs CT

all_CT_DM1_pos_women <- all_CT_DM1_pos_sex[,c("Feature","pval_women","Beta_women","p_adj_women",
                                              "LipidIon","Class","NCarbons","NBonds","Grade",
                                              "lipidsearch_mz_ppm","XCMS_rt","lipidsearch_rt",
                                              "diff_rt","confirmed")]
colnames(all_CT_DM1_pos_women)[2:4] <- c("pval", "Beta", "p_adj")

all_CT_DM1_pos_men <- all_CT_DM1_pos_sex[,c("Feature","pval_men","Beta_men","p_adj_men",
                                            "LipidIon","Class","NCarbons","NBonds","Grade",
                                            "lipidsearch_mz_ppm","XCMS_rt","lipidsearch_rt",
                                            "diff_rt","confirmed")]
colnames(all_CT_DM1_pos_men)[2:4] <- c("pval", "Beta", "p_adj")


all_CT_DM1_neg_women <- all_CT_DM1_neg_sex[,c("Feature","pval_women","Beta_women","p_adj_women",
                                              "LipidIon","Class","NCarbons","NBonds","Grade",
                                              "lipidsearch_mz_ppm","XCMS_rt","lipidsearch_rt",
                                              "diff_rt","confirmed")]
colnames(all_CT_DM1_neg_women)[2:4] <- c("pval", "Beta", "p_adj")


all_CT_DM1_neg_men <- all_CT_DM1_neg_sex[,c("Feature","pval_men","Beta_men","p_adj_men",
                                            "LipidIon","Class","NCarbons","NBonds","Grade",
                                            "lipidsearch_mz_ppm","XCMS_rt","lipidsearch_rt",
                                            "diff_rt","confirmed")]
colnames(all_CT_DM1_neg_men)[2:4] <- c("pval", "Beta", "p_adj")

# T2D vs CT

all_CT_DM2_pos_women <- all_CT_DM2_pos_sex[,c("Feature","pval_women","Beta_women","p_adj_women",
                                              "LipidIon","Class","NCarbons","NBonds","Grade",
                                              "lipidsearch_mz_ppm","XCMS_rt","lipidsearch_rt",
                                              "diff_rt","confirmed")]
colnames(all_CT_DM2_pos_women)[2:4] <- c("pval", "Beta", "p_adj")

all_CT_DM2_pos_men <- all_CT_DM2_pos_sex[,c("Feature","pval_men","Beta_men","p_adj_men",
                                            "LipidIon","Class","NCarbons","NBonds","Grade",
                                            "lipidsearch_mz_ppm","XCMS_rt","lipidsearch_rt",
                                            "diff_rt","confirmed")]
colnames(all_CT_DM2_pos_men)[2:4] <- c("pval", "Beta", "p_adj")


all_CT_DM2_neg_women <- all_CT_DM2_neg_sex[,c("Feature","pval_women","Beta_women","p_adj_women",
                                              "LipidIon","Class","NCarbons","NBonds","Grade",
                                              "lipidsearch_mz_ppm","XCMS_rt","lipidsearch_rt",
                                              "diff_rt","confirmed")]
colnames(all_CT_DM2_neg_women)[2:4] <- c("pval", "Beta", "p_adj")

all_CT_DM2_neg_men <- all_CT_DM2_neg_sex[,c("Feature","pval_men","Beta_men","p_adj_men",
                                            "LipidIon","Class","NCarbons","NBonds","Grade",
                                            "lipidsearch_mz_ppm","XCMS_rt","lipidsearch_rt",
                                            "diff_rt","confirmed")]
colnames(all_CT_DM2_neg_men)[2:4] <- c("pval", "Beta", "p_adj")

all_DM1_DM2_pos$uniqueID <- paste(all_DM1_DM2_pos$Feature, all_DM1_DM2_pos$LipidIon, sep = "_")
all_DM1_DM2_pos_women$uniqueID <- paste(all_DM1_DM2_pos_women$Feature, all_DM1_DM2_pos_women$LipidIon, sep = "_")
all_DM1_DM2_pos_men$uniqueID <- paste(all_DM1_DM2_pos_men$Feature, all_DM1_DM2_pos_men$LipidIon, sep = "_")
all_DM1_DM2_neg$uniqueID <- paste(all_DM1_DM2_neg$Feature, all_DM1_DM2_neg$LipidIon, sep = "_")
all_DM1_DM2_neg_women$uniqueID <- paste(all_DM1_DM2_neg_women$Feature, all_DM1_DM2_neg_women$LipidIon, sep = "_")
all_DM1_DM2_neg_men$uniqueID <- paste(all_DM1_DM2_neg_men$Feature, all_DM1_DM2_neg_men$LipidIon, sep = "_")

all_CT_DM1_pos$uniqueID <- paste(all_CT_DM1_pos$Feature, all_CT_DM1_pos$LipidIon, sep = "_")
all_CT_DM1_pos_women$uniqueID <- paste(all_CT_DM1_pos_women$Feature, all_CT_DM1_pos_women$LipidIon, sep = "_")
all_CT_DM1_pos_men$uniqueID <- paste(all_CT_DM1_pos_men$Feature, all_CT_DM1_pos_men$LipidIon, sep = "_")
all_CT_DM1_neg$uniqueID <- paste(all_CT_DM1_neg$Feature, all_CT_DM1_neg$LipidIon, sep = "_")
all_CT_DM1_neg_women$uniqueID <- paste(all_CT_DM1_neg_women$Feature, all_CT_DM1_neg_women$LipidIon, sep = "_")
all_CT_DM1_neg_men$uniqueID <- paste(all_CT_DM1_neg_men$Feature, all_CT_DM1_neg_men$LipidIon, sep = "_")

all_CT_DM2_pos$uniqueID <- paste(all_CT_DM2_pos$Feature, all_CT_DM2_pos$LipidIon, sep = "_")
all_CT_DM2_pos_women$uniqueID <- paste(all_CT_DM2_pos_women$Feature, all_CT_DM2_pos_women$LipidIon, sep = "_")
all_CT_DM2_pos_men$uniqueID <- paste(all_CT_DM2_pos_men$Feature, all_CT_DM2_pos_men$LipidIon, sep = "_")
all_CT_DM2_neg$uniqueID <- paste(all_CT_DM2_neg$Feature, all_CT_DM2_neg$LipidIon, sep = "_")
all_CT_DM2_neg_women$uniqueID <- paste(all_CT_DM2_neg_women$Feature, all_CT_DM2_neg_women$LipidIon, sep = "_")
all_CT_DM2_neg_men$uniqueID <- paste(all_CT_DM2_neg_men$Feature, all_CT_DM2_neg_men$LipidIon, sep = "_")

# binding everything
all_DM1_DM2_pos$Analysis <- "DM1_DM2"
all_DM1_DM2_neg$Analysis <- "DM1_DM2"
all_DM1_DM2_pos_women$Analysis <- "DM1_DM2_women"
all_DM1_DM2_neg_women$Analysis <- "DM1_DM2_women"
all_DM1_DM2_pos_men$Analysis <- "DM1_DM2_men"
all_DM1_DM2_neg_men$Analysis <- "DM1_DM2_men"
all_CT_DM1_pos$Analysis <- "DM1_CT"
all_CT_DM1_neg$Analysis <- "DM1_CT"
all_CT_DM1_pos_women$Analysis <- "DM1_CT_women"
all_CT_DM1_neg_women$Analysis <- "DM1_CT_women"
all_CT_DM1_pos_men$Analysis <- "DM1_CT_men"
all_CT_DM1_neg_men$Analysis <- "DM1_CT_men"
all_CT_DM2_pos$Analysis <- "DM2_CT"
all_CT_DM2_neg$Analysis <- "DM2_CT"
all_CT_DM2_pos_women$Analysis <- "DM2_CT_women"
all_CT_DM2_neg_women$Analysis <- "DM2_CT_women"
all_CT_DM2_pos_men$Analysis <- "DM2_CT_men"
all_CT_DM2_neg_men$Analysis <- "DM2_CT_men"

all_DM1_DM2_pos_id <- all_DM1_DM2_pos_id[all_DM1_DM2_pos_id$confirmed==2,]
all_DM1_DM2_pos_id_women <- all_DM1_DM2_pos_id_women[all_DM1_DM2_pos_id_women$confirmed==2,]
all_DM1_DM2_pos_id_men <- all_DM1_DM2_pos_id_men[all_DM1_DM2_pos_id_men$confirmed==2,]
all_DM1_DM2_neg_id <- all_DM1_DM2_neg_id[all_DM1_DM2_neg_id$confirmed==2,]
all_DM1_DM2_neg_id_women <- all_DM1_DM2_neg_id_women[all_DM1_DM2_neg_id_women$confirmed==2,]
all_DM1_DM2_neg_id_men <- all_DM1_DM2_neg_id_men[all_DM1_DM2_neg_id_men$confirmed==2,]

all_CT_DM1_pos_id <- all_CT_DM1_pos_id[all_CT_DM1_pos_id$confirmed==2,]
all_CT_DM1_pos_id_women <- all_CT_DM1_pos_id_women[all_CT_DM1_pos_id_women$confirmed==2,]
all_CT_DM1_pos_id_men <- all_CT_DM1_pos_id_men[all_CT_DM1_pos_id_men$confirmed==2,]
all_CT_DM1_neg_id <- all_CT_DM1_neg_id[all_CT_DM1_neg_id$confirmed==2,]
all_CT_DM1_neg_id_women <- all_CT_DM1_neg_id_women[all_CT_DM1_neg_id_women$confirmed==2,]
all_CT_DM1_neg_id_men <- all_CT_DM1_neg_id_men[all_CT_DM1_neg_id_men$confirmed==2,]

all_CT_DM2_pos_id <- all_CT_DM2_pos_id[all_CT_DM2_pos_id$confirmed==2,]
all_CT_DM2_pos_id_women <- all_CT_DM2_pos_id_women[all_CT_DM2_pos_id_women$confirmed==2,]
all_CT_DM2_pos_id_men <- all_CT_DM2_pos_id_men[all_CT_DM2_pos_id_men$confirmed==2,]
all_CT_DM2_neg_id <- all_CT_DM2_neg_id[all_CT_DM2_neg_id$confirmed==2,]
all_CT_DM2_neg_id_women <- all_CT_DM2_neg_id_women[all_CT_DM2_neg_id_women$confirmed==2,]
all_CT_DM2_neg_id_men <- all_CT_DM2_neg_id_men[all_CT_DM2_neg_id_men$confirmed==2,]

total.feats <- unique(
  c(paste(as.character(all_DM1_DM2_pos_id$Feature), all_DM1_DM2_pos_id$LipidIon, sep = "_"),
    paste(as.character(all_DM1_DM2_pos_id_women$Feature), all_DM1_DM2_pos_id_women$LipidIon, sep = "_"),
    paste(as.character(all_DM1_DM2_pos_id_men$Feature), all_DM1_DM2_pos_id_men$LipidIon, sep = "_"),
    paste(as.character(all_DM1_DM2_neg_id$Feature), all_DM1_DM2_neg_id$LipidIon, sep = "_"),
    paste(as.character(all_DM1_DM2_neg_id_women$Feature), all_DM1_DM2_neg_id_women$LipidIon, sep = "_"),
    paste(as.character(all_DM1_DM2_neg_id_men$Feature), all_DM1_DM2_neg_id_men$LipidIon, sep = "_"),
    paste(as.character(all_CT_DM1_pos_id$Feature), all_CT_DM1_pos_id$LipidIon, sep = "_"),
    paste(as.character(all_CT_DM1_pos_id_women$Feature), all_CT_DM1_pos_id_women$LipidIon, sep = "_"),
    paste(as.character(all_CT_DM1_pos_id_men$Feature), all_CT_DM1_pos_id_men$LipidIon, sep = "_"),
    paste(as.character(all_CT_DM1_neg_id$Feature), all_CT_DM1_neg_id$LipidIon, sep = "_"),
    paste(as.character(all_CT_DM1_neg_id_women$Feature), all_CT_DM1_neg_id_women$LipidIon, sep = "_"),
    paste(as.character(all_CT_DM1_neg_id_men$Feature), all_CT_DM1_neg_id_men$LipidIon, sep = "_"),
    paste(as.character(all_CT_DM2_pos_id$Feature), all_CT_DM2_pos_id$LipidIon, sep = "_"),
    paste(as.character(all_CT_DM2_pos_id_women$Feature), all_CT_DM2_pos_id_women$LipidIon, sep = "_"),
    paste(as.character(all_CT_DM2_pos_id_men$Feature), all_CT_DM2_pos_id_men$LipidIon, sep = "_"),
    paste(as.character(all_CT_DM2_neg_id$Feature), all_CT_DM2_neg_id$LipidIon, sep = "_"),
    paste(as.character(all_CT_DM2_neg_id_women$Feature), all_CT_DM2_neg_id_women$LipidIon, sep = "_"),
    paste(as.character(all_CT_DM2_neg_id_men$Feature), all_CT_DM2_neg_id_men$LipidIon, sep = "_")))

all_DM1_DM2_pos$Feature <- paste("pos", all_DM1_DM2_pos$Feature, sep = "_")
all_DM1_DM2_pos_women$Feature <- paste("pos", all_DM1_DM2_pos_women$Feature, sep = "_")
all_DM1_DM2_pos_men$Feature <- paste("pos", all_DM1_DM2_pos_men$Feature, sep = "_")
all_DM1_DM2_neg$Feature <- paste("neg", all_DM1_DM2_neg$Feature, sep = "_")
all_DM1_DM2_neg_women$Feature <- paste("neg", all_DM1_DM2_neg_women$Feature, sep = "_")
all_DM1_DM2_neg_men$Feature <- paste("neg", all_DM1_DM2_neg_men$Feature, sep = "_")

all_CT_DM1_pos$Feature <- paste("pos", all_CT_DM1_pos$Feature, sep = "_")
all_CT_DM1_pos_women$Feature <- paste("pos", all_CT_DM1_pos_women$Feature, sep = "_")
all_CT_DM1_pos_men$Feature <- paste("pos", all_CT_DM1_pos_men$Feature, sep = "_")
all_CT_DM1_neg$Feature <- paste("neg", all_CT_DM1_neg$Feature, sep = "_")
all_CT_DM1_neg_women$Feature <- paste("neg", all_CT_DM1_neg_women$Feature, sep = "_")
all_CT_DM1_neg_men$Feature <- paste("neg", all_CT_DM1_neg_men$Feature, sep = "_")

all_CT_DM2_pos$Feature <- paste("pos", all_CT_DM2_pos$Feature, sep = "_")
all_CT_DM2_pos_women$Feature <- paste("pos", all_CT_DM2_pos_women$Feature, sep = "_")
all_CT_DM2_pos_men$Feature <- paste("pos", all_CT_DM2_pos_men$Feature, sep = "_")
all_CT_DM2_neg$Feature <- paste("neg", all_CT_DM2_neg$Feature, sep = "_")
all_CT_DM2_neg_women$Feature <- paste("neg", all_CT_DM2_neg_women$Feature, sep = "_")
all_CT_DM2_neg_men$Feature <- paste("neg", all_CT_DM2_neg_men$Feature, sep = "_")

all_analyses <- rbind(all_DM1_DM2_pos[all_DM1_DM2_pos$uniqueID %in% total.feats,], 
                      all_DM1_DM2_neg[all_DM1_DM2_neg$uniqueID %in% total.feats,],
                      all_DM1_DM2_pos_women[all_DM1_DM2_pos_women$uniqueID %in% total.feats,], 
                      all_DM1_DM2_neg_women[all_DM1_DM2_neg_women$uniqueID %in% total.feats,],
                      all_DM1_DM2_pos_men[all_DM1_DM2_pos_men$uniqueID %in% total.feats,], 
                      all_DM1_DM2_neg_men[all_DM1_DM2_neg_men$uniqueID %in% total.feats,],
                      all_CT_DM1_pos[all_CT_DM1_pos$uniqueID %in% total.feats,], 
                      all_CT_DM1_neg[all_CT_DM1_neg$uniqueID %in% total.feats,],
                      all_CT_DM1_pos_women[all_CT_DM1_pos_women$uniqueID %in% total.feats,], 
                      all_CT_DM1_neg_women[all_CT_DM1_neg_women$uniqueID %in% total.feats,],
                      all_CT_DM1_pos_men[all_CT_DM1_pos_men$uniqueID %in% total.feats,], 
                      all_CT_DM1_neg_men[all_CT_DM1_neg_men$uniqueID %in% total.feats,],
                      all_CT_DM2_pos[all_CT_DM2_pos$uniqueID %in% total.feats,], 
                      all_CT_DM2_neg[all_CT_DM2_neg$uniqueID %in% total.feats,],
                      all_CT_DM2_pos_women[all_CT_DM2_pos_women$uniqueID %in% total.feats,], 
                      all_CT_DM2_neg_women[all_CT_DM2_neg_women$uniqueID %in% total.feats,],
                      all_CT_DM2_pos_men[all_CT_DM2_pos_men$uniqueID %in% total.feats,], 
                      all_CT_DM2_neg_men[all_CT_DM2_neg_men$uniqueID %in% total.feats,])

chexk <- ddply(all_analyses, "uniqueID", function(dx){
  data.frame(unique(dx$uniqueID), nrow(dx))
})

merge.pos <- read.csv("~/paper_dm_newpipeline/Final_results_Lleida_Moll/merge_pos_data.csv", 
                      row.names = 1)
merge.neg <- read.csv("~/paper_dm_newpipeline/Final_results_Lleida_Moll/merge_neg_data.csv", 
                      row.names = 1)

merge.pos[,c(2,6,7,9:11,14:ncol(merge.pos))] <- scale(merge.pos[,c(2,6,7,9:11,14:ncol(merge.pos))], 
                                                      center = TRUE, scale = TRUE)
merge.neg[,c(2,6,7,9:11,14:ncol(merge.neg))] <- scale(merge.neg[,c(2,6,7,9:11,14:ncol(merge.neg))], 
                                                      center = TRUE, scale = TRUE)

library(hash)
h <- hash() 
# set values
h[["DM1"]] <- "T1D"
h[["DM2"]] <- "T2D"
h[["CT"]] <- "Control"
h[["women"]] <- "Women"
h[["men"]] <- "Men"

FC_all_analyses <- ddply(all_analyses, "Analysis", function(dx){
  c1 <- unlist(strsplit(unique(dx$Analysis), split = "_"))[1]
  c2 <- unlist(strsplit(unique(dx$Analysis), split = "_"))[2]
  c1 <- h[[c1]]
  c2 <- h[[c2]]
  if (length(unlist(strsplit(unique(dx$Analysis), split = "_")))>2) c3 <- h[[unlist(strsplit(unique(dx$Analysis), split = "_"))[3]]]
  FC.comp <- ddply(dx, "Feature", function(dy){
    f <- unique(dy$Feature)
    if (!(length(unlist(strsplit(unique(dx$Analysis), split = "_")))>2)){
      if (grepl("pos", f)){
        j <- gsub(pattern = "pos_", replacement = "", f)
        FC <- 10**(mean(merge.pos[merge.pos$DM %in% c1,j],na.rm = TRUE)-mean(merge.pos[merge.pos$DM %in% c2,j],na.rm = TRUE))
        #FC <- (mean(merge.pos[merge.pos$DM %in% c1,j],na.rm = TRUE)/mean(merge.pos[merge.pos$DM %in% c2,j],na.rm = TRUE))
      } else {
        j <- gsub(pattern = "neg_", replacement = "", f)
        FC <- 10**(mean(merge.neg[merge.neg$DM %in% c1,j],na.rm = TRUE)-mean(merge.neg[merge.neg$DM %in% c2,j],na.rm = TRUE))
        #FC <- (mean(merge.neg[merge.neg$DM %in% c1,j],na.rm = TRUE)/mean(merge.neg[merge.neg$DM %in% c2,j],na.rm = TRUE))
      }
    } else {
      merge.pos.sex <- merge.pos[merge.pos$Sex %in% c3,]
      merge.neg.sex <- merge.neg[merge.neg$Sex %in% c3,]
      if (grepl("pos", f)){
        j <- gsub(pattern = "pos_", replacement = "", f)
        FC <- 10**(mean(merge.pos.sex[merge.pos.sex$DM %in% c1,j],na.rm = TRUE)-mean(merge.pos.sex[merge.pos.sex$DM %in% c2,j],na.rm = TRUE))
        #FC <- (mean(merge.pos.sex[merge.pos.sex$DM %in% c1,j],na.rm = TRUE)/mean(merge.pos.sex[merge.pos.sex$DM %in% c2,j],na.rm = TRUE))
      } else {
        j <- gsub(pattern = "neg_", replacement = "", f)
        FC <- 10**(mean(merge.neg.sex[merge.neg.sex$DM %in% c1,j],na.rm = TRUE)-mean(merge.neg.sex[merge.neg.sex$DM %in% c2,j],na.rm = TRUE))
        #FC <- (mean(merge.neg.sex[merge.neg.sex$DM %in% c1,j],na.rm = TRUE)/mean(merge.neg.sex[merge.neg.sex$DM %in% c2,j],na.rm = TRUE))
      }
    }
    
    dy$FC <- FC
    return(dy)
  })
  return(FC.comp)
},.parallel = T)

FC_all_analyses$Class <- as.factor(as.character(FC_all_analyses$Class))
FC_all_analyses$Class2 <- FC_all_analyses$Class

levels(FC_all_analyses$Class) <- c("Cer", "Others", "DG","Others","Others","Others",
                                   "Others","LPC","LPE", 
                                   "PC","PE","PI","SM","Others","TG")
FC_all_analyses$Class <- as.factor(as.character(FC_all_analyses$Class))

freq.class <- as.data.frame(table(as.numeric(FC_all_analyses$Class[!duplicated(FC_all_analyses$uniqueID)])))
freq.class$cumsum <- cumsum(freq.class$Freq)

data_breaks <- data.frame(start = c(0,freq.class$cumsum[-nrow(freq.class)])+0.5,  # Create data with breaks
                          end = freq.class$cumsum+0.5,
                          Class = levels(FC_all_analyses$Class))

FC_all_analyses <- ldply(1:length(unique(FC_all_analyses$Feature)), function(i){
  f <- unique(as.character(FC_all_analyses$Feature))[i]
  dx <- FC_all_analyses[FC_all_analyses$Feature %in% f,]
  dx$newFeat <- paste("F", i, sep = "")
  return(dx)
})

FC_all_analyses$uniqueID <- paste(FC_all_analyses$newFeat, FC_all_analyses$LipidIon, sep = "_")

analysis.names <- c(
  `DM1_CT` = "T1D vs CT",
  `DM1_CT_men` = "T1D vs CT (Men)",
  `DM1_CT_women` = "T1D vs CT (Women)",
  `DM1_DM2` = "T1D vs T2D",
  `DM1_DM2_men` = "T1D vs T2D (Men)",
  `DM1_DM2_women` = "T1D vs T2D (Women)",
  `DM2_CT` = "T2D vs CT",
  `DM2_CT_men` = "T2D vs CT (Men)",
  `DM2_CT_women` = "T2D vs CT (Women)"
)

FC_all_analyses$signif <- ""
FC_all_analyses$signif[FC_all_analyses$p_adj<sign_pval] <- "*"
FC_all_analyses$signif[FC_all_analyses$p_adj<0.01] <- "**"
FC_all_analyses$signif[FC_all_analyses$p_adj<0.001] <- "***"
FC_all_analyses$signif[FC_all_analyses$p_adj<0.0001] <- "****"

FC_all_analyses$nsignif <- unlist(llply(FC_all_analyses$signif, function(i) length(unlist(strsplit(i, split = "")))))
FC_all_analyses$Analysis <- factor(FC_all_analyses$Analysis, levels = c("DM1_DM2","DM1_DM2_men","DM1_DM2_women",
                                                                        "DM1_CT","DM1_CT_men","DM1_CT_women",
                                                                        "DM2_CT","DM2_CT_men","DM2_CT_women"))

FC_all_analyses$uniqueID <- unlist(plyr::llply(FC_all_analyses$uniqueID, function(i){
  if (grepl(pattern = "[+]", x = i)){
    paste0(unlist(strsplit(x = i, split = "[)][+]"))[1],")")
  } else{
    paste0(unlist(strsplit(x = i, split = "[)][-]"))[1],")")
  }
}))

bar.p <- ggplot() +
  scale_x_discrete(seq_len(length(total.feats)))+
  geom_rect(data = data_breaks,
            aes(xmin = start,
                xmax = end,
                ymin = -10,
                ymax = 7.5,
                fill = Class),
            alpha = 0.4)+
  scale_fill_manual(values = c("peachpuff", "skyblue", "plum2", "sandybrown",
                               "wheat2", "pink1", "deepskyblue2", "slateblue1","gray75", "palegreen3"))+ # gray75 lavenderblush2
  guides(fill = guide_legend(nrow = 1))+
  scale_y_continuous(limits = c(-10,7.5))+
  new_scale("fill")+
  geom_bar(data=FC_all_analyses, aes(x=reorder(uniqueID,as.numeric(Class)), y=log2(FC)),
           stat="identity", width=0.6, show.legend = F, fill = "gray40", alpha =0.6) + 
  facet_grid(~Analysis, labeller = labeller(Analysis = analysis.names))+
  scale_fill_manual(values = c("darksalmon", "gray60","cornflowerblue"))+
  theme_classic()+
  coord_flip() + 
  theme(legend.position = "bottom", axis.title.y = element_blank(), 
        axis.ticks.y = element_blank(), axis.text.y = element_text(size = 11),
        strip.text = element_text(size = 14.5), axis.text.x = element_text(size = 12), 
        legend.text = element_text(size = 16), legend.title = element_blank(), 
        axis.title.x = element_text(size = 16), axis.line.y = element_blank())+
  labs(y = "Fold-change")+
  geom_text(data = FC_all_analyses, aes(x=reorder(uniqueID,as.numeric(Class)), y = log2(FC) + (sign(log2(FC))*(1+(nsignif/2))), label = signif), 
            vjust = 0.8,size = 4.5)

ggsave("Results/FC_plot_ggsave_realCT.png", plot = bar.p, width = 23, height = 15, dpi = 300)

## sex differences
setdiff_feats_pos <- c(
  setdiff(unique(all_CT_DM2_pos_id_men$Feature), 
          unique(all_CT_DM2_pos_id_women$Feature)),
  setdiff(unique(all_CT_DM2_pos_id_women$Feature), 
          unique(all_CT_DM2_pos_id_men$Feature))
)

rbind_dm2.ct <- rbind(all_CT_DM2_pos_id_men, all_CT_DM2_pos_id_women)
rbind_dm2.ct <- rbind_dm2.ct[order(rbind_dm2.ct$Class),]
rbind_dm2.ct <- rbind_dm2.ct[rbind_dm2.ct$Feature %in% setdiff_feats_pos,]
setdiff_feats_pos <- rbind_dm2.ct$Feature

setdiff_feats_neg <- c(
  setdiff(unique(all_CT_DM2_neg_id_men$Feature), 
          unique(all_CT_DM2_neg_id_women$Feature)),
  setdiff(unique(all_CT_DM2_neg_id_women$Feature), 
          unique(all_CT_DM2_neg_id_men$Feature))
)

rbind_dm2.ct <- rbind(all_CT_DM2_neg_id_men, all_CT_DM2_neg_id_women)
rbind_dm2.ct <- rbind_dm2.ct[order(rbind_dm2.ct$Class),]
rbind_dm2.ct <- rbind_dm2.ct[rbind_dm2.ct$Feature %in% setdiff_feats_neg,]
setdiff_feats_neg <- rbind_dm2.ct$Feature

plots.sex.dm2.ct.pos <- plyr::llply(setdiff_feats_pos, function(j){
  dx <- all_CT_DM2_pos_sex[all_CT_DM2_pos_sex$Feature %in% j,]
  dx <- dx[dx$confirmed==2,]
  
  lipid.name <- paste(dx$LipidIon, collapse = "\n")
  lipid.name <- paste0(unlist(strsplit(lipid.name, split = "[)+]"))[1], ")")
  dx <- dx[1,]
  m.val <- max(merge.pos[merge.pos$DM %in% c("Control", "T2D"),j], na.rm = TRUE)
  ypos <- m.val+((m.val+2)/10)
  stats <- tibble::tribble(
    ~Sex,~.y., ~group1, ~group2, ~p,  ~y.position,~x, ~xmin, ~xmax,
    "Men",  j,  "Control", "T2D",   round(dx$p_adj_men,4), ypos,1, 0.7, 1.3,
    "Women", j,"Control", "T2D", round(dx$p_adj_women,4), ypos, 2,1.7, 2.3
  )
  stats$signif <- c(bquote(''^ns))
  stats$signif[stats$p<0.05] <- "'*'"
  stats$signif[stats$p<0.01] <- "'**'"
  stats$signif[stats$p<0.001] <- "'***'"
  stats$signif[stats$p<0.0001] <- "'****'"
  
  p <- ggplot(data = merge.pos[merge.pos$DM %in% c("Control", "T2D"),], 
              aes(x = Sex, y = !! sym(j))) + 
    geom_point(size=1, alpha = 0.6, aes(colour = DM),
               position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.7)) +
    geom_hpline(stat = "summary", fun = "mean", width = 0.25, aes(group = DM), 
                position = position_dodge(width = 0.7))+ 
    ggprism::add_pvalue(stats, 
                        xmin = "xmin", 
                        xmax = "xmax",
                        label = "signif", 
                        vjust = 0.25,
                        label.size = 4.5,
                        tip.length = 0.01, 
                        parse = TRUE,
                        bracket.size = 0.3) +
    theme_bw() +
    theme(axis.title.x = element_blank(), axis.text.x = element_text(size = 9), 
          axis.text.y = element_text(size = 9), axis.title.y =  element_text(size = 9),
          legend.position = "bottom", legend.title = element_text(size = 11), legend.text = element_text(size = 10)) +
    scale_colour_manual(labels = c("CT", "T2D"), values = c("darkorchid2","lightcoral")) +
    labs(y = lipid.name) + 
    ylim(min(merge.pos[merge.pos$DM %in% c("Control", "T2D"),j], na.rm = TRUE)-0.1, 
         max(merge.pos[merge.pos$DM %in% c("Control", "T2D"),j], na.rm = TRUE)+1) +
    guides(color = guide_legend(override.aes = list(size = 3)))
  p
},.parallel = T)


j <- setdiff_feats_neg
dx <- all_CT_DM2_neg_sex[all_CT_DM2_neg_sex$Feature %in% j,]
lipid.name <- paste(dx$LipidIon, collapse = "\n")
lipid.name <- paste0(unlist(strsplit(lipid.name, split = "[)+]"))[1], ")")
dx <- dx[1,]
m.val <- max(merge.neg[merge.neg$DM %in% c("Control", "T2D"),j], na.rm = TRUE)
ypos <- m.val+((m.val+2)/10)
stats <- tibble::tribble(
  ~Sex,~.y., ~group1, ~group2, ~p,  ~y.position,~x, ~xmin, ~xmax,
  "Men",  j,  "Control", "T2D",   round(dx$p_adj_men,4), ypos,1, 0.7, 1.3,
  "Women", j,"Control", "T2D", round(dx$p_adj_women,4), ypos, 2,1.7, 2.3
)
stats$signif <- c(bquote(''^ns))
stats$signif[stats$p<0.05] <- "'*'"
stats$signif[stats$p<0.01] <- "'**'"
stats$signif[stats$p<0.001] <- "'***'"
stats$signif[stats$p<0.0001] <- "'****'"

p <- ggplot(data = merge.neg[merge.neg$DM %in% c("Control", "T2D"),], 
            aes(x = Sex, y = !! sym(j))) + 
  geom_point(size=1, alpha = 0.6, aes(colour = DM),
             position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.7)) +
  geom_hpline(stat = "summary", fun = "mean", width = 0.25, aes(group = DM), 
              position = position_dodge(width = 0.7))+ 
  ggprism::add_pvalue(stats, 
                      xmin = "xmin", 
                      xmax = "xmax",
                      label = "signif", 
                      vjust = 0.25,
                      label.size = 4.5,
                      tip.length = 0.01, 
                      parse = TRUE,
                      bracket.size = 0.3) +
  theme_bw() +
  theme(axis.title.x = element_blank(), axis.text.x = element_text(size = 9), 
        axis.text.y = element_text(size = 9), axis.title.y =  element_text(size = 9),
        legend.position = "bottom", legend.title = element_text(size = 11), legend.text = element_text(size = 10)) +
  scale_colour_manual(labels = c("CT", "T2D"), values = c("darkorchid2","lightcoral")) +
  labs(y = lipid.name) + 
  ylim(min(merge.neg[merge.neg$DM %in% c("Control", "T2D"),j], na.rm = TRUE)-0.1, 
       max(merge.neg[merge.neg$DM %in% c("Control", "T2D"),j], na.rm = TRUE)+1) +
  guides(color = guide_legend(override.aes = list(size = 3)))

dm2.ct.sex <- ggarrange(ggarrange(plots.sex.dm2.ct.pos[[1]] + theme(legend.position="none"), 
                                  plots.sex.dm2.ct.pos[[2]] + theme(legend.position="none"),
                                  plots.sex.dm2.ct.pos[[3]] + theme(legend.position="none"), 
                                  plots.sex.dm2.ct.pos[[4]] + theme(legend.position="none"),
                                  plots.sex.dm2.ct.pos[[5]] + theme(legend.position="none"), 
                                  nrow = 1, ncol = 5, widths = c(0.2,0.2,0.2,0.2,0.2)),
                        ggarrange(NULL, plots.sex.dm2.ct.pos[[6]]+ theme(legend.position="none"), 
                                  plots.sex.dm2.ct.pos[[7]]+ theme(legend.position="none"), 
                                  plots.sex.dm2.ct.pos[[8]]+ theme(legend.position="none"), 
                                  p+ theme(legend.position="none"),
                                  NULL, nrow = 1, ncol = 6, widths = c(0.1,0.2,0.2,0.2,0.2,0.1)), 
                        nrow=2, legend.grob = get_legend(plots.sex.dm2.ct.pos[[1]]), legend = "bottom") 

ggsave("Results/sexplot_dm2_ct.png", plot = dm2.ct.sex, width = 9, height = 4, dpi = 400)

all_CT_DM1_pos_id_men$Class <- as.character(all_CT_DM1_pos_id_men$Class)
all_CT_DM1_pos_id_women$Class <- as.character(all_CT_DM1_pos_id_women$Class)

setdiff_feats_pos <- c(
  setdiff(unique(all_CT_DM1_pos_id_men$Feature), 
          unique(all_CT_DM1_pos_id_women$Feature)),
  setdiff(unique(all_CT_DM1_pos_id_women$Feature), 
          unique(all_CT_DM1_pos_id_men$Feature))
)

rbind_dm1.ct <- rbind(all_CT_DM1_pos_id_men, all_CT_DM1_pos_id_women)
rbind_dm1.ct <- rbind_dm1.ct[order(rbind_dm1.ct$Class),]
rbind_dm1.ct <- rbind_dm1.ct[rbind_dm1.ct$Feature %in% setdiff_feats_pos,]
setdiff_feats_pos <- rbind_dm1.ct$Feature

setdiff_feats_neg <- c(
  setdiff(unique(all_CT_DM1_neg_id_men$Feature), 
          unique(all_CT_DM1_neg_id_women$Feature)),
  setdiff(unique(all_CT_DM1_neg_id_women$Feature), 
          unique(all_CT_DM1_neg_id_men$Feature))
)

rbind_dm1.ct <- rbind(all_CT_DM1_neg_id_men, all_CT_DM1_neg_id_women)
rbind_dm1.ct <- rbind_dm1.ct[order(rbind_dm1.ct$Class),]
rbind_dm1.ct <- rbind_dm1.ct[rbind_dm1.ct$Feature %in% setdiff_feats_neg,]
setdiff_feats_neg <- rbind_dm1.ct$Feature

plots.sex.dm1.ct.pos <- plyr::llply(setdiff_feats_pos,.inform = T, function(j){
  dx <- all_CT_DM1_pos_sex[all_CT_DM1_pos_sex$Feature %in% j,]
  dx <- dx[dx$confirmed==2,]
  lipid.name <- paste(dx$LipidIon, collapse = "\n")
  lipid.name <- paste0(unlist(strsplit(lipid.name, split = "[)+]"))[1], ")")
  dx <- dx[1,]
  m.val <- max(merge.pos[merge.pos$DM %in% c("Control", "T1D"),j], na.rm = TRUE)
  ypos <- m.val+((m.val+2)/10)
  stats <- tibble::tribble(
    ~Sex,~.y., ~group1, ~group2, ~p,  ~y.position,~x, ~xmin, ~xmax,
    "Men",  j,  "Control", "T1D",   round(dx$p_adj_men,4), ypos,1, 0.7, 1.3,
    "Women", j,"Control", "T1D", round(dx$p_adj_women,4), ypos, 2,1.7, 2.3
  )
  stats$signif <- c(bquote(''^ns))
  stats$signif[stats$p<0.05] <- "'*'"
  stats$signif[stats$p<0.01] <- "'**'"
  stats$signif[stats$p<0.001] <- "'***'"
  stats$signif[stats$p<0.0001] <- "'****'"
  
  p <- ggplot(data = merge.pos[merge.pos$DM %in% c("Control", "T1D"),], 
              aes(x = Sex, y = !! sym(j))) + 
    geom_point(size=1, alpha = 0.6, aes(colour = DM),
               position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.7)) +
    geom_hpline(stat = "summary", fun = "mean", width = 0.25, aes(group = DM), 
                position = position_dodge(width = 0.7))+ 
    ggprism::add_pvalue(stats, 
                        xmin = "xmin", 
                        xmax = "xmax",
                        label = "signif", 
                        vjust = 0.2,
                        label.size = 4.5,
                        tip.length = 0.01, 
                        parse = TRUE,
                        bracket.size = 0.3) +
    theme_bw() +
    theme(axis.title.x = element_blank(), axis.text.x = element_text(size = 8), 
          axis.text.y = element_text(size = 8), axis.title.y =  element_text(size = 8),
          legend.position = "bottom") +
    scale_colour_manual(labels = c("CT", "T1D"), values = c("darkorchid2","deepskyblue3")) +
    labs(y = lipid.name) + 
    ylim(min(merge.pos[merge.pos$DM %in% c("Control", "T1D"),j], na.rm = TRUE)-0.1, 
         max(merge.pos[merge.pos$DM %in% c("Control", "T1D"),j], na.rm = TRUE)+1.5)+
    guides(color = guide_legend(override.aes = list(size = 3)))
  p
},.parallel = T)


plots.sex.dm1.ct.neg <- plyr::llply(setdiff_feats_neg, function(j){
  dx <- all_CT_DM1_neg_sex[all_CT_DM1_neg_sex$Feature %in% j,]
  dx <- dx[dx$confirmed==2,]
  lipid.name <- paste(dx$LipidIon, collapse = "\n")
  lipid.name <- paste0(unlist(strsplit(lipid.name, split = "[)+]"))[1], ")")
  dx <- dx[1,]
  m.val <- max(merge.neg[merge.neg$DM %in% c("Control", "T1D"),j], na.rm = TRUE)
  ypos <- m.val+((m.val+2)/10)
  stats <- tibble::tribble(
    ~Sex,~.y., ~group1, ~group2, ~p,  ~y.position,~x, ~xmin, ~xmax,
    "Men",  j,  "Control", "T1D",   round(dx$p_adj_men,4), ypos,1, 0.7, 1.3,
    "Women", j,"Control", "T1D", round(dx$p_adj_women,4), ypos, 2,1.7, 2.3
  )
  stats$signif <- c(bquote(''^ns))
  stats$signif[stats$p<0.05] <- "'*'"
  stats$signif[stats$p<0.01] <- "'**'"
  stats$signif[stats$p<0.001] <- "'***'"
  stats$signif[stats$p<0.0001] <- "'****'"
  
  p <- ggplot(data = merge.neg[merge.neg$DM %in% c("Control", "T1D"),], 
              aes(x = Sex, y = !! sym(j))) + 
    geom_point(size=1, alpha = 0.6, aes(fill = DM, colour = DM),
               position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.7)) +
    geom_hpline(stat = "summary", fun = "mean", width = 0.25, aes(group = DM), 
                position = position_dodge(width = 0.7))+ 
    ggprism::add_pvalue(stats, 
                        xmin = "xmin", 
                        xmax = "xmax",
                        label = "signif", 
                        vjust = 0.2,
                        label.size = 4.5,
                        tip.length = 0.01, 
                        parse = TRUE,
                        bracket.size = 0.3) +
    theme_bw() +
    theme(axis.title.x = element_blank(), axis.text.x = element_text(size = 8), 
          axis.text.y = element_text(size = 8), axis.title.y =  element_text(size = 8),
          legend.position = "bottom") +
    scale_colour_manual(labels = c("CT", "T1D"), values = c("darkorchid2","deepskyblue3")) +
    labs(y = lipid.name) + 
    ylim(min(merge.neg[merge.neg$DM %in% c("Control", "T1D"),j], na.rm = TRUE)-0.1, 
         max(merge.neg[merge.neg$DM %in% c("Control", "T1D"),j], na.rm = TRUE)+1.5)+
    guides(color = guide_legend(override.aes = list(size = 3)))
  p
},.parallel = T)

combsex <- ggpubr::ggarrange(plotlist = list(plots.sex.dm1.ct.pos[[1]],
                                             plots.sex.dm1.ct.pos[[2]],
                                             plots.sex.dm1.ct.pos[[3]],
                                             plots.sex.dm1.ct.pos[[4]],
                                             plots.sex.dm1.ct.pos[[5]],
                                             plots.sex.dm1.ct.pos[[13]],
                                             plots.sex.dm1.ct.pos[[14]],
                                             plots.sex.dm1.ct.neg[[6]],
                                             plots.sex.dm1.ct.neg[[7]],
                                             plots.sex.dm1.ct.neg[[8]],
                                             plots.sex.dm1.ct.pos[[6]],
                                             plots.sex.dm1.ct.pos[[7]],
                                             plots.sex.dm1.ct.pos[[8]],
                                             plots.sex.dm1.ct.pos[[9]],
                                             plots.sex.dm1.ct.pos[[10]],
                                             plots.sex.dm1.ct.pos[[12]],
                                             plots.sex.dm1.ct.neg[[2]],
                                             plots.sex.dm1.ct.neg[[3]],
                                             plots.sex.dm1.ct.neg[[4]],
                                             plots.sex.dm1.ct.neg[[5]]), 
                             ncol = 5, nrow = 4, 
                             legend = "bottom", common.legend = T)

ggsave("Results/sexplot_combined_t1d.png", plot = combsex, width = 9, height = 8, dpi = 400)




