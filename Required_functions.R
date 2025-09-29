
# Function to prepare LipidSearch identification data
lipidsearch_prep <- function(df){
  df$Grade <- gsub(x = df$Grade, pattern = "NA", replacement = "")
  colnames(df)[colnames(df) %in% "XCMS_ID"] <- "Feature"
  df <- data.frame(uniqueID = paste(df$Feature, df$LipidIon, sep = "_"), df)
  df$diff_rt <- df$lipidsearch_rt-df$XCMS_rt
  
  df <- ldply(unique(df$uniqueID),.inform = T, function(i){
    dx <- df[df$uniqueID %in% i,]
    dx$confidence <- 0
    if (nrow(dx)>1){
      dx$confidence[which.min(abs(dx$diff_rt))] <- dx$confidence[which.min(abs(dx$diff_rt))]+1
      dx$confidence[which.min(abs(dx$lipidsearch_mz_ppm))] <- dx$confidence[which.min(abs(dx$lipidsearch_mz_ppm))]+1
      if (sum(dx$Rej.)>0)  dx$confidence[dx$Rej.==1] <- dx$confidence[dx$Rej.==1]-1
    }
    dx
  })
  df.conf <- ldply(unique(df$uniqueID),.inform = T, function(i){
    dx <- df[df$uniqueID %in% i,]
    dx[which.max(dx$confidence),]
  })
  df.conf$confirmed <- 1
  df.conf <- ldply(1:nrow(df.conf), function(i){
    dx <- df.conf[i,]
    grade <- ifelse(grepl(pattern = "A", dx$Grade) | grepl(pattern = "B", dx$Grade), 1, 0)
    ppm <- ifelse(abs(dx$lipidsearch_mz_ppm)<5, 1, 0)
    rt <- ifelse(abs(dx$diff_rt)<5, 1, 0)
    if (sum(grade, ppm, rt)==3) {
      dx$confirmed<-2
    } else if ((sum(ppm,rt)==2) & (grade==0) & (grepl(pattern = "C", dx$Grade))){
      dx$confirmed <- 1
    } else {
      dx$confirmed <- 0
    }
    return(dx)
  })
  
  lipid.props <- ldply(df.conf$uniqueID, function(i){
    dx <- df.conf[df.conf$uniqueID %in% i,]
    ion <- as.character(dx$LipidGroup)
    props <- unlist(strsplit(unlist(strsplit(ion, split = "[(]")), split = "[)]"))[2]
    n1 <- unlist(strsplit(props, split = "[:]"))[1]
    n1 <- gsub(pattern = "d", replacement = "", x = n1)
    n1 <- gsub(pattern = "m", replacement = "", x = n1)
    n1 <- gsub(pattern = "t", replacement = "", x = n1)
    nCarbons <- as.numeric(n1)
    n2 <- unlist(strsplit(props, split = "[:]"))[2]
    n2 <- gsub(pattern = "[+]O", replacement = "", x = n2)
    n2 <- gsub(pattern = "p", replacement = "", x = n2)
    nBonds <- as.numeric(n2)
    data.frame(uniqueID = i, LipidGroup = ion, NCarbons = nCarbons, NBonds = nBonds)
  },.parallel = TRUE)
  
  lipidsearch <- merge(lipid.props[,-2], df.conf, by = "uniqueID")
  
  return(lipidsearch)
}

## Function to perform the linear models analysis
## group_subjects: always put the basseline level at the beginning, 
## possibilities are: c(c("T1D", "T2D"), c("CT", "T1D"), c("CT", "T2D"))
linear_analysis <- function(data, feat_cols, confounders_str, 
                            group_subjects, doPar = TRUE){
  data_spec <- data[data$DM %in% group_subjects,]
  data_spec$DM <- as.factor(as.character(data_spec$DM))
  pvals_df <- plyr::ldply(feat_cols, function(j){
    f <- paste(j, "~", confounders_str)
    f <- as.formula(f)
    model <- lm(f, data = data_spec)
    s <- summary(model)
    var_name <- paste0("DM", group_subjects[2])
    pval <- s$coefficients[var_name,4]
    Beta <- s$coefficients[var_name,1]
    data.frame(Feature = j, pval = pval, Beta = Beta)
  },.parallel = doPar)
  pvals_df$p_adj <- p.adjust(pvals_df$pval, "fdr")
  return(pvals_df)
}

## Function to perform the linear models with interactions analysis
linear_analysis_interactions <- function(data, feat_cols, confounders_str, 
                                         group_subjects, doPar = TRUE){
  data_spec <- data[data$DM %in% group_subjects,]
  data_spec$DM <- as.factor(as.character(data_spec$DM))
  pvals_df <- plyr::ldply(feat_cols, function(j){
    f <- paste(j, "~", confounders_str, "+Sex:DM")
    f <- as.formula(f)
    model <- lm(f, data = data_spec)
    s <- summary(model)
    var_name <- paste0("DM", group_subjects[2])
    
    pval_men <- s$coefficients[var_name,4]
    Beta_men <- s$coefficients[var_name,1]
    
    Beta_women <- sum(s$coefficients[var_name,1],s$coefficients[paste0("SexWomen:",var_name),1])
    a <- rep(0,nrow(na.omit(s$coefficients)))
    a[rownames(s$coefficients) %in% var_name] <- 1
    a[rownames(s$coefficients) %in% paste0("SexWomen:",var_name)] <- 1
    se <- c(sqrt(t(a) %*% vcov(model) %*% a))
    tvalue_women <- Beta_women/se
    pval_women <- pt(abs(tvalue_women), df = s$df[2], lower.tail = F)*2
    data.frame(Feature = j, pval_men = pval_men, Beta_men = Beta_men,
               pval_women = pval_women, Beta_women = Beta_women)
  },.parallel = doPar)
  
  pvals_df$p_adj_men <- p.adjust(pvals_df$pval_men, "fdr")
  pvals_df$p_adj_women <- p.adjust(pvals_df$pval_women, "fdr")
  return(pvals_df)
}

## geom_hpline function
geom_hpline <- function(mapping = NULL, data = NULL,
                        stat = "identity", position = "identity",
                        ...,
                        na.rm = FALSE,
                        show.legend = NA,
                        inherit.aes = TRUE) {
  layer(
    data = data,
    mapping = mapping,
    stat = stat,
    geom = GeomHpline,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(
      na.rm = na.rm,
      ...
    )
  )
}

#' @rdname geom_hpline
#' @format NULL
#' @usage NULL
#' @export
GeomHpline <- ggproto("GeomHpline", GeomSegment,
                      required_aes = c("x", "y"),
                      non_missing_aes = c("size", "colour", "linetype", "width"),
                      default_aes = aes(
                        width = 0.5, colour = "black", size = 0.5, linetype = 1,
                        alpha = NA
                      ),
                      
                      draw_panel = function(self, data, panel_params, coord, arrow = NULL, arrow.fill = NULL,
                                            lineend = "butt", linejoin = "round", na.rm = FALSE) {
                        data <- mutate(data, x = x - width/2, xend = x + width, yend = y)
                        ggproto_parent(GeomSegment, self)$draw_panel(
                          data, panel_params, coord, arrow = arrow, arrow.fill = arrow.fill,
                          lineend = lineend, linejoin = linejoin, na.rm = na.rm
                        )
                      }
)

sign_features <- function(pvals_data, variable_pval, analysis_name, ionization, 
                          sign_pval = 0.05, do.Par = TRUE){
  sign_df <- pvals_data[(pvals_data[,variable_pval]<sign_pval) & (!(is.na(pvals_data$LipidIon))),]
  sign_df <- sign_df[sign_df$confirmed==2,]
  sign_df$Analysis <- analysis_name
  sign_df$Ionization <- ionization
  sign_df <- sign_df[,c("Feature", "Class", "LipidIon", "Analysis", 
                        variable_pval, "Ionization")]
  colnames(sign_df)[5] <- "q-value"
  sign_df <- ldply(unique(as.character(sign_df$Feature)), function(i){
    dx <- sign_df[sign_df$Feature %in% i,]
    if (nrow(dx)>1){
      dy <- dx[1,]
      dy$LipidIon <- paste(unique(dx$LipidIon), collapse = "|")
      dy$Class <- paste(unique(dx$Class), collapse = "|")
      return(dy)
    } else{
      return(dx)
    }
  },.parallel = do.Par)
}
