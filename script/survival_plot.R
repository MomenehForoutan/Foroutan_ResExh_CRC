
##======================= Author =======================
## ------------ Momeneh (Sepideh) Foroutan -------------
## ------------ First genearted: April 2019 --------------
## ---------- Last updated: 12th of Fe 2020 --------------
##======================================================

# This function depends on the below libraries 
library(survival)
library(ggfortify)
library(RColorBrewer)

#=================== Survival analysis =================
# This function looks at the associations between different variables and survival outcome. These variables can be one of the below options that stratifies samples for survival analysis:

# **expr** : expression of a gene, will be split based on 33%-tile and 66%-tile (e.g. low, medium, high)\n
# **score** : score of a single signatre, will be split based on 33%-tile and 66%-tile (low, medium, high)\n
# **covariate** : A continouse covariate (e.g. age), will be split based on 33%-tile and 66%-tile (low, medium, high)\n
# **score_expr** : stratifies samples based on scores from a signature (high and low) and expression of a gene (high and low)\n
# **covariate_expr** : startifies samples according to covariate (age; high and low) and expression of a gene (high and low)\n
# **score_covariate** : stratifies samples according to scores from a single signature (high and low) and covariate (age; high and low)\n
# **expr_expr** : stratifies samples according to expression of two genes (high gene1/high gene2, high gene1/low gene2, etc)\n
# **score_score** : stratifies samples according to scores obtained from two signatures (high score1/high score2, high score1/low score2, etc)

##---------------------- INPUT ARGUMENTS ------------------------

# 1. **data**: a `SummarizedExperiment` object; for example, here we use the TCGA SKCL data that we downloaded using TCGAbiolink package (see above). This is a specific data object in R that stores expression data as well as several meta-data. Therefore, this function can not take a data frame as input at this stage. The `SummarizedExperiment` object needs to further have an `assay` slot called "logFPKM". There must be an "external_gene_name" column in the rowData that has gene names in the same format as the genes that we provide in the gene arguments of the function. To learn more about SummarizedExperiment object and how to construct it, please see https://bioconductor.org/packages/release/bioc/vignettes/SummarizedExperiment/inst/doc/SummarizedExperiment.html.
# 2. **stratify**: A character value of one of the above listed options for stratification (e.g. "expr", "score_expr").
# 3. **scores**: A data frame with maximum of three columns: one column needs to be called "sample" which has the sample names consistent with sample names in expression data (first argument of teh function), and minimum one or maximum two columns of signature scores, which have "score" as part of their column names. An error will be given if the data has more than two score columns. This argument can be `NULL` if you are not inetersted in the relationship between scores and survival outcome.
# 4. **gene**: A character vector containing the names of maximum of two genes. This argument can be `NULL` if you are not inetersted in the relationship between genes and survival outcome.
# 5. **covariate**: Name of the column for a covariate; This is the age factor by default. This argument can be `NULL` if you are not inetersted in the relationship between the covariate and survival outcome.
# 6. **timeCol**: The name of time column to be used in the survival analysis
# 7. **eventCol**: The name of the event column to be used in survival analysis (e.g. vital_status)
# 8. **nGroup**: The number of groups for each stratification. Can be 2 or 3. For example, a value of 2 (default) generates two groups of samples with high and low expression for a desired gene/score/covariate, while 3 would stratify samples into three groups of high, low, and medium.
# 9. **confInt**: Boolean; if TRUE, the confidence intervals of survival curves are plotted.
## survival analysis
## rnaseq_fpkm is a SummarizedExperiment object downloaded using TCGAbioLink package


## YOU MAY NEED TO DO THIS FOR TCGA DATA bfeore running this function
## As there are many NAs in the days_to_death column, we replace
## the NAs with info from the days_to_last_follow_up

#  timeSurvived <- colData(currentData)$days_to_death
# colData(currentData)$finalTime <-
#   ifelse(is.na(timeSurvived),
#          colData(currentData)$days_to_last_follow_up,
#          timeSurvived)
# 
# currentData <-
#   currentData[, complete.cases(colData(currentData)$finalTime)]


survival_plot <- function(data = exprData,
                          stratify = "score_score",
                          annot = newAnnot,
                          scoreCol =  c("TRM TGFb IL2 Sel Com", "Mes (Byers)"),
                          gene = c("ITGAE", "ZNF683"),
                          covariate = "age_at_initial_pathologic_diagnosis",
                          isCategoricalCov = FALSE,
                          timeCol = "OS.time",
                          eventCol = "OS",
                          nGroup = 2,
                          confInt = F, 
                          ylabel = "Survival",
                          cols = c(brewer.pal(9, "Set1")[c(2, 3, 4, 5, 7, 8)],
                                    brewer.pal(8, "Dark2")[c(8, 1, 4, 6)]), 
                          nColLegend = 1,
                          plotType = "autoplot") {
  
  annot[, timeCol] <- gsub("#N/A", NA, annot[, timeCol])
  
  comSamples <- intersect(rownames(annot), colnames(data))
  data <- data[, comSamples]
  annot <- annot[comSamples, ]
  
  
  data <- data[, complete.cases(annot[, timeCol])]
  annot <- annot[complete.cases(annot[, timeCol]),]
  
  annot[, timeCol] <- as.numeric(annot[, timeCol] )
  
  
  ##---------------------------------------- Scores
  if (!is.null(scoreCol)) {
    
    currentscoreCol <- scoreCol[1]
    median_score <- median(annot[, currentscoreCol])
    lowQ_score <- as.numeric(quantile(annot[, currentscoreCol], prob = 0.33))
    upQ_score <- as.numeric(quantile(annot[, currentscoreCol], prob = 0.66))
    
    # annot$scores_2status <-
    #   ifelse(
    #     annot[, currentscoreCol] >= median_score,
    #     paste("High ", currentscoreCol),
    #     paste("Low ", currentscoreCol)
    #   )
    
    annot$scores_2status[annot[, currentscoreCol] >= median_score] <- paste("High ", currentscoreCol)
    annot$scores_2status[annot[, currentscoreCol] < median_score] <- paste("Low ", currentscoreCol)
    
    annot$scores_3status[annot[, currentscoreCol] >= upQ_score] <-
      paste("High ", currentscoreCol)
    annot$scores_3status[annot[, currentscoreCol] <= lowQ_score] <-
      paste("Low ", currentscoreCol)
    annot$scores_3status[annot[, currentscoreCol] < upQ_score &
                           annot[, currentscoreCol] > lowQ_score] <-
      paste("Medium ", currentscoreCol)
    
    
    if (length(scoreCol) == 2) {
      currentscoreCol <- scoreCol[2]
      median_score <- median(annot[, currentscoreCol])
      lowQ_score <-
        quantile(annot[, currentscoreCol], prob = 0.33)
      upQ_score <-
        quantile(annot[, currentscoreCol], prob = 0.66)
      
      annot$scores2_2status <-
        ifelse(
          annot[, currentscoreCol] >= median_score,
          paste("High ", currentscoreCol),
          paste("Low ", currentscoreCol)
        )
      
      annot$scores2_3status[annot[, currentscoreCol] >= upQ_score] <-
        paste("High ", currentscoreCol)
      annot$scores2_3status[annot[, currentscoreCol] <= lowQ_score] <-
        paste("Low ", currentscoreCol)
      annot$scores2_3status[annot[, currentscoreCol] < upQ_score &
                              annot[, currentscoreCol] > lowQ_score] <-
        paste("Medium ", currentscoreCol)
    }
    if (length(scoreCol) > 2) {
      stop(paste0("You must specify maximum of 2 score columns at a time"))
    }
    ## save this new annotation data as sample annotation for the data
    # colData(data) <- newAnnot
  }
  
  ##-------------------------------------- Covariate
  
  if (!is.null(covariate)) {
   
    
    annot <- annot[complete.cases(annot[, covariate]), ]
    badcols <-
      c(
        "not reported",
        "NA",
        "Indeterminate",
        "[Not Applicable]",
        "[Not Available]",
        "[Discrepancy]",
        "[Unknown]",
        "Not Evaluable"
      )
    annot <- annot[ ! annot[, covariate] %in% badcols, ]
    
    comSamples <- intersect(row.names(annot), colnames(data))
    annot <- annot[comSamples, ]
    data <- data[, comSamples]
    
    
    # newAnnot <- colData(data)
    if(isCategoricalCov){
      annot[, covariate] <- as.factor(annot[, covariate])
    }
    else if(! isCategoricalCov) {
      annot[, covariate] <- as.numeric(annot[, covariate])
      median_cov <- median(annot[, covariate], na.rm = T)
      lowQ_cov <-
        as.numeric(quantile(annot[, covariate], prob = 0.33, na.rm = T))
      upQ_cov <-
        as.numeric(quantile(annot[, covariate], prob = 0.66, na.rm = T))
      
      annot$cov_2status <-
        ifelse(annot[, covariate] >= median_cov,
               "High covariate",
               "Low covariate")
      
      annot$cov_3status[annot[, covariate] >= upQ_cov] <-
        "High covariate"
      annot$cov_3status[annot[, covariate] <= lowQ_cov] <-
        "Low covariate"
      annot$cov_3status[annot[, covariate] < upQ_cov &
                          annot[, covariate] > lowQ_cov] <-
        "Medium covariate"
    
    }
  }
  
  
  currentData <- data
  
  ##------------------------------------- Gene expression
  
  if (!is.null(gene)) {
    if (sum(rownames(data) %in% gene) < 1) {
      stop(paste0(gene, " does not present in the row names of the expression data"))
    }
    if (length(gene) > 3) {
      stop(paste0("Please provide maximum of 3 genes at a time"))
    }
    
    
    currentGene <- gene[1]
    currentGeneIndx <- which(rownames(currentData) == currentGene)
    # newAnnot <- colData(currentData)
    
    ## calculate median and 33%-tile and 66%-tile of gene expression
    
    median_expr <- median(as.numeric(currentData[ currentGeneIndx, ]))
    lowQ_expr <- as.numeric(quantile(currentData[ currentGeneIndx, ], prob = 0.33))
    upQ_expr <- as.numeric(quantile(currentData[ currentGeneIndx, ], prob = 0.66))
    
    # annot$expr1_2status <-
    #   ifelse(
    #     currentData[ currentGeneIndx, ] >= median_expr,
    #     paste0("High ", currentGene),
    #     paste0("Low ", currentGene)
    #   )
    
    annot$expr1_2status[currentData[ currentGeneIndx, ] >= median_expr] <- paste0("High ", currentGene)
    annot$expr1_2status[currentData[ currentGeneIndx, ] < median_expr] <- paste0("Low ", currentGene)
    
    # annot$expr1_3status[as.numeric(currentData[ currentGeneIndx, ]) >= upQ_expr] <-
    #   paste0("High ", currentGene)
    
    annot$expr1_3status[currentData[ currentGeneIndx, ] >= upQ_expr] <-
      paste0("High ", currentGene)

    annot$expr1_3status[currentData[ currentGeneIndx, ] <= lowQ_expr] <-
      paste0("Low ", currentGene)
    annot$expr1_3status[currentData[ currentGeneIndx, ] < upQ_expr &
                         currentData[ currentGeneIndx, ] > lowQ_expr] <-
      paste0("Medium ", currentGene)
    
    ## save this new annotation as sample annotation for the data
    # colData(currentData) <- newAnnot
    
    if (length(gene) > 1) {
      currentGene <- gene[2]
      # newAnnot <- colData(currentData)
      currentGeneIndx <- which(rownames(currentData) == currentGene)
     
       ## calculate median and 33%-tile and 66%-tile of gene expression
      median_expr <-
        median(currentData[currentGeneIndx, ])
      lowQ_expr <-
        as.numeric(quantile(currentData[currentGeneIndx, ], prob = 0.33))
      upQ_expr <-
        as.numeric(quantile(currentData[currentGeneIndx, ], prob = 0.66))
      
      annot$expr2_2status <-
        ifelse(
          currentData[currentGeneIndx, ] >= median_expr,
          paste0("High ", currentGene),
          paste0("Low ", currentGene)
        )
      
      annot$expr2_3status[currentData[currentGeneIndx, ] >= upQ_expr] <-
        paste0("High ", currentGene)
      annot$expr2_3status[currentData[currentGeneIndx, ] <= lowQ_expr] <-
        paste0("Low ", currentGene)
      annot$expr2_3status[currentData[currentGeneIndx, ] < upQ_expr &
                              currentData[currentGeneIndx, ] > lowQ_expr] <-
        paste0("Medium ", currentGene)
      
      if (length(gene) == 3) {
        currentGene <- gene[3]
        currentGeneIndx <- which(rownames(currentData) == currentGene)
        
        ## calculate median and 33%-tile and 66%-tile of gene expression
        median_expr <-
          median(currentData[currentGeneIndx, ])
        lowQ_expr <-
          as.numeric(quantile(currentData[currentGeneIndx, ], prob = 0.33))
        upQ_expr <-
          as.numeric(quantile(currentData[currentGeneIndx, ], prob = 0.66))
        
        annot$expr3_2status <-
          ifelse(
            currentData[currentGeneIndx, ] >= median_expr,
            paste0("High ", currentGene),
            paste0("Low ", currentGene)
          )
        
        annot$expr3_3status[currentData[currentGeneIndx, ] >= upQ_expr] <-
          paste0("High ", currentGene)
        annot$expr3_3status[currentData[currentGeneIndx, ] <= lowQ_expr] <-
          paste0("Low ", currentGene)
        annot$expr3_3status[currentData[currentGeneIndx, ] < upQ_expr &
                              currentData[currentGeneIndx, ] > lowQ_expr] <-
          paste0("Medium ", currentGene)
        
      }
    }
    
    # colData(currentData) <- newAnnot
  }
  
  # currentData <-
  #   currentData[, complete.cases(annot[, timeCol])]
  
  
  ##------------------------------------- Check for stratification type
  if (stratify == "expr") {
    currentStrata <- paste0("expr1_", nGroup, "status")
    mainTitle <- gene[1]
  }
  if(stratify == "score"){
    currentStrata <- paste0("scores_", nGroup, "status")
    mainTitle <- scoreCol[1]
  }
  if(stratify == "covariate"){
    if(isCategoricalCov){
      currentStrata <- covariate
    }
    else if(!isCategoricalCov){
      currentStrata <- paste0("cov_", nGroup, "status")
    }
    mainTitle <- covariate
  }
  if (stratify == "score_expr") {
    currentSt_score <- paste0("scores_", nGroup, "status")
    currentSt_expr  <- paste0("expr1_", nGroup, "status")
    annot$score_expr <-
      paste0(annot[, currentSt_score],
             " / ",
             annot[, currentSt_expr])
    currentStrata <- "score_expr"
    mainTitle <- paste(scoreCol[1], " &\n", gene[1])
  }
  if (stratify == "covariate_expr") {
    if(is.null(covariate) | is.null(gene)){
      stop("Make sure both covriate and gene are provided")
    }
    
    if(isCategoricalCov){
      currentSt_cov <- covariate
    } else if(!isCategoricalCov){
      currentSt_cov   <- paste0("cov_", nGroup, "status")
    }
    
    #currentSt_cov  <- paste0("cov_", nGroup, "status")
    currentSt_expr <- paste0("expr1_", nGroup, "status")
    ## Remove samples with NA annotation as covariate:
    currentData <- currentData[ , complete.cases(annot[, currentSt_cov])]
    annot$cov_expr <-
      paste0(annot[, currentSt_cov],
             " / ",
             annot[, currentSt_expr])
    currentStrata <- "cov_expr"
    mainTitle <- paste(covariate, " &\n", gene[1])
  }
  if (stratify == "score_covariate") {
    
    if(isCategoricalCov){
      currentSt_cov <- covariate
    } else if(!isCategoricalCov){
      currentSt_cov   <- paste0("cov_", nGroup, "status")
    }
    
    currentSt_score <- paste0("scores_", nGroup, "status")

    annot$score_cov <-
      paste0(annot[, currentSt_score],
             " / ",
             annot[, currentSt_cov])
    currentStrata <- "score_cov"
    mainTitle <- paste(scoreCol[1], " &\n", covariate)
  }
  if (stratify == "expr_expr") {
    currentSt_expr1 <- paste0("expr1_", nGroup, "status")
    currentSt_expr2 <- paste0("expr2_", nGroup, "status")
    annot$expr_expr <-
      paste0(annot[, currentSt_expr1],
             " / ",
             annot[, currentSt_expr2])
    currentStrata <- "expr_expr"
    mainTitle <- paste(gene[1], " &\n", gene[2])
  }
  if (stratify == "score_score") {
    currentSt_score1 <- paste0("scores_", nGroup, "status")
    currentSt_score2 <- paste0("scores2_", nGroup, "status")
    annot$score_score <-
      paste0(
        annot[, currentSt_score1],
        " / ",
        annot[, currentSt_score2])
    currentStrata <- "score_score"
    mainTitle <- paste(scoreCol[1], " &\n", scoreCol[2])
  }
  
  if (stratify == "expr_expr_expr") {
    currentSt_expr1 <- paste0("expr1_", nGroup, "status")
    currentSt_expr2 <- paste0("expr2_", nGroup, "status")
    currentSt_expr3 <- paste0("expr3_", nGroup, "status")
    annot$expr_expr_expr <-
      paste0(annot[, currentSt_expr1],
             " / ",
             annot[, currentSt_expr2],
             " / ",
             annot[, currentSt_expr3])
    currentStrata <- "expr_expr_expr"
    mainTitle <- paste(gene[1], " &\n", gene[2], " & ", gene[3])
  }
  
  
  ##----------------------------- Fit survival model
  # annot2 <- annot
  # annot2[, timeCol] <- gsub("#N/A", NA, annot2[, timeCol])
  # annot2 <- annot2[complete.cases(annot2[, timeCol]), ]
  
  
  tt <- data.frame(table(annot[, currentStrata]))
  tt$Var1 <- as.character(tt$Var1)
  tt$Freq <- as.character(tt$Freq)
  
  
  for (i in 1:nrow(tt)) {
    annot$currentStrata_n[annot[, currentStrata] == tt$Var1[i]] <-
      paste0(tt$Var1[i], " (", tt$Freq[i], ")")
  }
  
  annot <- annot[, c("currentStrata_n", timeCol, eventCol)]
  
  
  annot[, timeCol] <- as.numeric(annot[, timeCol])
  annot[, eventCol] <- as.numeric(annot[, eventCol])
  
  fitValues <- survfit(Surv(time = annot[, timeCol],
                            event = annot[, eventCol]) ~
                         annot$currentStrata_n)

  ss <- survdiff(Surv(time =  annot[, timeCol],
                      event = annot[, eventCol]
                      ) ~
                   annot$currentStrata_n)
# 
#   ss <- survdiff(Surv(annot[, timeCol],
#                       as.numeric(as.factor(
#                         annot[, eventCol]
#                       )) - 1) ~
#                    annot$currentStrata_n)
  
  ##------------------------------- Calculate p-value
  ## This does not adjust for any covariates, unless the covariate option is included
  
  pval <-  ifelse (is.na(ss), next, (round(1 - pchisq(
    ss$chisq, length(ss$n) - 1
  ), 6)))[[1]]
  
  if(pval < 0.01){
    pval_add <- "p-value (log-rank) < 0.01"
  }
  else{
    pval_add <- paste0("p-value (log-rank) = ", round(pval, 2))
  }

  
  ##------------------------------ Plot survival curve
  
  if(plotType == "autoplot"){
  p <-  autoplot(fitValues, surv.size = 1.5, conf.int = confInt) +
    scale_color_manual(values = cols) +
    scale_fill_manual(values = cols) +
    ggtitle(paste0(mainTitle,
                   # " (Chisq = ", round(ss$chisq, 3),
                   "\n", pval_add)) +
    ylab(ylabel) +
    xlab("Time") +
    theme_bw() +
    theme(legend.position = "bottom") +
    guides(color = guide_legend(ncol = nColLegend), 
           fill = guide_legend(ncol = nColLegend))
  
  }
  else if (plotType == "ggsurvplot"){
    p <- ggsurvplot (
      fitValues,
      data = annot,
      fun = "pct",
      pval = TRUE,
      # pval.method = TRUE,  ## Log Rank
      # test.for.trend = T,  ## when we have more than two groups
      conf.int = confInt,
      surv.median.line = "hv",
      # linetype = "strata",
      palette = cols,
      xlab = "Time",
      legend.title = mainTitle,
      # legend.labs = c("High score", "Low score"),
      # legend = c(.2, .2),
      # break.time.by = 4,
      # risk.table = TRUE,
      # tables.height = 0.2,
      # tables.theme = theme_cleantable(),
      # risk.table.y.text.col = TRUE,
      # risk.table.y.text = TRUE
    )
  }
  p_pval <- list(plot = p, pval = pval)
  return(p_pval)
  
}






survival_pairPlots <- function(data1 = selScoreExprPair_tcga,
                               data2 = selScoreExprPair_mar,
                               strataName = unique(comOS_scoreExprPairs)[1],
                               strataNameColumn = Comb_Name,
                               strataColumn = Comb_2status,
                               mainTitle1 = "TCGA\n",
                               mainTitle2 = "Marisa\n",
                               timeCol1 = "OS.time",
                               eventCol1 = "OS",
                               timeCol2 = "OS.time",
                               eventCol2 = "OS",
                               nColLegend = 1,
                               confInt = F,
                               ylabel1 = "OS",
                               ylabel2 = "OS",
                               CoxPval = " < 0.05", 
                               cols = cols, 
                               currentTheme = currentTheme
                               
){
  
  
  i <-  strataName
  annot <- data1[data1[, strataNameColumn] == i,]
  annot[, timeCol1] <- as.numeric(annot[, timeCol1])
  annot[, eventCol1] <- as.numeric(annot[, eventCol1])
  
  
  tt <- data.frame(table(annot[, strataColumn]))
  tt$Var1 <- as.character(tt$Var1)
  tt$Freq <- as.character(tt$Freq)
  
  
  for (i in 1:nrow(tt)) {
    annot$currentStrata_n[annot[, strataColumn] == tt$Var1[i]] <-
      paste0(tt$Var1[i], " (", tt$Freq[i], ")")
  }
  
  annot <- annot[, c("currentStrata_n", timeCol1, eventCol1)]
  
  # currentPval <- cox_PFI_2Comb_AgeStage_Extremes$`Pr(>|z|)`[
  #   cox_PFI_2Comb_AgeStage_Extremes$VariableName == i
  # ]
  # if(currentPval < 0.01) {
  #   currentPval <- " < 0.01"
  # }
  # else{
  #   currentPval <-  paste0(" = ", round(currentPval, 2))
  # }
  
  fitValues <- survfit(Surv(time = annot[, timeCol1],
                            event = annot[, eventCol1]) ~
                         annot$currentStrata_n)

  
  ss <- survdiff(Surv(time =  annot[, timeCol1],
                      event = annot[, eventCol1]) ~
                   annot$currentStrata_n)
  
  currentPval <-  CoxPval
  mainTitle1 <- mainTitle1
  
  ##------------------------------- Calculate p-value
  ## This does not adjust for any covariates, unless the covariate option is included
  
  pval <-  ifelse (is.na(ss), next, (round(1 - pchisq(
    ss$chisq, length(ss$n) - 1
  ), 6)))[[1]]
  
  if (pval < 0.01) {
    pval <- " < 0.01"
  } else if(pval >= 0.01) {
    pval <- paste0(" = ", round(pval, 2))
  }
  

  p1 <- autoplot(fitValues, surv.size = 1.5, conf.int = confInt) +
    scale_color_manual(values = cols) +
    scale_fill_manual(values = cols) +
    ggtitle(
      paste0(
        mainTitle1,
        "p (log-rank)",
        pval,
        "; p (Cox)",
        currentPval
        # "\nNumber of samples in each group = ",
        # as.numeric(table(annot$Comb_2status)[1])
      )) +
    ylab(ylabel1) +
    xlab("Time") +
    # theme_bw() +
    currentTheme +
    theme(
      legend.position = "bottom",
      plot.title = element_text(
        face = "plain",
        size = rel(textSize),
        hjust = 0.5
      )
    ) +
    scale_y_continuous(limits=c(0, 1)) +
    guides(
      color = guide_legend(ncol = nColLegend, title = NULL),
      fill = guide_legend(ncol = nColLegend, title = NULL)
    )
  
  
  
  ##-------------------- Marisa
  
  i <-  strataName
  annot <-
    data2[data2[, strataNameColumn] == i,]
  
  annot[, timeCol2] <- as.numeric(annot[, timeCol2])
  annot[, eventCol2] <-
    as.numeric(annot[, eventCol2])
  
  
  tt <- data.frame(table(annot[, strataColumn]))
  tt$Var1 <- as.character(tt$Var1)
  tt$Freq <- as.character(tt$Freq)
  
  
  for (i in 1:nrow(tt)) {
    annot$currentStrata_n[annot[, strataColumn] == tt$Var1[i]] <-
      paste0(tt$Var1[i], " (", tt$Freq[i], ")")
  }
  
  annot <- annot[, c("currentStrata_n", timeCol2, eventCol2)]
  
  # currentPval <- cox_PFI_2Comb_AgeStage_Extremes$`Pr(>|z|)`[
  #   cox_PFI_2Comb_AgeStage_Extremes$VariableName == i
  # ]
  # if(currentPval < 0.01) {
  #   currentPval <- " < 0.01"
  # }
  # else{
  #   currentPval <-  paste0(" = ", round(currentPval, 2))
  # }
  
  currentPval <-  CoxPval

  fitValues <-
    survfit(Surv(time = annot[, timeCol2],
                 event = annot[, eventCol2]) ~
              annot$currentStrata_n)
  
  mainTitle2 <- mainTitle2
  
  ss <- survdiff(Surv(time =  annot[, timeCol2],
                      event = annot[, eventCol2]) ~
                   annot$currentStrata_n)
  
  ##------------------------------- Calculate p-value
  ## This does not adjust for any covariates, unless the covariate option is included
  
  pval <-
    ifelse (is.na(ss), next, (round(
      1 - pchisq(ss$chisq, length(ss$n) - 1), 6
    )))[[1]]
  
  if (pval < 0.01) {
    pval <- " < 0.01"
  }
  else {
    pval <- paste0(" = ", round(pval, 2))
  }
  
  
  
  p2 <-
    autoplot(fitValues, surv.size = 1.5, conf.int = confInt) +
    scale_color_manual(values = cols) +
    scale_fill_manual(values = cols) +
    ggtitle(
      paste0(
        mainTitle2,
        "p (log-rank)",
        pval,
        "; p (Cox)",
        currentPval
        # "\nNumber of samples in each group = ",
        # as.numeric(table(annot$Comb_2status)[1])
      )) +
    ylab(ylabel2) +
    xlab("Time") +
    # theme_bw() +
    currentTheme +
    theme(
      legend.position = "bottom",
      plot.title = element_text(
        face = "plain",
        size = rel(textSize),
        hjust = 0.5
      )
    ) +
    scale_y_continuous(limits=c(0, 1)) +
    guides(
      color = guide_legend(ncol = nColLegend, title = NULL),
      fill = guide_legend(ncol = nColLegend, title = NULL)
    )
  
  grid.arrange(p1, p2,
               ncol = 2)
  
}





