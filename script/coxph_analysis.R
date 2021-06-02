library(survival)
library(foreach)



##----- perform Cox analysis and return either the stats or the models themselves
# 
coxph_test <-
  function(data = scoreAnnot,
           timeCol = timeCol,
           eventCol = eventCol,
           # sigCol = "Signature",
           variableCol = scoreCol,
           nStrata = 4, ## how many groups in the main variable
           returnCovarStat = F, ## whether or not also return the stats for the covariates
           covarCol = c("Age", "Stage"),
           returnCoxModels = F) {

    if (!is.null(covarCol)) {
      covarInfo <- paste(covarCol, collapse = " + ")

      currentFormula <-
        as.formula(paste(
          'Surv(',
          timeCol,
          ', ',
          eventCol,
          ')~',
          variableCol,
          ' + ',
          covarInfo
        ))

      coxFit <-
        coxph(currentFormula,
              data = data)

      if(returnCoxModels){
        return(coxFit)
      }

      else if(!returnCoxModels){

        resNames <-
          c(
            "variable",
            colnames(summary(coxFit)$coef),
            colnames(summary(coxFit)$conf.int)[2:4],
            names(summary(coxFit)$concordance),
            names(summary(coxFit)$sctest)
          )

        resData <-
          data.frame(matrix(ncol = length(resNames), nrow = 0))
        colnames(resData) <- resNames


        ## Get stats for all variables including covariates
        if (returnCovarStat) {
          for (i in rownames(summary(coxFit)$coef)) {
            currentCoxRes <-  c(
              summary(coxFit)$coef[i, ],
              ## first column is exp(coef) which we extracted in the first line
              summary(coxFit)$conf.int[i, 2:4],
              ## C is concordance
              summary(coxFit)$concordance,
              ## p-value is the overall statistical test
              summary(coxFit)$sctest
            )

            currentCoxRes <-
              data.frame(variable = i,
                         t(data.frame(currentCoxRes)),
                         check.names = F)


            resData <- rbind(resData, currentCoxRes)
          }
          resData$formula <-
            as.character(currentFormula)[3]
        }
        ## Get three first line of the results as we have 4 strata
        if (!returnCovarStat) {
          for (i in rownames(summary(coxFit)$coef)[1:(nStrata - 1)]) {
            currentCoxRes <-  c(
              summary(coxFit)$coef[i, ],
              ## first column is exp(coef) which we extracted in the first line
              summary(coxFit)$conf.int[i, 2:4],
              ## C is concordance
              summary(coxFit)$concordance,
              ## p-value is the overall statistical test
              summary(coxFit)$sctest
            )

            currentCoxRes <-
              data.frame(variable = i,
                         t(data.frame(currentCoxRes)),
                         check.names = F)


            resData <- rbind(resData, currentCoxRes)
          }

          resData$formula <-
            as.character(currentFormula)[3]
        }
        resData$variable <-
          gsub(variableCol, "", resData$variable)
        return(resData)


      }
    }


    if (is.null(covarCol)) {

      currentFormula <-
        as.formula(paste(
          'Surv(',
          timeCol,
          ', ',
          eventCol,
          ')~',
          variableCol
        ))

      coxFit <-
        coxph(currentFormula,
              data = data)

      if(returnCoxModels){
        return(coxFit)
      }
      else if(!returnCoxModels){

        resNames <-
          c(
            "variable",
            colnames(summary(coxFit)$coef),
            colnames(summary(coxFit)$conf.int)[2:4],
            names(summary(coxFit)$concordance),
            names(summary(coxFit)$sctest)
          )

        resData <-
          data.frame(matrix(ncol = length(resNames), nrow = 0))
        colnames(resData) <- resNames


        for (i in rownames(summary(coxFit)$coef)) {
          currentCoxRes <-  c(
            summary(coxFit)$coef[i, ],
            ## first column is exp(coef) which we extracted in the first line
            summary(coxFit)$conf.int[i, 2:4],
            ## C is concordance
            summary(coxFit)$concordance,
            ## p-value is the overall statistical test
            summary(coxFit)$sctest
          )

          currentCoxRes <-
            data.frame(variable = i,
                       t(data.frame(currentCoxRes)),
                       check.names = F)


          resData <- rbind(resData, currentCoxRes)
        }

        resData$formula <-
          as.character(currentFormula)[3]

        resData$variable <-
          gsub(variableCol, "", resData$variable)
        return(resData)

      }

    }
  }



##----- If we want to extract contrasts in the model:

coxph_test_contrasts <-
  function(data = scoreAnnot,
           timeCol = timeCol,
           eventCol = eventCol,
           # sigCol = "Signature",
           variableCol = scoreCol,
           nStrata = 4, ## how many groups in the main variable
           returnCovarStat = F, ## whether or not also return the stats for the covariates
           covarCol = c("Age", "Stage"),
           returnCoxModels = F) {
    
    if (!is.null(covarCol)) {
      covarInfo <- paste(covarCol, collapse = " + ")
      
      currentFormula <-
        as.formula(paste(
          'Surv(',
          timeCol,
          ', ',
          eventCol,
          ')~',
          variableCol,
          ' + ',
          covarInfo
        ))
      
      ## added this to consider the contrasts
      data[, variableCol] <- as.factor(data[, variableCol])
      contrasts(data[, variableCol]) <- contr.helmert(4)
      
      # mm <- matrix(
      #   c(3,-1,-1,-1 , 
      #     -1,3,-1,-1, 
      #     -1,-1,3,-1, 
      #     -1,-1,-1,3), 
      #   ncol = 4)
      # 
      # rownames(mm) <- levels(data[, variableCol])
      # 
      # contrasts(data[, variableCol]) <- mm
      
      # contrasts(data[, variableCol])
      #                                             [,1] [,2] [,3]
      # High NK_Res_Bulk_Str & High NK_Exh_Bulk_Str   -1   -1   -1
      # High NK_Res_Bulk_Str & Low NK_Exh_Bulk_Str     1   -1   -1
      # Low NK_Res_Bulk_Str & High NK_Exh_Bulk_Str     0    2   -1
      # Low NK_Res_Bulk_Str & Low NK_Exh_Bulk_Str      0    0    3
      
      coxFit <-
        coxph(currentFormula,
              data = data)
      
      if(returnCoxModels){
        return(coxFit)
      }
      
      else if(!returnCoxModels){
        
        resNames <-
          c(
            "variable",
            colnames(summary(coxFit)$coef),
            colnames(summary(coxFit)$conf.int)[2:4],
            names(summary(coxFit)$concordance),
            names(summary(coxFit)$sctest)
          )
        
        resData <-
          data.frame(matrix(ncol = length(resNames), nrow = 0))
        colnames(resData) <- resNames
        
        
        ## Get stats for all variables including covariates
        if (returnCovarStat) {
          for (i in rownames(summary(coxFit)$coef)) {
            currentCoxRes <-  c(
              summary(coxFit)$coef[i, ],
              ## first column is exp(coef) which we extracted in the first line
              summary(coxFit)$conf.int[i, 2:4],
              ## C is concordance
              summary(coxFit)$concordance,
              ## p-value is the overall statistical test
              summary(coxFit)$sctest
            )
            
            currentCoxRes <-
              data.frame(variable = i,
                         t(data.frame(currentCoxRes)),
                         check.names = F)
            
            
            resData <- rbind(resData, currentCoxRes)
          }
          resData$formula <-
            as.character(currentFormula)[3]
        }
        ## Get three first line of the results as we have 4 strata
        if (!returnCovarStat) {
          for (i in rownames(summary(coxFit)$coef)[1:(nStrata - 1)]) {
            currentCoxRes <-  c(
              summary(coxFit)$coef[i, ],
              ## first column is exp(coef) which we extracted in the first line
              summary(coxFit)$conf.int[i, 2:4],
              ## C is concordance
              summary(coxFit)$concordance,
              ## p-value is the overall statistical test
              summary(coxFit)$sctest
            )
            
            currentCoxRes <-
              data.frame(variable = i,
                         t(data.frame(currentCoxRes)),
                         check.names = F)
            
            
            resData <- rbind(resData, currentCoxRes)
          }
          
          resData$formula <-
            as.character(currentFormula)[3]
          
          ## add proper contrast names:
          
          if(nStrata == 4){
            resData$contrast.helmert[1] <- 
              paste0(
                row.names(contrasts(data[, variableCol]))[2], " - ", 
                row.names(contrasts(data[, variableCol]))[1])
            
            resData$contrast.helmert[2] <- 
              paste0(
                row.names(contrasts(data[, variableCol]))[3], " - ", 
                row.names(contrasts(data[, variableCol]))[1], " AND ",
                  row.names(contrasts(data[, variableCol]))[2]   )
            
            resData$contrast.helmert[3] <- 
              paste0(
                row.names(contrasts(data[, variableCol]))[3], " - ", 
                row.names(contrasts(data[, variableCol]))[1], " AND ",
                row.names(contrasts(data[, variableCol]))[2], " AND ",
                row.names(contrasts(data[, variableCol]))[3]  )
          }

        }
        # resData$variable <-
        #   gsub(variablvCol, "", resData$variable)
        return(resData)
        
        
      }
    }
    
    
    if (is.null(covarCol)) {
      
      currentFormula <-
        as.formula(paste(
          'Surv(',
          timeCol,
          ', ',
          eventCol,
          ')~',
          variableCol
        ))

      
      data[, variableCol] <- as.factor(data[, variableCol])
      contrasts(data[, variableCol]) <- contr.helmert(4)
      
      
      coxFit <-
        coxph(currentFormula,
              data = data)
      
      if(returnCoxModels){
        return(coxFit)
      }
      else if(!returnCoxModels){
        
        resNames <-
          c(
            "variable",
            colnames(summary(coxFit)$coef),
            colnames(summary(coxFit)$conf.int)[2:4],
            names(summary(coxFit)$concordance),
            names(summary(coxFit)$sctest)
          )
        
        resData <-
          data.frame(matrix(ncol = length(resNames), nrow = 0))
        colnames(resData) <- resNames
        
        
        for (i in rownames(summary(coxFit)$coef)) {
          currentCoxRes <-  c(
            summary(coxFit)$coef[i, ],
            ## first column is exp(coef) which we extracted in the first line
            summary(coxFit)$conf.int[i, 2:4],
            ## C is concordance
            summary(coxFit)$concordance,
            ## p-value is the overall statistical test
            summary(coxFit)$sctest
          )
          
          currentCoxRes <-
            data.frame(variable = i,
                       t(data.frame(currentCoxRes)),
                       check.names = F)
          
          
          resData <- rbind(resData, currentCoxRes)
        }
        
        resData$formula <-
          as.character(currentFormula)[3]
        
        ## add proper contrast names
        if(nStrata == 4){
          resData$contrast.helmert[1] <- 
            paste0(
              row.names(contrasts(data[, variableCol]))[2], " - ", 
              row.names(contrasts(data[, variableCol]))[1])
          
          resData$contrast.helmert[2] <- 
            paste0(
              row.names(contrasts(data[, variableCol]))[3], " - ", 
              row.names(contrasts(data[, variableCol]))[1], " AND ",
              row.names(contrasts(data[, variableCol]))[2]   )
          
          resData$contrast.helmert[3] <- 
            paste0(
              row.names(contrasts(data[, variableCol]))[3], " - ", 
              row.names(contrasts(data[, variableCol]))[1], " AND ",
              row.names(contrasts(data[, variableCol]))[2], " AND ",
              row.names(contrasts(data[, variableCol]))[3]  )
        }
        
        # resData$variable <-
        #   gsub(variableCol, "", resData$variable)
        return(resData)
        
      }
      
    }
  }



##----------------- When we have multiple main variable names to test for: 
coxph_multitest <- function( 
  data = wide2CombExprLong,
  variableNameCol = "Comb_Name", ## column with name of different variables
  variableCol = "Comb_2status", ## column with the groups for each variable
  timeCol = "OS.time",
  eventCol = "OS",
  nStrata = 4,
  returnCovarStat = F,
  covarCol = c("Age", "Stage"),
  nCores = 10,
  returnCoxModels = FALSE, 
  IncludeContrasts = FALSE){
  
  if(! IncludeContrasts){
    cox_res <-
      # lapply(unique(data[, variableNameCol]),
      parallel::mclapply(unique(data[, variableNameCol]),
                         function(x) {
                           coxph_test(
                             data = data[data[, variableNameCol] == x,],
                             timeCol = timeCol,
                             eventCol = eventCol,
                             variableCol = variableCol,
                             nStrata = nStrata,
                             returnCovarStat = returnCovarStat,
                             returnCoxModels = returnCoxModels,
                             covarCol = covarCol
                           )
                           
                         }, mc.cores = nCores)
  }

  if(IncludeContrasts){
    cox_res <-
      # lapply(unique(data[, variableNameCol]),
      parallel::mclapply(unique(data[, variableNameCol]),
                         function(x) {
                           coxph_test_contrasts(
                             data = data[data[, variableNameCol] == x,],
                             timeCol = timeCol,
                             eventCol = eventCol,
                             variableCol = variableCol,
                             nStrata = nStrata,
                             returnCovarStat = returnCovarStat,
                             returnCoxModels = returnCoxModels,
                             covarCol = covarCol
                           )
                           
                         }, mc.cores = nCores)
  }
  

  names(cox_res) <- unique(data[, variableNameCol])
  
  if(returnCoxModels){
    return(cox_res)
  }
  else if(!returnCoxModels){
    cox_res <- do.call(rbind, cox_res)
    
    ## add group names:
    cox_res$VariableName <- sapply(rownames(cox_res), function(y){
      unlist(strsplit(y, "\\."))[1]
    })
    
    return(cox_res)
  }
  # cox_res <- lapply(cox_res, function(y){
  #   y$VariableName <- variableNameCol
  # }) 

}







coxph_test_scores <-
  function(scoreAnnotData = scoreAnnotWide,
           timeCol = "OS.time",
           eventCol = "OS",
           scoreCol = c("Gut_TRM_Up", "TGFb_EMT_Up"),
           covarCol = c("Age_status", "Stage_2group")) {
    
    currentScoreAnnot <-  scoreAnnotData[, c(timeCol,
                                         eventCol,
                                         scoreCol,
                                         covarCol)]
    
    currentScoreAnnot <- currentScoreAnnot[complete.cases(currentScoreAnnot), ]
    
  ## generate Gut_TRM_Up + TGFb_EMT_Up
    scoreInfo <- paste(scoreCol, collapse = " + ")
  ## generate Gut_TRM_Up * TGFb_EMT_Up
    scoreInteraction <- paste(scoreCol, collapse = " * ")
    
  ## add  
    if (! is.null(covarCol)){
      covarInfo <- paste(covarCol, collapse = " + ")
      
      currentFormula <-
        as.formula(
          paste(
            'Surv(',
            timeCol,
            ', ',
            eventCol,
            ')~',
            scoreInfo,
            ' + ',
            scoreInteraction,
            ' + ',
            covarInfo
          )
        )
    }
    else{
      currentFormula <-
        as.formula(
          paste(
            'Surv(',
            timeCol,
            ', ',
            eventCol,
            ')~',
            scoreInfo,
            ' + ',
            scoreInteraction
          )
        )
    }
    
    
    coxFit <-
      coxph(
        currentFormula,
        data = currentScoreAnnot
      )
        
    
    resNames <-
      c(
        "variable",
        colnames(summary(coxFit)$coef),
        colnames(summary(coxFit)$conf.int)[2:4],
        names(summary(coxFit)$concordance),
        names(summary(coxFit)$sctest)
      )
    
    resData <- data.frame(matrix(ncol = length(resNames), nrow = 0))
    colnames(resData) <- resNames
    
    for(i in rownames(summary(coxFit)$coef)){
      currentCoxRes <-  c(
        summary(coxFit)$coef[i,],
        ## first column is exp(coef) which we extracted in the first line
        summary(coxFit)$conf.int[i, 2:4],
        ## C is concordance
        summary(coxFit)$concordance,
        ## p-value is the overall statistical test
        summary(coxFit)$sctest
      )
      
      currentCoxRes <- data.frame(variable = i, t(data.frame(currentCoxRes)), check.names = F)
  
      
      resData <- rbind(resData, currentCoxRes)
      
    }
    resData$formula <- as.character(currentFormula)[3]
    
    return(resData)
  }
                   
                   

  

coxph_Fits <-
  function(scoreAnnot = scoreAnnot,
           timeCol = timeCol,
           eventCol = eventCol,
           sigCol = "Signature",
           scoreCol = scoreCol,
           covarCol = c("Age", "Stage")) {
    coxFits <-
      foreach (i = unique(scoreAnnot[, sigCol]),
               .packages = "survival") %do% {
                 currentScoreAnnot <-  scoreAnnot[scoreAnnot[, sigCol] == i, ]
                 
                 if (length(covarCol) == 1) {
                   coxFit <-
                     coxph(
                       Surv(as.numeric(currentScoreAnnot[, timeCol]), currentScoreAnnot[, eventCol]) ~ 
                         currentScoreAnnot[, scoreCol] + 
                         currentScoreAnnot[, covarCol],
                       data = currentScoreAnnot
                     )
                 }
                 
                 else if (length(covarCol) == 2) {
                   coxFit <-
                     coxph(
                       Surv(as.numeric(currentScoreAnnot[, timeCol]), currentScoreAnnot[, eventCol]) ~ 
                         currentScoreAnnot[, scoreCol] +
                         currentScoreAnnot[, covarCol[1]] + 
                         currentScoreAnnot[, covarCol[2]],
                       data = currentScoreAnnot
                     )
                 }
                 else if (is.null(covarCol)) {
                   coxFit <-
                     coxph(
                       Surv(as.numeric(currentScoreAnnot[, timeCol]), currentScoreAnnot[, eventCol]) ~
                         currentScoreAnnot[, scoreCol],
                       data = currentScoreAnnot
                     )
                 }
                 else if (length(covarCol) > 2) {
                   stop("Please provide maximun of 2 covariate")
                 }

                 return(coxFit)
               }
    
   
    
    names(coxFits) <- unique(scoreAnnot[, sigCol])
    coxFits[["Covariate"]] <- paste(covarCol, collapse = "_")
    
    return(coxFits)
    
  }

