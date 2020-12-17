library(survival)
library(foreach)

coxph_test <-
  function(scoreAnnot = scoreAnnot,
           timeCol = timeCol,
           eventCol = eventCol,
           sigCol = "Signature",
           scoreCol = scoreCol,
           covarCol = c("Age", "Stage")) {
    coxStats <-
      foreach (i = unique(scoreAnnot[, sigCol]),
               .packages = "survival",
               .combine = "rbind") %do% {
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
                 # else (length(covarCol) > 1){
                 #   SurvFormula <- as.formula(
                 #     paste("Surv(timeCol, eventCol ~ scoreCol +",
                 #           paste(covarCol, collapse = "+"), 
                 #           ", data = get(currentScoreAnnot)"))
                 # 
                 #   as.formula(paste("Surv(t1, e) ~ ",
                 #                    paste(c2, collapse= "+")))
                 # }
                 
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
                 
                 ##--- export only the results for the scores levels (low vs high)
                 currentCoxRes <-  c(
                   summary(coxFit)$coef[1,],
                   ## first column is exp(coef) which we extracted in the first line
                   summary(coxFit)$conf.int[1, 2:4],
                   ## C is concordance
                   summary(coxFit)$concordance,
                   ## p-value is the overall statistical test
                   summary(coxFit)$sctest
                 )
                 return(currentCoxRes)
               }
    
    coxStats <- data.frame(coxStats, check.names = F)
    
    coxStats$Signature <- unique(scoreAnnot[, sigCol])
    coxStats$Covariates <- paste(covarCol, collapse = "_")
    
    return(coxStats)
    
  }




##------ added on the 
## commented on Oct 2020

## data has columns for categorical score groups and 

# coxph_test <-
#   function(data = scoreAnnot,
#            timeCol = timeCol,
#            eventCol = eventCol,
#            # sigCol = "Signature",
#            variableCol = scoreCol,
#            nStrata = 4,
#            returnCovarStat = F,
#            covarCol = c("Age", "Stage")) {
#     
#     # coxStats <-
#     #   foreach (i = unique(data[, sigCol]),
#     #            .packages = "survival",
#     #            .combine = "rbind") %do% {
#                  
#     # currentData <-  data[data[, sigCol] == i,]
#     
#                  
#                  if (!is.null(covarCol)) {
#                    covarInfo <- paste(covarCol, collapse = " + ")
#                    
#                    currentFormula <-
#                      as.formula(paste(
#                        'Surv(',
#                        timeCol,
#                        ', ',
#                        eventCol,
#                        ')~',
#                        variableCol,
#                        ' + ',
#                        covarInfo
#                      ))
#                    
#                    coxFit <-
#                      coxph(currentFormula,
#                            data = data)
#                    
#                    
#                    resNames <-
#                      c(
#                        "variable",
#                        colnames(summary(coxFit)$coef),
#                        colnames(summary(coxFit)$conf.int)[2:4],
#                        names(summary(coxFit)$concordance),
#                        names(summary(coxFit)$sctest)
#                      )
#                    
#                    resData <-
#                      data.frame(matrix(ncol = length(resNames), nrow = 0))
#                    colnames(resData) <- resNames
#                    
#                    
#                    ## Get stats for all variables including covariates
#                    if (returnCovarStat) {
#                      for (i in rownames(summary(coxFit)$coef)) {
#                        currentCoxRes <-  c(
#                          summary(coxFit)$coef[i, ],
#                          ## first column is exp(coef) which we extracted in the first line
#                          summary(coxFit)$conf.int[i, 2:4],
#                          ## C is concordance
#                          summary(coxFit)$concordance,
#                          ## p-value is the overall statistical test
#                          summary(coxFit)$sctest
#                        )
#                        
#                        currentCoxRes <-
#                          data.frame(variable = i,
#                                     t(data.frame(currentCoxRes)),
#                                     check.names = F)
#                        
#                        
#                        resData <- rbind(resData, currentCoxRes)
#                      }
#                      resData$formula <-
#                        as.character(currentFormula)[3]
#                    }
#                      ## Get three first line of the results as we have 4 strata
#                      if (!returnCovarStat) {
#                        for (i in rownames(summary(coxFit)$coef)[1:(nStrata - 1)]) {
#                          currentCoxRes <-  c(
#                            summary(coxFit)$coef[i, ],
#                            ## first column is exp(coef) which we extracted in the first line
#                            summary(coxFit)$conf.int[i, 2:4],
#                            ## C is concordance
#                            summary(coxFit)$concordance,
#                            ## p-value is the overall statistical test
#                            summary(coxFit)$sctest
#                          )
#                          
#                          currentCoxRes <-
#                            data.frame(variable = i,
#                                       t(data.frame(currentCoxRes)),
#                                       check.names = F)
#                          
#                          
#                          resData <- rbind(resData, currentCoxRes)
#                        }
#                        
#                        resData$formula <-
#                          as.character(currentFormula)[3]
#                      }
#                    resData$variable <-
#                      gsub(variableCol, "", resData$variable)
#                    return(resData)
#                    
#                  }
#     
#     
#     if (is.null(covarCol)) {
#       
#       currentFormula <-
#         as.formula(paste(
#           'Surv(',
#           timeCol,
#           ', ',
#           eventCol,
#           ')~',
#           variableCol
#         ))
#       
#       coxFit <-
#         coxph(currentFormula,
#               data = data)
#       
#       
#       resNames <-
#         c(
#           "variable",
#           colnames(summary(coxFit)$coef),
#           colnames(summary(coxFit)$conf.int)[2:4],
#           names(summary(coxFit)$concordance),
#           names(summary(coxFit)$sctest)
#         )
#       
#       resData <-
#         data.frame(matrix(ncol = length(resNames), nrow = 0))
#       colnames(resData) <- resNames
#       
# 
#         for (i in rownames(summary(coxFit)$coef)) {
#           currentCoxRes <-  c(
#             summary(coxFit)$coef[i, ],
#             ## first column is exp(coef) which we extracted in the first line
#             summary(coxFit)$conf.int[i, 2:4],
#             ## C is concordance
#             summary(coxFit)$concordance,
#             ## p-value is the overall statistical test
#             summary(coxFit)$sctest
#           )
#           
#           currentCoxRes <-
#             data.frame(variable = i,
#                        t(data.frame(currentCoxRes)),
#                        check.names = F)
#           
#           
#           resData <- rbind(resData, currentCoxRes)
#         }
#         
#         resData$formula <-
#           as.character(currentFormula)[3]
#         
#       resData$variable <-
#         gsub(variableCol, "", resData$variable)
#       return(resData)
#       
#     }
#   }












##----- the below code tried to return the cox models

coxph_test <-
  function(data = scoreAnnot,
           timeCol = timeCol,
           eventCol = eventCol,
           # sigCol = "Signature",
           variableCol = scoreCol,
           nStrata = 4,
           returnCovarStat = F,
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

##----------------- When we have multiple variable names to test for: 
coxph_multitest <- function( 
  data = wide2CombExprLong,
  variableNameCol = "Comb_Name",
  variableCol = "Comb_2status",
  timeCol = "OS.time",
  eventCol = "OS",
  nStrata = 4,
  returnCovarStat = F,
  covarCol = c("Age", "Stage"),
  nCores = 10,
  returnCoxModels = FALSE){
  
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




# l <- for( x in unique(data[, variableNameCol])){
#   coxph_test(
#     data = data[data[, variableNameCol] == x,],
#     timeCol = timeCol,
#     eventCol = eventCol,
#     variableCol = variableCol,
#     nStrata = nStrata,
#     returnCovarStat = returnCovarStat,
#     covarCol = covarCol
#   )
# }

                   
                   
                   # if (is.null(covarCol)) {
                   #   currentFormula <-
                   #     as.formula(paste('Surv(',
                   #                      timeCol,
                   #                      ', ',
                   #                      eventCol,
                   #                      ')~',
                   #                      scoreCol))
                   #   
                   #   coxFit <-
                   #     coxph(currentFormula,
                   #           data = currentScoreAnnot)
                   #   
                   #   
                   #   resNames <-
                   #     c(
                   #       "variable",
                   #       colnames(summary(coxFit)$coef),
                   #       colnames(summary(coxFit)$conf.int)[2:4],
                   #       names(summary(coxFit)$concordance),
                   #       names(summary(coxFit)$sctest)
                   #     )
                   #   
                   #   resData <-
                   #     data.frame(matrix(ncol = length(resNames), nrow = 0))
                   #   colnames(resData) <- resNames
                   #   
                   #   
                   #   ## Get stats for all variables including covariates
                   #   
                   #   for (i in rownames(summary(coxFit)$coef)) {
                   #     currentCoxRes <-  c(
                   #       summary(coxFit)$coef[i, ],
                   #       ## first column is exp(coef) which we extracted in the first line
                   #       summary(coxFit)$conf.int[i, 2:4],
                   #       ## C is concordance
                   #       summary(coxFit)$concordance,
                   #       ## p-value is the overall statistical test
                   #       summary(coxFit)$sctest
                   #     )
                   #   }
                   #   currentCoxRes <-
                   #     data.frame(variable = i,
                   #                t(data.frame(currentCoxRes)),
                   #                check.names = F)
                   #   
                   #   
                   #   resData <- rbind(resData, currentCoxRes)
                   #   
                   # }
                   # resData$formula <-
                   #   as.character(currentFormula)[3]
                   # resData$variable <-
                   #   gsub(scoreCol, "", resData$variable)
                   # return(resData)
                 # }
                 
  #               
  #                }
  #              }
  #   return(coxStats)
  # }
  #                      
                     
               
                       
        





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
                   
                   

  
# 
# covariates <- c( "Age", "Age_status", "Stage", "Stage_2group", "bmi", "race", "MSI_Group", "Score_2status")
# 
# univ_formulas <- sapply(covariates,
#                         function(x) as.formula(paste('Surv(time, status)~', x)))
# 
# univ_models <- lapply( univ_formulas, function(x){coxph(x, data = lung)})
# # Extract data 
# univ_results <- lapply(univ_models,
#                        function(x){ 
#                          x <- summary(x)
#                          p.value<-signif(x$wald["pvalue"], digits=2)
#                          wald.test<-signif(x$wald["test"], digits=2)
#                          beta<-signif(x$coef[1], digits=2);#coeficient beta
#                          HR <-signif(x$coef[2], digits=2);#exp(beta)
#                          HR.confint.lower <- signif(x$conf.int[,"lower .95"], 2)
#                          HR.confint.upper <- signif(x$conf.int[,"upper .95"],2)
#                          HR <- paste0(HR, " (", 
#                                       HR.confint.lower, "-", HR.confint.upper, ")")
#                          res<-c(beta, HR, wald.test, p.value)
#                          names(res)<-c("beta", "HR (95% CI for HR)", "wald.test", 
#                                        "p.value")
#                          return(res)
#                          #return(exp(cbind(coef(x),confint(x))))
#                        })
# res <- t(as.data.frame(univ_results, check.names = FALSE))
# as.data.frame(res)


# coxph_univar <-
#   function(data = scoreAnnot[scoreAnnot$Signature == "MSI", ],
#            timeCol = timeCol,
#            eventCol = eventCol,
#            sigName = "MSI",
#            covarCol = covariates) {
#     
#     univ_formulas <- sapply(covariates,
#                             function(x) as.formula(paste('Surv(', timeCol, ', ', eventCol, ')~', x)))
#     # univ_formulas <- sapply(covariates,
#     #                         function(x) as.formula(paste('Surv(time, status)~', x)))
#     
#     univ_models <- lapply( univ_formulas, function(x){coxph(x, data = data)})
#     
#     univ_results <- lapply(univ_models,
#                            function(x) {
#                              x <- summary(x)
#                              p.value <-
#                                signif(x$wald["pvalue"], digits = 2)
#                              wald.test <-
#                                signif(x$wald["test"], digits = 2)
#                              beta <-
#                                signif(x$coef[1], digits = 2)
#                              #coeficient beta
#                              HR <-
#                                signif(x$coef[2], digits = 2)
#                              #exp(beta)
#                              HR.CI.lower <-
#                                signif(x$conf.int[, "lower .95"], 2)
#                              HR.CI.upper <-
#                                signif(x$conf.int[, "upper .95"], 2)
#                              HR_CIs <- paste0(HR, " (",
#                                               HR.CI.lower, "-", HR.CI.upper, ")")
#                              # Signature <- sigName
#                              res <- c(beta, HR, HR_CIs, HR.CI.lower, HR.CI.upper, wald.test, p.value)
#                              names(res) <-
#                                c("beta", "HR", "HR (95% CI for HR)", 
#                                  "HR.CI.lower", "HR.CI.upper", "wald.test",
#                                  "p.value")
#                              return(res)
#                              #return(exp(cbind(coef(x),confint(x))))
#                            })
#     res <- t(as.data.frame(univ_results, check.names = FALSE))
#     as.data.frame(res)
#     
#     return(res)
#   }    
#     


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


# Plots of scaled Schoenfeld residuals against transformed time for each covariate in a model
# fit to the recidivism data. The solid line is a smoothing spline fit to the plot, with the broken lines
# representing a ± 2-standard-error band around the fit.

# ll <- coxph_Fits(scoreAnnot = scoreAnnot,
#            timeCol = timeCol,
#            eventCol = eventCol,
#            sigCol = "Signature",
#            scoreCol = scoreCol,
#            covarCol = c("Age", "Stage"))
# 
# 
# plot(cox.zph(ll$B.CELL))
# plot(cox.zph(ll$EpiB))
# plot(cox.zph(ll$MSI))



# survData = annScoreClin2[ ! annScoreClin2$Score_status == "Intermediate", ]
# timeCol = "days_to_last_follow_up"
# eventCol = "Overall_Survival_Status"
# covarCol = colnames(annScoreClin2)[10:ncol(annScoreClin2)-1]
# groupCol = "Score_status"
# 
# 
# coxStats = foreach(dr = 1:length(unique(survData$drugs)), .packages="survival", .combine = "rbind") %do% {
#   drName = unique(survData$drugs)[dr]
#   currentData = survData[survData$drugs == drName, ]
#   
#   coxFit = coxph(Surv(as.numeric(currentData[,timeCol]), currentData[,eventCol]) ~ currentData[, groupCol] + 
#       # sample_type +
#       # Contribution +                                         
#       tumor_stage + 
#       race +                                               
#       ethnicity +
#       histological_type +                                  
#       PAM50Call_RNAseq +                                  
#       lab_proc_her2_neu_immunohistochemistry_receptor_status +
#       breast_carcinoma_estrogen_receptor_status +         
#       breast_carcinoma_progesterone_receptor_status +   
#       es.TotalScore  +                                    
#       ms.TotalScore +                                     
#       ts.TotalScore +                                    
#       age_at_initial_pathologic_diagnosis   ,
#     data = currentData)
#   
#   ##----------- If I want to export only the results for the drug scores levels (low vs high)
#   coxRes = c(summary(coxFit)$coef[1,],
#     summary(coxFit)$conf.int[1,2:4],
#     summary(coxFit)$concordance,
#     summary(coxFit)$sctest)
#   
#   coxRes3 = coxRes
#   
#   ##----------- If I want to export everything which is significant:
#   
#   # currentCoef = summary(coxFit)$coef
#   # currentCI = summary(coxFit)$conf.int[, 2:4]
#   #   coxRes1 = cbind(currentCoef, 
#   #              currentCI)
#   #   coxRes1 = coxRes1[complete.cases(coxRes1), ]
#   #   coxRes2 = cbind(coxRes1, 
#   #                   C = summary(coxFit)$concordance[1], 
#   #                   `se(C)` = summary(coxFit)$concordance[2], 
#   #                   test = summary(coxFit)$sctest[1],
#   #                   df = summary(coxFit)$sctest[2],
#   #                   pvalue = summary(coxFit)$sctest[3])
#   #   
#   #   ## filter for significant associations:
#   #   coxRes3 = coxRes2[ coxRes2[, grepl("Pr(>|z|)", colnames(coxRes2))] < 0.05, ] 
#   #   coxRes3 = data.frame(coxRes3, check.names = F)
#   #   if(nrow(coxRes3)< 1){
#   #     coxRes3[1, ] = NA
#   #   } 
#   #   coxRes3$drug = drName
#   
#   
#   return(coxRes3)
#   
# }
# 
# ##----------- If I want to export only the results for the drug scores levels (low vs high)
# coxStats = data.frame(coxStats, check.names = F)
# coxStats$drugs = unique(survData$drugs)
# write.csv(coxStats, paste0(outTest, "Survival_Cox_TCGA_drugScores_HighLow.csv"), row.names = F)
# 
# 
# 
# 
# 
# 
# currentClin = clinER[, c("OS_dmfs", "OSbin_dmfs", "ER_final", "Age")]
# row.names(currentClin) = clinER$SampleID
# 
# 
# ## -- exprData: genes in cols and sample names in rows
# ## -- survData: sample names in rows, and survival and covariate information in cols
# ## -- timeCol: column name that has time for survival
# ## -- eventCol: column name that has event status for survival
# ## -- column names for teh covariates to be included
# 
# survivalTest = function(exprData = as.matrix(uncorExpr[, -6]),
#                         survData = currentClin,
#                         timeCol = "OS_dmfs",
#                         eventCol = "OSbin_dmfs",
#                         covarCol = c("Age", "ER_final"),
#                         plotKM = F){
#   
#   ## remove samples without annotation:
#   survData = survData[complete.cases(survData), ]
#   
#   ## subset expression based on those samples that have OS annotation:
#   exprData = exprData[row.names(survData), ]
# 
#   
#   ## calculate median values for all genes and add that to the expression data:
#   medExpr = apply(exprData, 2, median)
#   exprData = rbind(exprData , medExpr)
#   
#   ## replace expression values using
#   newData = apply(exprData, 2, function(x){
#     sapply(x[-length(x)], function(y){
#       ifelse(y > x[length(x)], "High", "Low")
#     })
#   })
#   
#   row.names(newData) = row.names(exprData[- nrow(exprData),])
#   newData = data.frame(newData, check.names = F)
# 
#   
#   coxStats = foreach(g = 1:ncol(newData), .packages="survival", .combine = "rbind") %do% {
#     
#     coxFit = coxph(Surv(survData[,timeCol], survData[,eventCol]) ~ newData[,g] +
#                      survData[,covarCol[1]] +
#                      survData[,covarCol[2]],
#                    data = survData)
#     
#     coxRes = c(summary(coxFit)$coef[1,],
#                summary(coxFit)$conf.int[1,2:4],
#                summary(coxFit)$concordance,
#                summary(coxFit)$sctest)
#     return(coxRes)
#     
#   }
#   row.names(coxStats) = colnames(newData)
#   
#   return(as.data.frame(coxStats))
# }
# 
# 
# ##------------------------------------------------------- Example:
# 
# uncorStats = survivalTest(exprData = mat,
#                           survData = currentClin,
#                           timeCol = "OS_dmfs",
#                           eventCol = "OSbin_dmfs",
#                           covarCol = c("Age", "ER_final"))
# 
# 
# If I do something similar:
# 
#  
# # Columns/variables that might influence the survial
#  
# covCol = colnames(mydata)[c(6,11:19,22,23)]
# sapply(covCol, function(x){
#  cov = paste(x, collapse = "+")
#  coxFit = coxph(formula = as.formula(paste("mySurv~pred+Age+Gender",cov,sep="+")), data = mydata)
#      
#     coxRes = c(summary(coxFit)$coef[1,],
#                summary(coxFit)$conf.int[1,2:4],
#                summary(coxFit)$concordance,
#                summary(coxFit)$sctest)
#     coxRes
# })
# 
# Therefore I got some measurements for each new column, I'm not sure how to interpret this though. Perhaps discuss this next time we meet?
