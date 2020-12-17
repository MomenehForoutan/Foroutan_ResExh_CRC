


visReduction_feature <-
  function(object = tt,
    reduction = "umap",
    dim1Name = "UMAP_1",
    dim2Name = "UMAP_2",
    objAssay = "SCT",
    objSlot = "data",
    # dataType = "data",
    nrow = 2,
    ncol = 3,
    features = gg,
    textSize = 12,
    pointSize = 0.3,
    pointAlpha = 1,
    plotTitle = "UMAP plot",
    plotHexbin = FALSE,
    actionGene = "median") {
   
     # ss <- object
    
    features <- features[features  %in% rownames(object)]
    if (!length(features) > 0) {
      stop("Genes do not exist in the data")
    }
    
    
    currentRedData <- object@reductions[[reduction]]@cell.embeddings
    # currentExpr <- object[[assaySlot]]@data
    # currentExpr <- as.matrix(GetAssayData(object))
    comFeatures <- unique(as.character(intersect(features, rownames(object))))
    object <- object[comFeatures, ]
    
    currentExpr <- as.matrix(GetAssayData(object[[objAssay]], slot = objSlot))
    currentExpr <-
      currentExpr[rownames(currentExpr) %in% features, , drop = F]
    
    currentDataGenes <-
      data.frame(cbind(currentRedData, t(currentExpr)))
    
    if (plotHexbin) {
      plotList <-  lapply(features, function(x) {
        maxExpr <- max(currentExpr[x, ])
        plot_hexbin_gene(object,
          type = objSlot,
          gene = x,
          action = actionGene) +
          scale_color_manual(values = cols) +
          ggtitle(paste0(x, " (max = ", round(maxExpr, 2), ")")) +
          theme_dark() +
          theme(
            legend.position = "none",
            axis.title = element_blank(),
            axis.text = element_blank(),
            strip.text = element_text(size = textSize)
          ) +
          guides(fill = guide_legend(
            title = annot,
            override.aes = list(size = textSize - 10, stroke = 1.2)
          ))
        
      })
      
      glist <- lapply(plotList, ggplotGrob)
      
      ggpubr::ggarrange(plotlist = glist,
        nrow = nrow,
        ncol = ncol)
    }
    
    if (!plotHexbin) {
      
      ## Reformat the data from wide to long so that we can use ggplot
      currentDataGenesLong <-
        tidyr::gather(currentDataGenes,
          Genes,
          LogExpr,
          3:ncol(currentDataGenes))
      
      ## Reorder cells based on expression and add max values for gene expression in the title
      currentDataGenesLongArrange <- currentDataGenesLong %>%
        group_by(Genes) %>%
        arrange(LogExpr) %>%
        mutate(maxVal = round(max(LogExpr), 2)) %>%
        mutate(minVal = 0) %>%
        mutate(GenesMax = paste0(Genes, " (max = ", maxVal, ")")) %>%
        ungroup() %>%
        data.frame()
      
      
      p <- currentDataGenesLongArrange %>% group_by(GenesMax) %>%
        do(
          plots = ggplot(data = .) + aes_string(dim1Name , dim2Name, colour = "LogExpr") +
            geom_point(size = pointSize, alpha = pointAlpha) +
            facet_wrap(~ GenesMax) +
            # scale_color_gradient(low = "darkorchid4", high = "yellow") +
            scale_color_viridis_c(100) +
            # scale_color_gradient(low = "gray80", high = "darkblue") +
            ggtitle(plotTitle) +
            theme_dark() +
            # theme_bw() +
            theme(
              legend.position = "none",
              axis.title = element_blank(),
              axis.text = element_blank(),
              strip.text = element_text(size = textSize )
            )
        )
      # return(p)
      glist <- lapply(p$plots, ggplotGrob)
      
    }
    
    # gridExtra::marrangeGrob(glist, nrow = nrow, ncol = ncol)
    
    multiPage <-
      ggpubr::ggarrange(plotlist = glist,
        nrow = nrow,
        ncol = ncol)
    return(multiPage)
    
  }





visReduction_meta <-
  function(object = tt,
    reduction = "umap",
    dim1Name = "UMAP_1",
    dim2Name = "UMAP_2",
    annot = "Non.malignant",
    textSize = 1.6,
    pointSize = 0.3,
    pointAlpha = 1,
    plotTitle = "UMAP plot",
    plotHexbin = FALSE,
    actionMeta = "majority", 
    cols = brewer.pal(9, "Set1")[c(2, 8, 5, 4, 3, 9)]) {
    ##----- if we have annotation file:
    # ss <- object
    if (!annot  %in%  colnames(object@meta.data)) {
      stop(
        "The selected name does not exist in the meta-data columns; check object@metadata for available columns"
      )
    }
    
    currentRedData <- object@reductions[[reduction]]@cell.embeddings
    
    theme_base <-
      theme_bw() +
      # theme_dark() +
      theme(
        title = element_text(size = textSize * 4),
        # legend.position = c(0.72, 0.25),
        legend.title = element_text(size = textSize * 5, face = "italic"),
        legend.text = element_text(size = textSize * 4),
        axis.text.x = element_text(size = textSize * 4),
        axis.text.y = element_text(size = textSize * 4),
        axis.title = element_text(size = textSize * 5)
      )
    
    if (plotHexbin) {
      p <-
        plot_hexbin_meta(object,
          annot,
          action = actionMeta,
          xlab = dim1Name,
          ylab = dim2Name) +
        scale_fill_manual(values = cols) +
        ggtitle(plotTitle) +
        theme_base +
        guides(color = guide_legend(
          title = annot,
          override.aes = list(size = textSize + 0.5, stroke = 1.2)
        ))
    }
    else if (!plotHexbin) {
      currentAnnot <- object@meta.data[[annot]]
      
      currentDataAnnot <-
        data.frame(annot = currentAnnot, currentRedData)
      colnames(currentDataAnnot)[1] <- annot
      
      p <-
        ggplot(currentDataAnnot,
          aes_string(dim1Name , dim2Name, colour = annot)) +
        geom_point(size = (pointSize), alpha = pointAlpha) +
        scale_color_manual(values = cols) +
        ggtitle(plotTitle) +
        theme_base +
        guides(color = guide_legend(
          title = annot,
          override.aes = list(size = textSize + 0.5, stroke = 1.2)
        ))
      
    }
    
    return(p)
  }
    









# visReduction <-
#   function(object = tt,
#     # assaySlot = "logTPM",
#     reduction = "umap",
#     dim1Name = "UMAP_1",
#     dim2Name = "UMAP_2",
#     nrow = 2,
#     ncol = 3,
#     features = gg,
#     # annot = "Non.malignant",
#     annot = NULL,
#     # cols = brewer.pal(8, "Set1")[c(5, 3, 1, 2)],
#     cols = brewer.pal(9, "Set1")[c(2, 8, 5, 4, 3, 9)],
#     # cols = brewer.pal(8, "Set1")[c(2, 5, 3, 1)]
#     textSize = 1.6,
#     pointSize = 0.3,
#     pointAlpha = 1,
#     plotTitle = "UMAP plot",
#     plotHexbin = FALSE,
#     actionMeta = "majority",
#     actionGene = "median") {
#     
#     ss <- object
#     
#     if (!is.null(features)) {
#       features <- features[features  %in% rownames(ss)]
#       if (!length(features) > 0) {
#         stop("Genes do not exist in the data")
#       }
#     }
#     
#     
#     currentRedData <- ss@reductions[[reduction]]@cell.embeddings
#     # currentExpr <- ss[[assaySlot]]@data
#     currentExpr <- as.matrix(GetAssayData(ss))
#     currentExpr <-
#       currentExpr[rownames(currentExpr) %in% features, ]
#     
#     
#     currentDataGenes <-
#       data.frame(cbind(currentRedData, t(currentExpr)))
#     
#     ##----- if we have annotation file:
#     if (!is.null(annot)) {
#       theme_base <-
#         theme_bw() +
#         # theme_dark() +
#         theme(
#           title = element_text(size = textSize * 4),
#           # legend.position = c(0.72, 0.25),
#           legend.title = element_text(size = textSize * 5, face = "italic"),
#           legend.text = element_text(size = textSize * 4),
#           axis.text.x = element_text(size = textSize * 4),
#           axis.text.y = element_text(size = textSize * 4),
#           axis.title = element_text(size = textSize * 5)
#         )
#       
#       if (plotHexbin) {
#         p.1 <-
#           plot_hexbin_meta(tt,
#             annot,
#             action = actionMeta,
#             xlab = dim1Name,
#             ylab = dim2Name) +
#           scale_fill_manual(values = cols) +
#           ggtitle(plotTitle) +
#           theme_base +
#           guides(color = guide_legend(
#             title = annot,
#             override.aes = list(size = textSize + 0.5, stroke = 1.2)
#           ))
#       }
#       else if (!plotHexbin) {
#         currentAnnot <- ss@meta.data[[annot]]
#         
#         currentDataGenesAnnot <-
#           data.frame(cbind(currentAnnot, currentDataGenes))
#         colnames(currentDataGenesAnnot)[1] <- annot
#         
#         currentDataGenesLong <-
#           tidyr::gather(currentDataGenesAnnot,
#             Genes,
#             LogExpr,
#             4:ncol(currentDataGenesAnnot))
#         
#         p.1 <-
#           ggplot(currentDataGenesLong,
#             aes_string(dim1Name , dim2Name, colour = annot)) +
#           geom_point(size = (pointSize - 0.2), alpha = pointAlpha) +
#           scale_color_manual(values = cols) +
#           ggtitle(plotTitle) +
#           theme_base +
#           guides(color = guide_legend(
#             title = annot,
#             override.aes = list(size = textSize + 0.5, stroke = 1.2)
#           ))
#         
#       }
#       
#       return(p.1)
#     }
#     else{
#       if (plotHexbin) {
#         plotList <-  lapply(features, function(x) {
#           maxExpr <- max(currentExpr[x,])
#           plot_hexbin_gene(tt,
#             type = "data",
#             gene = x,
#             action = actionGene) +
#             scale_color_manual(values = cols) +
#             ggtitle(paste0(x, " (max = ", round(maxExpr, 2), ")")) +
#             theme_dark() +
#             theme(
#               legend.position = "none",
#               axis.title = element_blank(),
#               axis.text = element_blank(),
#               strip.text = element_text(size = textSize * 8)
#             ) +
#             guides(fill = guide_legend(
#               title = annot,
#               override.aes = list(size = textSize + 0.5, stroke = 1.2)
#             ))
#           
#         })
#         
#         glist <- lapply(plotList, ggplotGrob)
#         
#         ggpubr::ggarrange(plotlist = glist,
#           nrow = nrow,
#           ncol = ncol)
#         
#         
#         # gridExtra::marrangeGrob(glist, nrow = nrow, ncol = ncol)
#         
#       }
#       else if (!plotHexbin) {
#         currentDataGenesAnnot <-
#           data.frame(cbind(rep(NA, nrow(
#             currentDataGenes
#           )), currentDataGenes))
#         
#         ## Reformat the data from wide to long so that we can use ggplot
#         currentDataGenesLong <-
#           tidyr::gather(currentDataGenesAnnot,
#             Genes,
#             LogExpr,
#             4:ncol(currentDataGenesAnnot))
#         
#         ## Reorder cells based on expression and add max values for gene expression in the title
#         currentDataGenesLongArrange <- currentDataGenesLong %>%
#           group_by(Genes) %>%
#           arrange(LogExpr) %>%
#           mutate(maxVal = round(max(LogExpr), 2)) %>%
#           mutate(minVal = 0) %>%
#           mutate(GenesMax = paste0(Genes, " (max = ", maxVal, ")")) %>%
#           ungroup() %>%
#           data.frame()
#         
#         
#         p <- currentDataGenesLongArrange %>% group_by(GenesMax) %>%
#           do(
#             plots = ggplot(data = .) + aes_string(dim1Name , dim2Name, colour = "LogExpr") +
#               geom_point(size = pointSize) +
#               facet_wrap( ~ GenesMax, ncol = 4) +
#               scale_color_gradient(low = "darkorchid4", high = "yellow") +
#               ggtitle(plotTitle) +
#               theme_dark() +
#               theme(
#                 legend.position = "none",
#                 axis.title = element_blank(),
#                 axis.text = element_blank(),
#                 strip.text = element_text(size = textSize * 8)
#               )
#           )
#         # return(p)
#         glist <- lapply(p$plots, ggplotGrob)
#         
#       }
#       
#       # gridExtra::marrangeGrob(glist, nrow = nrow, ncol = ncol)
#       
#       multiPage <-
#         ggpubr::ggarrange(plotlist = glist,
#           nrow = nrow,
#           ncol = ncol)
#       return(multiPage)
#       
#     }
#     
#   }
# 
# 
# # gridExtra::grid.arrange(
# #   p$plots[[1]],
# #   p$plots[[2]],
# #   p$plots[[3]],
# #   p$plots[[4]],
# #   p$plots[[5]],
# #   p$plots[[6]],
# #   # p$plots[[7]],
# #   # p$plots[[8]],
# #   nrow = nrow
# # )
