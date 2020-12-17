## This is a script for ploting heatmap of gene expression values using teh ComplexHeatmap package in R.

##--------------- Author: Sepideh Foroutan July 2019 ----------------

##--------------- INPUTS:
## genes : character vector of gene names/IDs to be plotted in teh heatmap
## exprData = expression matrix (in log format) with gene names/IDs in rows and sample names in columns
## subsetSample = a character vector of sample names to subset the data based on them (a subset of column nanmes)

## annotData = an annotation data frame, whose row names are the same as the column names of the expression data
## annotColumns = a character vector of column names in the annotation data to be plotted as sample annotation in the heatmap (a subset of column names in the annotation data)
## annotColList = a named list of colours for the sample annotations (e.g. if we select "Patients" and "Treatment" as annotation columns, then it is a list with two elements called "Patients" and "Treatment" that has the colour values: e.g. annotColList <- list(Patients = c("red", "blue"), Treatment = c("black", "gray"))
## annotLegened_ncol = number of columns for annotation legend

## geneAnnot = data frame that has ONLY the columns that we want to generate gene annotations for (e.g. up and down sets)
## geneColList = named color list for gene annotations
## heatmapCol = color of heatmap; default: viridis::viridis(100)

## geneColList = Colour for gene annotation
## clustRows = T,
## clustCols = T,
## showRowNames = T,
## showColNames = T,
## textSize = 10,
## plotTitle = "Heatmap"

##-------------- OUTPUT:
## A Complext heatmap of scaled gene expression values


##----------------------------------- FUNCTION ---------------------------
plotHeatmapExpr <-
  function(genes,
    exprData = logCPM,
    subsetSample = NULL,
    annotData = NULL,
    annotColumns = c("PatientID", "Tissue", "Response"),
    annotColList = currentCol,
    annotLegened_ncol = 1,
    geneAnnot = NULL,
    geneColList = NULL,
    ## heatmapCol = circlize::colorRamp2(c(minExpr, 0, maxExpr), c("yellow3" , "white" , "navy")),
    # heatmapCol = "default",
    heatmapCol = viridis::viridis(100),
    clustRows = T,
    clustCols = T,
    showRowNames = T,
    showColNames = T,
    scaleData = T,
    heatmapLegendTitle = "Scaled logRPKM",
    textSize = 10,
    plotTitle = "Heatmap") {
    
    ##-------------------------- Subset data for the genes and samples 
    ##                            then scale the expression matrix:
    genes <- as.character(genes)
    gg <- genes[genes %in% rownames(exprData)]
    
    if (length(gg) == 0) {
      stop("Selected genes do not present in the expression data")
    }
    if(! all(rownames(annotData) == colnames(exprData))){
      stop("Make sure annotation data is in the same order as the expression data and has the same row names as the column names in the expression data")
    }


    if(! is.null(subsetSample)){
      annotData <- annotData[subsetSample, , drop = F]
      exprData <- exprData[, subsetSample, drop = F]
    }
    exprData <- exprData[gg, ]
    
    if(scaleData){
      ddgene <- t(scale(t(exprData)))
    }
    else{
      ddgene <- exprData
    }

    
    minVal <- min(ddgene)
    maxVal <- max(ddgene)
    
    if(heatmapCol == "default"){
      heatmapCol <-  circlize::colorRamp2(c(minVal, 0, maxVal), c("yellow3" , "white" , "navy"))
    }
    
    
    if(! is.null(annotData)){
      ss <-
        sapply(names(annotColList), function(x)
          list(
            list(
              legend_direction = "horizontal",
              ncol = annotLegened_ncol,
              title_gp = gpar(fontsize = textSize + 2, fontface = "italic"),
              labels_gp = gpar(fontsize = textSize), 
              grid_height = unit(0.8, "cm")
            )
          ))
      
      ss <- ss[annotColumns]
      
      sampleAnnotHm <-
        HeatmapAnnotation(df = annotData[, annotColumns, drop = F],
                          col = annotColList,
                          annotation_name_side = "left", 
                          annotation_legend_param = ss)
    }

    if(! is.null(geneAnnot)){
      geneAnnot <- geneAnnot[gg, , drop = F]
      geneHmAnn <- HeatmapAnnotation(
        df = geneAnnot,
        # col =   list(Genes = geneCol),
        col =   geneColList,
        which = "row",
        name = "\nGenes",
        annotation_legend_param = list(
          legend_direction = "vertical",
          ncol = 1,
          title_gp = gpar(fontsize = textSize + 2, fontface = "italic"),
          labels_gp = gpar(fontsize = textSize),
          grid_height = unit(0.5, "cm")
        )
      )
    }
      
      # geneHm <- Heatmap(geneAnnot,
      #   name = "Genes",
      #   col = geneCol,
      #   width = unit(5, "mm"), 
      #   show_row_names = showRowNames)
    if(! is.null(annotData) & ! is.null(geneAnnot)){  
      expr_hm <- ComplexHeatmap::Heatmap(
        ddgene,
        col = heatmapCol,
        # name = "logExpr",
        column_title = plotTitle,
        cluster_rows = clustRows,
        cluster_columns = clustCols,
        show_row_names = showRowNames,
        show_column_names = showColNames,
        row_names_gp = gpar(fontsize = textSize),
        column_names_gp = gpar(fontsize = textSize),
        
        # column_title_gp = gpar(fontsize = textSize), 
        heatmap_legend_param = list(
          title_position = "topcenter",
          title = heatmapLegendTitle,
          title_gp = gpar(fontsize = textSize + 2, fontface = "italic"),
          labels_gp = gpar(fontsize = textSize),
          color_bar = "continuous", 
          legend_direction = "vertical",
          # legend_direction = "vertical",
          # legend_width = unit(5, "cm"),
          legend_height = unit(4, "cm"),
          at = c(minVal, 0, maxVal), 
          labels = c(round(minVal, 2), 0, round(maxVal, 2))),
          # legend_height = unit(3, "cm")
        # ),
        top_annotation = sampleAnnotHm,
        right_annotation = geneHmAnn,
        # left_annotation = geneHmAnn,
        na_col = "gray90"
      )
      
    } else if(! is.null(annotData) & is.null(geneAnnot)){
      expr_hm <- ComplexHeatmap::Heatmap(
        ddgene,
        col = heatmapCol, 
        # name = "logExpr",
        column_title = plotTitle,
        cluster_rows = clustRows,
        cluster_columns = clustCols,
        show_row_names = showRowNames,
        show_column_names = showColNames,
        row_names_gp = gpar(fontsize = textSize),
        column_names_gp = gpar(fontsize = textSize),
        # column_title_gp = gpar(fontsize = textSize), 
        heatmap_legend_param = list(
          title_position = "lefttop",
          title = heatmapLegendTitle,
          title_gp = gpar(fontsize = textSize + 2, fontface = "italic"),
          labels_gp = gpar(fontsize = textSize),
          color_bar = "continuous", 
          legend_direction = "horizontal", 
          legend_width = unit(5, "cm"),
          at = c(minVal, 0, maxVal), 
          labels = c(round(minVal, 2), 0, round(maxVal, 2))),
        top_annotation = sampleAnnotHm, 
        na_col = "gray90"
      ) 
      
    }
    else if(is.null(annotData) & ! is.null(geneAnnot)){
      expr_hm <- ComplexHeatmap::Heatmap(
        ddgene,
        col = heatmapCol,
        # name = "logExpr",
        column_title = plotTitle,
        cluster_rows = clustRows,
        cluster_columns = clustCols,
        show_row_names = showRowNames,
        show_column_names = showColNames,
        row_names_gp = gpar(fontsize = textSize),
        column_names_gp = gpar(fontsize = textSize),
        # column_title_gp = gpar(fontsize = textSize), 
        heatmap_legend_param = list(
          title_position = "topcenter",
          title = paste0(heatmapLegendTitle, "\n"),
          title_gp = gpar(fontsize = textSize + 2, fontface = "italic"),
          labels_gp = gpar(fontsize = textSize),
          color_bar = "continuous", 
          legend_direction = "vertical", 
          # legend_width = unit(5, "cm"),
          legend_height = unit(4, "cm"),
          at = c(minVal, 0, maxVal), 
          labels = c(round(minVal, 2), 0, round(maxVal, 2))),
        # top_annotation = sampleAnnotHm, 
        right_annotation = geneHmAnn,
        na_col = "gray90"
      ) 
      
    }
    else if( is.null(annotData) & is.null(geneAnnot)){
      expr_hm <- ComplexHeatmap::Heatmap(
        ddgene,
        col = heatmapCol,
        # name = "logExpr",
        column_title = plotTitle,
        cluster_rows = clustRows,
        cluster_columns = clustCols,
        show_row_names = showRowNames,
        show_column_names = showColNames,
        row_names_gp = gpar(fontsize = textSize),
        column_names_gp = gpar(fontsize = textSize),
        # column_title_gp = gpar(fontsize = textSize), 
        heatmap_legend_param = list(
          title_position = "lefttop",
          title = heatmapLegendTitle, 
          title_gp = gpar(fontsize = textSize + 2, fontface = "italic"),
          labels_gp = gpar(fontsize = textSize),
          color_bar = "continuous", 
          legend_direction = "horizontal", 
          legend_width = unit(5, "cm"),
          at = c(minVal, 0, maxVal), 
          labels = c(round(minVal, 2), 0, round(maxVal, 2))),
        # top_annotation = sampleAnnotHm, 
        na_col = "gray90"
      ) 
    }
    
    ComplexHeatmap::draw(
      expr_hm,
      heatmap_legend_side = "bottom",
      annotation_legend_side = "right",
      legend_title_gp = gpar(fontsize = textSize, fontface = "italic")
    )
    
  }
