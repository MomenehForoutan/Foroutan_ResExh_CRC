### The code below have been taken from online sources
## See https://github.com/kstreet13/slingshot/issues/73

## Point on curve function ----
points_on_curve <- function(curve, lambda, ...) {
  UseMethod("points_on_curve", curve)
}

points_on_curve.principal_curve <- function(curve, lambda, ...) {
  if (nrow(curve$s) == length(curve$lambda)) { # didn't use approx_points
    S <- apply(curve$s, 2, function(sjj) {
      return(approx(
        x = curve$lambda[curve$ord],
        y = sjj[curve$ord],
        xout = lambda, ties = "ordered"
      )$y)
    })
  } else {
    if (all(curve$ord == seq_along(curve$lambda))) { # used approx_points
      curvelambda <- seq(min(curve$lambda), max(curve$lambda), length.out = nrow(curve$s))
      S <- apply(curve$s, 2, function(sjj) {
        return(approx(
          x = curvelambda,
          y = sjj,
          xout = lambda, ties = "ordered"
        )$y)
      })
    }
  }
  return(S)
}

points_on_curve.SlingshotDataSet <- function(curve, lambda, ...) {
  locs <- lapply(slingCurves(curve), function(crv) {
    points_on_curve(crv, lambda, ...)
  })
  locs <- do.call('rbind', locs)
  colnames(locs) <- paste0("dim", seq_len(ncol(locs)))
  return(as.data.frame(locs))
}

### Extend ggplot function

#' Plot the gene in reduced dimension space
#'
#' @param sds The output from a lineage computation
#' @param col The assignation of each cell to a label. If none is provided, 
#' default to the cluster labels provided when creating the \code{\link{SlingshotDataSet}}
#' @param title Title for the plot.
#' @param lineSize Size of the curve lineages. Default to 1.
#' @param ... Other options passed to \code{\link{geom_point}}
#' @return A \code{\link{ggplot}} object
#' @examples
#' data('slingshotExample')
#' sds <- slingshot(rd, cl)
#' gg_plot(sds)
#' 
#' ## Change point size and transparency
#' gg_plot(sds, size = 2, alpha = .5)
#' 
#' ## Use grey background
#' 
#' gg_plot(sds) + theme_grey()
#' 
#' ## Color by gene expression
#' gene_count <- sample(0:10, nrow(reducedDims(sds)), replace = TRUE) 
#' gg_plot(sds, col = gene_count)
#' 
#' ## Add a marker of pseudotime
#' gg_plot(sds) + geom_point(data = points_on_curve(sds, 10), size = 3)
#' @importFrom slingshot slingPseudotime slingCurves reducedDim slingClusterLabels
#' @import ggplot2
#' @export

gg_plot <- function(sds, col = NULL, title = NULL, lineSize = 1, ...) {
  rd <- reducedDim(sds)
  
  if (is.null(col)) {
    cl <- slingClusterLabels(sds)
    if ("matrix" %in% is(cl)) {
      cl <- apply(cl, 1, which.max)
      cl <- as.character(cl)
    }
  } else {
    cl <- col
  }
  
  # Getting the main plot
  df <- data.frame(dim1 = rd[, 1], dim2 = rd[, 2], col = cl)
  p <- ggplot(df, aes(x = dim1, y = dim2, col = col)) +
    geom_point(...) +
    theme_classic() +
    labs(title = title, col = "")
  
  # Adding the curves
  for (i in seq_along(slingCurves(sds))) {
    curve_i <- slingCurves(sds)[[i]]
    curve_i <- curve_i$s[curve_i$ord, ]
    colnames(curve_i) <- c("dim1", "dim2")
    p <- p + geom_path(data = as.data.frame(curve_i), col = "black", size = 1)
  }
  return(p)
}