#' Plot variable features
#'
#' @description
#' This function creates a dotplot of Average expression vs. Standardized Variance for each
#' feature, highlighting the most variable ones of a seurat object.
#'
#' @param seurat_object Seurat object to use.
#' @param top_n Number of top variable feature to highlight. Default to 20.
#' @param cols 2-value vectors that specify non-variable/variable status colors. Default to black for non-variable
#' and blue for variable.
#'
#' @returns A ggplot object.
#'
#' @concept DimRed
#'
#' @examples
#' data("pbmc_small")
#' plot_variable_features(pbmc_small)
#'
#' # Specify the number of variable features to highlight
#' plot_variable_features(pbmc_small, top_n = 10)
#'
#' @export
plot_variable_features <- function(seurat_object,
                                   top_n = 20,
                                   cols = c("black", "blue")) {

  # Extract top_n variable genes
  top_variable <- utils::head(Seurat::VariableFeatures(seurat_object), top_n)

  # Create the plot
  var_plot <- suppressMessages(Seurat::LabelPoints(plot = Seurat::VariableFeaturePlot(seurat_object, cols = cols),
                                           points = top_variable, repel = T,
                                           xnudge = 0, ynudge = 0) +
                                 scale_x_continuous(labels = scales::comma, trans = "log10"))

  return(var_plot)
}



#' Plot splitted dimplots
#'
#' @description
#' This function creates a collection of dimplots: one generic grouped by the variable
#' set, and one for each unique value of the variable set, arranging them into a single plot.
#'
#' @details
#' Each subplot (one for each category of the group_by variable) will have the dots of the corresponding category colored
#' as specified in cols_vector (if specified) and all the other dots colored in gray, in order to see how the samples of that
#' category spread throughout the DimPlot compared to all the others.
#'
#' @param seurat_object Seurat object to use.
#' @param reduction Name of the reduction to use for plotting. Default to "pca".
#' @param group_by Name of the variable to use for grouping cells into dimplots. Default to "orig.ident".
#' @param dims 2-values vector with the dimension to use for plotting. Default to 1:2.
#' @param cols_vector Optional named vector that specifies the colors to use for each value of the group_by variable.
#' Default to NULL, using viridis as default palette.
#' @param heights 2-values vector to pass to ggarrange(height). It sets the height proportion between the main DimPlot and the
#' panel of the splitted dimplots. Default to c(1, 2).
#' @param axis_prefix Prefix to add to dims to label the axes. Default to "PCA".
#'
#' @returns A ggarrange object.
#'
#' @concept DimRed
#'
#' @examples
#' data("pbmc_small")
#' plot_splitted_dimplots(pbmc_small)
#'
#' # Specify the reduction, the axis prefix and the group_by variable
#' plot_splitted_dimplots(pbmc_small,
#'                        reduction = "tsne",
#'                        axis_prefix = "tsne",
#'                        group_by = "RNA_snn_res.1")
#'
#' @export
plot_splitted_dimplots <- function(seurat_object,
                                   reduction = "pca",
                                   group_by = "orig.ident",
                                   dims = 1:2,
                                   cols_vector = NULL,
                                   heights = c(1, 2),
                                   axis_prefix = "PCA") {

  # Get unique values of group_by column
  group_by_levels <- sort(unique(seurat_object@meta.data[, group_by]))

  # Create cols_vector if NA
  if (is.null(cols_vector)) {
    cols_vector <- viridisLite::viridis(length(group_by_levels))
    names(cols_vector) <- group_by_levels
  }

  # Create main plot
  main_plot <- DimPlot(seurat_object, reduction = reduction, group.by = group_by, pt.size = 0.8, dims = dims, shuffle = T, seed = 111) +
    scale_color_manual(values = cols_vector[names(cols_vector) %in% group_by_levels], name = group_by) +
    labs(title = group_by, x = paste(axis_prefix, "1"), y = paste(axis_prefix, "2")) +
    theme(legend.position = "bottom")

  # Create subplots
  ## Create empty list
  sub_list <- list()

  ## Create a loop
  for (group_level in group_by_levels) {

    ### Create vector of colors for this plot making all "gray" except the one active
    temp_colors <- cols_vector
    temp_colors[names(temp_colors) != group_level] <- "#D6D6D6"

    ### Create the plot
    sub_list[[group_level]] <- DimPlot(seurat_object, reduction = reduction, group.by = group_by, shuffle = T, seed = 222) +
      scale_color_manual(values = temp_colors[names(temp_colors) %in% group_by_levels]) +
      labs(title = group_level, x = paste(axis_prefix, "1"), y = paste(axis_prefix, "2")) +
      NoLegend()
  }

  # Arrange the list
  ## Decide ncol and nrow
  ncol <- ifelse(length(group_by_levels) <= 4,
                 yes = length(group_by_levels),
                 no = ceiling(sqrt(length(group_by_levels))))
  ncol <- ifelse(test = ncol <= 4, yes = ncol, no = 4) # max 4

  nrow <- ceiling(length(group_by_levels) / ncol)

  sub_arrange <- ggarrange(plotlist = sub_list, ncol = ncol, nrow = nrow, align = "hv")

  # Arrange everything
  all_arrange <- ggarrange(main_plot, sub_arrange, ncol = 1, nrow = 2, heights = heights)

  return(all_arrange)
}



#' Calculate reduction optimal dimensions
#'
#' @description
#' This function calculates the optimal dimensions to use for UMAP, tSNE, clustering etc starting from the values of a
#' pre-calculated dimensional reduction (e.g. pca).
#'
#' @details
#' This function will first calculate stdev associated with each embedding of a reduction, if not already calculated.
#' It then retrieve the % of variation associated with each PC, as well as the comulative sum; and finally it prints
#' both the dimension that exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5 (pc90) and the
#' first dimension that exhibits a difference with the next < 0.1 (pc01). Those two values are then stored in the \@misc slot of the
#' corresponding reduction.
#'
#'
#' @param seurat_object Seurat object to use.
#' @param reduction Name of the reduction to use. Default to "pca".
#'
#' @returns A Seurat object.
#'
#' @concept DimRed
#'
#' @examples
#' data("pbmc_small")
#'
#' # Specify the reduction, the axis prefix and the group_by variable
#' pbmc_small <- calc_dims(pbmc_small)
#'
#' @export
calc_dims <- function(seurat_object,
                      reduction = "pca") {

  # Calculate stdev if absent
  if (is_empty(seurat_object[[reduction]]@stdev)) {
    seurat_object[[reduction]]@stdev <- as.numeric(apply(seurat_object[[reduction]]@cell.embeddings, 2, stats::sd))
  }

  # Determine percent of variation associated with each PC
  pct <- seurat_object[[reduction]]@stdev / sum(seurat_object[[reduction]]@stdev) * 100

  # Calculate cumulative percents for each PC
  cumu <- cumsum(pct)

  # Determine which PC exhibits cumulative percent greater than 90% and % variation
  # associated with the PC as less than 5
  co1 <- which(cumu > 90 & pct < 5)[1]

  # Determine which PC is the first that exhibits a difference with the next < 0.1
  co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1),
              decreasing = T)[1] + 1

  # Store in misc
  seurat_object[[reduction]]@misc$pc90 <- co1
  seurat_object[[reduction]]@misc$pc01 <- co2

  # 7 Print
  cat(paste("The dimension that exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5 is", co1), "\n")
  cat(paste("The first dimension that exhibits a difference with the next < 0.1 is", co2), "\n")
  cat(paste0("This values are stored in seurat_object[[", reduction, "]]@misc."), "\n")

  return(seurat_object)

}
