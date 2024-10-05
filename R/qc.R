#' Calculate QC metrics
#'
#' @description
#' This function calculates QC metrics (log10GenesPerUMI, mitochondrial gene ratio, ribosomal gene ratio and hemoglobin gene ratio)
#' for the provided seurat object.
#'
#' @param seurat_object Seurat object of which calculate QC metrics.
#' @param assay Assay to use. By default it takes the current active assay.
#' @param mito Pattern to use to match for mitochondrial genes. Default to "^MT-".
#' @param ribo Pattern to use to match for ribosomal genes. Default to "^RP\[SL\]".
#' @param hemo Pattern to use to match for hemoglobin genes. Default to "^HB\[^(P)\]".
#'
#' @returns A Seurat object with 4 new columns in meta.data:
#' * log10GenesPerUMI
#' * mitoRatio
#' * riboRatio
#' * hemoRatio
#'
#' @concept QC
#'
#' @examples
#' data("pbmc_small")
#' pbmc_small <- calculate_qc_metrics(pbmc_small)
#'
#' @export
calculate_qc_metrics <- function(seurat_object,
                                 assay = NULL,
                                 mito = "^MT-",
                                 ribo = "^RP[SL]",
                                 hemo = "^HB[^(P)]") {

  # Check assay
  assay <- assay %||% DefaultAssay(object = seurat_object)

  # Calculate the log10 number of genes per UMI
  seurat_object$log10GenesPerUMI <- log10(seurat_object@meta.data[, paste0("nFeature_", assay)]) /
                                      log10(seurat_object@meta.data[, paste0("nCount_", assay)])

  # Calculate fraction of mitochondrial genes
  seurat_object$mitoRatio <- PercentageFeatureSet(seurat_object, pattern = mito, assay = assay) / 100

  # Calculate fraction of ribosomal genes
  seurat_object$riboRatio <- PercentageFeatureSet(seurat_object, pattern = ribo, assay = assay) / 100

  # Calculate fraction of hemoglobin genes
  seurat_object$hemoRatio <- PercentageFeatureSet(seurat_object, pattern = hemo, assay = assay) / 100

  # Store command in object
  seurat_object <- SeuratObject::LogSeuratCommand(seurat_object)

  return(seurat_object)

}



#' Create a density plot for a QC metric
#'
#' @description
#' This function creates a density plot for the given qc metric (or any continuous variable
#' of the Seurat object meta.data).
#'
#' @param seurat_object Seurat object to use to create the plot.
#' @param metric Name of the metric (meta.data column) to plot. Default to "nCount_RNA".
#' @param title Title to add to the plot. Default to "nCount distribution".
#' @param fill Name of the column to use to pass to aes(fill) in geom_density.
#' It is used to split the data into groups and draw a density plot for each value of it. Default to "orig.ident".
#' @param fill_vector Optional named vector used to assign fill colors to the different group of fill variable.
#' Default to NULL.
#' @param alpha Value to pass to alpha in geom_density. Default to 0.2.
#' @param threshold Optional value of the threshold line to plot as a vertical line. Default to NULL, not plotted.
#' @param x_trans Transformation to pass to scale_x_continuous(trans). Default to "identity".
#'
#' @returns A ggplot list object with the plot.
#'
#' @concept QC
#'
#' @examples
#' data("pbmc_small")
#' pbmc_small <- calculate_qc_metrics(pbmc_small)
#'
#' # Basic plot
#' qc_density_plot(pbmc_small)
#'
#' # Specify fill colors
#' qc_density_plot(pbmc_small, fill = "groups", fill_vector = c("g1" = "blue", "g2" = "red"))
#'
#' @export
qc_density_plot <- function(seurat_object,
                            metric = "nCount_RNA",
                            title = "nCount distribution",
                            fill = "orig.ident",
                            fill_vector = NULL,
                            alpha = 0.2,
                            threshold = NULL,
                            x_trans = "identity") {

  # Create the plot
  density_plot <- ggplot(seurat_object@meta.data) +
    geom_density(aes_string(x = metric, fill = fill),
                 alpha = alpha) +
    labs(title = title, y = "") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    scale_x_continuous(labels = scales::comma,
                       expand = expansion(mult = c(0,0)),
                       trans = x_trans) +
    theme_classic() +
    theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
          axis.title = element_text(face = "bold", size = 12),
          axis.text = element_text(size = 10),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          legend.position = "right")

  # Add threshold line if provided
  if (!is.null(threshold)) {

    # Calculate the max value of density
    max_density <- ggplot_build(density_plot)$layout$panel_scales_y[[1]]$range$range[2]

    # Add to plot
    density_plot <- density_plot +
      geom_vline(xintercept = threshold, color = "red", size = 1.2) +
      geom_text(aes(x = threshold, max_density, label = "Threshold"),
                color = "red",
                hjust = -0.1, size = 4)
  }

  # Add fill color vector if provided
  if (!is.null(fill_vector)) {
    density_plot <- density_plot +
      scale_fill_manual(values = fill_vector)
  }


  return(density_plot)
}



#' Create all density plots for a QC metrics
#'
#' @description
#' This function is a wrapper around qc_density_plot. It creates density plots for different QC metrics
#' (nCount, nFeature, log10GenesPerUMI, mitochondrial gene ratio, ribosomal gene ratio and hemoglobin gene ratio) and wrap them up
#' into a single plot.
#'
#' @details For nCount and nFeature density plots, the default x_trans is set to log10.
#'
#' @param seurat_object Seurat object to use to create the plot.
#' @param assay Assay to use to select nCount and nFeature meta.data. By default it takes the current active assay.
#' @param fill Name of the column to use to pass to aes(fill) in geom_density.
#' It is used to split the data into groups and draw a density plot for each value of it. Default to "orig.ident".
#' @param fill_vector Optional named vector used to assign fill colors to the different group of fill variable.
#' Default to NULL.
#' @param alpha Value to pass to alpha in geom_density. Default to 0.2.
#' @param count_thresh Value of the count threshold line to pass to qc_density_plot(threshold). Default to 1000.
#' @param feature_thresh Value of the feature threshold line to pass to qc_density_plot(threshold). Default to 1000.
#' @param log10GenesPerUMI_thresh Value of the log10GenesPerUMI threshold line to pass to qc_density_plot(threshold). Default to 0.8
#' @param mito_thresh Value of the mitoRatio threshold line to pass to qc_density_plot(threshold). Default to 0.15.
#' @param ribo_thresh Value of the riboRatio threshold line to pass to qc_density_plot(threshold). Default to 0.05.
#' @param hemo_thresh Value of the hemoRatio threshold line to pass to qc_density_plot(threshold). Default to 0.NULL.
#'
#' @returns A ggarrange object with the plots.
#'
#' @concept QC
#'
#' @examples
#' data("pbmc_small")
#' pbmc_small <- calculate_qc_metrics(pbmc_small)
#'
#' # Basic plot
#' plot_qc_metrics(pbmc_small)
#'
#' # Specify other thresholds
#' plot_qc_metrics(pbmc_small, count_thresh = 750,
#'                 feature_thresh = 900,
#'                 mito_thresh = 0.2)
#'
#' @export
plot_qc_metrics <- function(seurat_object,
                            assay = NULL,
                            fill = "orig.ident",
                            fill_vector = NULL,
                            alpha = 0.2,
                            count_thresh = 1000,
                            feature_thresh = 1000,
                            log10GenesPerUMI_thresh = 0.8,
                            mito_thresh = 0.15,
                            ribo_thresh = 0.05,
                            hemo_thresh = NULL) {

  # Check if calculate_qc_metrics have been run
  if (!any((grepl("calculate_qc_metrics", names(seurat_object@commands))))) {
    stop("Please run calculate_qc_metrics on this object before using this function.")
  }

  # Check assay
  assay <- assay %||% DefaultAssay(object = seurat_object)

  # Create plot for number of counts
  count_plot <- qc_density_plot(seurat_object = seurat_object,
                                metric = paste0("nCount_", assay),
                                title = "nCount/cell distribution",
                                fill = fill,
                                fill_vector = fill_vector,
                                alpha = alpha,
                                threshold = count_thresh,
                                x_trans = "log10")

  # Create plot for number of features
  feature_plot <- qc_density_plot(seurat_object = seurat_object,
                                  metric = paste0("nFeature_", assay),
                                  title = "nFeature/cell distribution",
                                  fill = fill,
                                  fill_vector = fill_vector,
                                  alpha = alpha,
                                  threshold = feature_thresh,
                                  x_trans = "log10")

  # Create plot for complexity score
  complexity_plot <- qc_density_plot(seurat_object = seurat_object,
                                     metric = "log10GenesPerUMI",
                                     title = "Log10 Features/Counts distribution",
                                     fill = fill,
                                     fill_vector = fill_vector,
                                     alpha = alpha,
                                     threshold = log10GenesPerUMI_thresh)

  # Create plot for mitochondrial genes
  mito_plot <- qc_density_plot(seurat_object = seurat_object,
                               metric = "mitoRatio",
                               title = "Mitochondrial genes expression distribution",
                               fill = fill,
                               fill_vector = fill_vector,
                               alpha = alpha,
                               threshold = mito_thresh)

  # Create plot for ribosomal genes
  ribo_plot <- qc_density_plot(seurat_object = seurat_object,
                               metric = "riboRatio",
                               title = "Ribosomal genes expression distribution",
                               fill = fill,
                               fill_vector = fill_vector,
                               alpha = alpha,
                               threshold = ribo_thresh)

  # Create plot for hemoglobin genes
  hemo_plot <- qc_density_plot(seurat_object = seurat_object,
                               metric = "hemoRatio",
                               title = "Hemoglobin genes expression distribution",
                               fill = fill,
                               fill_vector = fill_vector,
                               alpha = alpha,
                               threshold = hemo_thresh)

  # Merge all plots
  all_plots <- ggarrange(count_plot, feature_plot, complexity_plot, mito_plot, ribo_plot, hemo_plot,
                         ncol = 3, nrow = 2,
                         align = "hv" ,
                         common.legend = T)

  return(all_plots)
}



#' Create a correlation plot between 2 different QC metrics.
#'
#' @description
#' This function creates a correlation plot between two given qc metric (or any continuous variable
#' of the Seurat object meta.data).
#'
#' @param seurat_object Seurat object to use to create the plot.
#' @param metric_x Name of the metric (meta.data column) to plot on x axis. Default to "nCount_RNA".
#' @param metric_y Name of the metric (meta.data column) to plot on y axis. Default to "nFeature_RNA".
#' @param title Title to add to the plot. Default to "".
#' @param shape Name of the column to use to pass to aes(shape) in geom_point.
#' It is used to change the shape of the dots based on values of a meta.data column. Default to "orig.ident".
#' @param color Name of the column to use to pass to aes(color) in geom_point.
#' It is used to change the color of the dots based on values of a meta.data column. Default to "mitoRatio".
#' @param x_threshold Optional value of the threshold line to plot as a vertical line. Default to NULL, not plotted.
#' @param y_threshold Optional value of the threshold line to plot as an horizontal line. Default to NULL, not plotted.
#' @param x_trans Transformation to pass to scale_x_continuous(trans). Default to "identity".
#' @param y_trans Transformation to pass to scale_y_continuous(trans). Default to "identity".
#'
#' @returns A ggplot list object with the plot.
#'
#' @concept QC
#'
#' @examples
#' data("pbmc_small")
#' pbmc_small <- calculate_qc_metrics(pbmc_small)
#'
#' # Basic plot
#' qc_correlation_plot(pbmc_small)
#'
#' # Specify x and y metrics
#' qc_correlation_plot(pbmc_small, metric_x = "mitoRatio",
#'                     metric_y = "riboRatio")
#'
#' @export
qc_correlation_plot <- function(seurat_object,
                                metric_x = "nCount_RNA",
                                metric_y = "nFeature_RNA",
                                title = "",
                                shape = "orig.ident",
                                color = "mitoRatio",
                                x_threshold = NULL,
                                y_threshold = NULL,
                                x_trans = "identity",
                                y_trans = "identity") {

  # Create the plot
  correlation_plot <- ggplot(seurat_object@meta.data,
                             aes_string(x = metric_x, y = metric_y)) +
    geom_point(aes_string(color = color, shape = shape)) +
    labs(title = title) +
    scale_y_continuous(labels = scales::comma,
                       expand = expansion(mult = c(0, 0.05)),
                       trans = y_trans) +
    scale_x_continuous(labels = scales::comma,
                       expand = expansion(mult = c(0,0.05)),
                       trans = x_trans) +
    scale_color_gradient(low = "gray90", high = "black") +
    theme_classic() +
    theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
          axis.title = element_text(face = "bold", size = 12),
          axis.text = element_text(size = 10),
          legend.position = "right")


  # Add x_threshold line if provided
  if (!is.null(x_threshold)) {

    # Calculate the max value of density
    max_y <- max(seurat_object@meta.data[, metric_y])

    # Add to plot
    correlation_plot <- correlation_plot +
      geom_vline(xintercept = x_threshold, color = "red", size = 1) +
      geom_text(aes(x = x_threshold, y = max_y, label = "Threshold"),
                color = "red",
                hjust = -0.1, size = 4)
  }

  # Add y_threshold line if provided
  if (!is.null(y_threshold)) {

    # Calculate the max value of density
    max_x <- max(seurat_object@meta.data[, metric_x])

    # Add to plot
    correlation_plot <- correlation_plot +
      geom_hline(yintercept = y_threshold, color = "red", size = 1) +
      geom_text(aes(y = y_threshold, x = max_x, label = "Threshold"),
                color = "red",
                vjust = 1.3,
                hjust = 1,
                size = 4)
  }

  return(correlation_plot)
}



#' Plot top expressed features.
#'
#' @description
#' This function creates a boxplot representing the top n features, based on median value of % of expression in cells.
#'
#' @param seurat_object Seurat object to use to create the plot.
#' @param assay Name of the assay to use for extrapolating the data. Default to "RNA".
#' @param slot Name of the slot to use. Default to "counts".
#' @param n Number of top features to count. Default to 10.
#'
#' @returns A ggplot list object with the plot.
#'
#' @concept QC
#'
#' @examples
#' data("pbmc_small")
#'
#' # Basic plot
#' plot_top_n_genes(pbmc_small)
#'
#' # Specify n
#' plot_top_n_genes(pbmc_small, n = 25)
#'
#' @export
plot_top_n_genes <- function(seurat_object,
                             assay = "RNA",
                             slot = "count",
                             n = 10) {

  # Get the most expressed genes
  mat <- GetAssayData(seurat_object, slot = slot, assay = assay)

  # Convert to sparse if not already
  if (!inherits(mat, "dgCMatrix")) {
    mat <- as(mat, "dgCMatrix")
  }

  col_sums <- Matrix::colSums(mat)
  mat <- sweep(mat, 2, col_sums, FUN = "/") * 100

  most_expressed <- order(apply(mat, 1, median), decreasing = T)[n:1]

  # Reshape for ggplot
  mat <- stack(as.data.frame(Matrix::t(mat[most_expressed, ])))

  # Create the plot
  most_expressed_box <- ggplot(mat) +
    geom_boxplot(aes_string(x = "ind", y = "values", fill = "ind"), show.legend = F, outlier.size = 0.2) +
    labs(x = "", y = "% counts/cell", title = paste("Top", n, "features")) +
    scale_y_continuous(expand = expansion(mult = c(0, 0))) +
    coord_flip() +
    theme_classic() +
    theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
          axis.title = element_text(face = "bold", size = 12),
          axis.text = element_text(size = 10),
          legend.position = "right")

  return(most_expressed_box)
}



