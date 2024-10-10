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
    theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 8),
          axis.title = element_text(face = "bold", size = 6),
          axis.text = element_text(size = 5),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          legend.position = "right")

  # Add threshold line if provided
  if (!is.null(threshold)) {

    # Calculate the max value of density
    max_density <- ggplot_build(density_plot)$layout$panel_scales_y[[1]]$range$range[2]

    # Add to plot
    density_plot <- density_plot +
      geom_vline(xintercept = threshold, color = "red", size = 0.8) +
      geom_text(aes(x = threshold, max_density, label = "Threshold"),
                color = "red",
                hjust = -0.1, size = 2)
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
#' @param return_list Whether to return a list of ggplot objects (TRUE) or a ggarrange object with all the plots (FALSE). Default to FALSE.
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
                            hemo_thresh = NULL,
                            return_list = FALSE) {

  # Check if calculate_qc_metrics have been run
  if (!any((grepl("calculate_qc_metrics", names(seurat_object@commands))))) {
    stop("Please run calculate_qc_metrics on this object before using this function.")
  }

  # Get assay
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

  # Return list if to return
  if (return_list) {return(list("count" = count_plot,
                                "feature" = feature_plot,
                                "lo10UmiPerCell" = complexity_plot,
                                "mito" = mito_plot,
                                "ribo" = ribo_plot,
                                "hemo" = hemo_plot))}

  # Merge all plots
  all_plots <- ggarrange(plotlist = list(count_plot, feature_plot, complexity_plot, mito_plot, ribo_plot, hemo_plot),
                         align = "hv" ,
                         common.legend = T)

  return(all_plots)
}



#' Create a correlation plot between 2 different QC metrics
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
                       expand = expansion(mult = c(0.05, 0.05)),
                       trans = y_trans) +
    scale_x_continuous(labels = scales::comma,
                       expand = expansion(mult = c(0.05, 0.05)),
                       trans = x_trans) +
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



#' Calculate QC metrics outliers with MADs
#'
#' @description
#' This function calculates, for each qc metric previously calculated with [calculate_qc_metrics()] using the MAD approach
#' (check [scater::isOutlier()] for more info).
#' By default, nCount and nFeature for the current assay, as well as mitoRatio and riboRatio are examinated. To see default values,
#' run [view_default_qc_mads()].
#' For each metric, it adds a new column in meta.data called metric_outlier which mark whether a cell is an outlier (TRUE) for that metric or not (FALSE);
#' it also creates a dataframe in seurat_object@misc$qc_thresholds containing the calculated threshold values.
#'
#' @section Add custom metrics:
#' To add custom metrics, provide a list as argument to `extra` in this format:
#' list("name_of_meta.data_column" = list("nmads" = c(lower_value, upper_value), "log" = boolean)). Both nmads and log are REQUIRED and cannot be omitted.
#' lower_value or upper_value can be set to `NA` if that threshold should not be calculated.
#'
#' @section Override/skip default values:
#' To override default values, provide the corresponding list in extra argument as if it is an extra column to calculate. It will override default parameters for that metric.
#' To skip a default metric, set it to NULL as metric element list (e.g `extra = list("mitoRatio" = NULL)`)
#'
#'
#' @param seurat_object Seurat object to use.
#' @param assay Assay to use to select nCount and nFeature meta.data. By default it takes the current active assay.
#' @param batch Name of the column to pass to `isOutlier(batch)`. Default to "orig.ident".
#' @param extra Optional list containing custom metrics to calculate or default values to override.
#'
#' @returns A Seurat object.
#'
#' @concept QC
#'
#' @examples
#' data("pbmc_small")
#' pbmc_small <- calculate_qc_metrics(pbmc_small)
#'
#' # Default with no customization
#' pbmc_small <- calculate_qc_mad_outliers(pbmc_small)
#'
#' # Setting batch and omit riboRatio calculation
#' pbmc_small <- calculate_qc_mad_outliers(pbmc_small,
#'                                         batch = "groups",
#'                                         extra = list("riboRatio" = NULL))
#'
#'
#' @export
calculate_qc_mad_outliers <- function(seurat_object,
                                      assay = NULL,
                                      batch = "orig.ident",
                                      extra = list()) {

  # Get assay
  assay <- assay %||% DefaultAssay(object = seurat_object)

  # Create empty dataframe to store threshold values
  thresholds_df <- data.frame("metric" = character(),
                              "type" = character(),
                              "batch" = character(),
                              "value" = numeric()
  )

  # Create parameter list
  parameter_list <- create_parameter_list(assay = assay, extra = extra)

  # Loop through parameters
  for (parameter_name in unique(names(parameter_list))) {

    # Skip if null
    if (!is_null(parameter_list[[parameter_name]])) {

      # Create temporary boolean vector that stores F, so that when | with isOutlier results
      # will be T only for real outliers
      tmp_out <- rep(F, nrow(seurat_object@meta.data))

      # Skip if higher mads is NA
      if (!is.na(parameter_list[[parameter_name]]$nmads[2])) {
        tmp_out_higher <- isOutlier(seurat_object@meta.data[, parameter_name],
                                    nmads = parameter_list[[parameter_name]]$nmads[2],
                                    type = "higher",
                                    log = parameter_list[[parameter_name]]$log,
                                    batch = seurat_object@meta.data[, batch])

        tmp_out <- tmp_out | tmp_out_higher

        # Add to threshold df
        tmp_thresholds <- attr(tmp_out_higher, "thresholds")
        thresholds_df <- rbind(thresholds_df,
                               data.frame("metric" = rep(parameter_name, dim(tmp_thresholds)[2]),
                                          "type" = c("higher"),
                                          "batch" = colnames(tmp_thresholds),
                                          "value" = tmp_thresholds[2, ]
                                          )
                               )

      }

      # Skip if lower mads is NA
      if (!is.na(parameter_list[[parameter_name]]$nmads[1])) {
        tmp_out_lower <- isOutlier(seurat_object@meta.data[, parameter_name],
                                   nmads = parameter_list[[parameter_name]]$nmads[1],
                                   type = "lower",
                                   log = parameter_list[[parameter_name]]$log,
                                   batch = seurat_object@meta.data[, batch])

        tmp_out <- tmp_out | tmp_out_lower

        # Add to threshold df
        tmp_thresholds <- attr(tmp_out_lower, "thresholds")
        thresholds_df <- rbind(thresholds_df,
                               data.frame("metric" = rep(parameter_name, dim(tmp_thresholds)[2]),
                                          "type" = c("lower"),
                                          "batch" = colnames(tmp_thresholds),
                                          "value" = tmp_thresholds[1, ]
                                          )
                               )
      }

      # Create column in seurat object
      seurat_object@meta.data[, paste0(parameter_name, "_outlier")] <- tmp_out

    }


  }

  # Insert thresholds df in seurat object
  seurat_object@misc$qc_thresholds <- thresholds_df

  return(seurat_object)
}



#' Create parameter list for MAD outlier
#'
#' @description
#' This function creates a list with the parameter to use for MAD outlier QC calculation.
#'
#' @param assay Assay to use to select nCount and nFeature meta.data. By default it takes the current active assay.
#' @param extra Optional list containing custom metrics to calculate or default values to override.
#'
#' @returns A list with elemets like this: "name_of_meta.data_column" = list("nmads" = c(lower_value, upper_value), "log" = boolean).
#'
#' @concept QC
create_parameter_list <- function(assay = "RNA",
                                  extra = list()) {

  # Set default parameters
  defaults <- list(
     "nCount" = list(
      "nmads" = c(2, 3),
      "log" = T
    ),
     "nFeature" = list(
      "nmads" = c(2, 3),
      "log" = T
    ),
    mitoRatio = list(
      "nmads" = c(NA, 2),
      "log" = T
    ),
    riboRatio = list(
      "nmads" = c(2, NA),
      "log" = T
    )
  )

  # Correct names of first 2 elements of defaults
  names(defaults)[1:2] <- c(paste0("nCount_", assay), paste0("nFeature_", assay))

  return(c(extra,
           defaults))
}



#' View default QC MAD parameters
#'
#' @description
#' This function is used to see the default values for [calculate_qc_mad_outliers()] metrics.
#'
#' @returns A tibble with columns:
#' * metric: name of the metric
#' * lower.nmads: value of the lower nmads used to calculate outliers
#' * higher.nmads: value of the higher nmads used to calculate outliers
#' * log: whether outlier detection is done on log-transformed data
#'
#' @concept QC
#'
#' @examples
#' view_default_qc_mads()
#'
#' @export
view_default_qc_mads <- function(){

  parameter_list <- create_parameter_list()

  print(as.data.frame(unlist(parameter_list)) %>%
    tibble::rownames_to_column("what") %>%
    tidyr::separate(col = "what", into = c("metric", "type"), sep = "\\.") %>%
    dplyr::mutate("type" = dplyr::case_when(
      type == "nmads1" ~ "lower.nmads",
      type == "nmads2" ~ "higher.nmads",
      TRUE ~ .data$type
    )) %>%
    tidyr::pivot_wider(names_from = .data$type, values_from = 3) %>%
    dplyr::mutate(log = as.logical(log))
  )

}



#' Plot qc metrics outliers
#'
#' @description
#' This function is a wrapper around [plot_violin_outliers()]. It creates a list of violin plots for the desired qc metrics (or meta.data column).
#'
#' @param seurat_object Seurat object to use.
#' @param split_by Name of the column used to group_by the violins (usually should correspond to the one used as `batch` in [calculate_qc_mad_outliers()]).
#' Default to "orig.ident".
#' @param fill_vector Optional named vector used to assign fill colors to the different group of fill variable.
#' Default to NULL.
#' @param metrics Named vector with metrics to plot as names and the desiretd x_trans as value.
#' Default to metrics = c("nCount_RNA" = "log10", "nFeature_RNA" = "log10", "mitoRatio" = "identity", "riboRatio" = "identity")
#' @param return_list Whether to return a list of ggplot objects (TRUE) or a ggarrange object with all the plots (FALSE). Default to FALSE.
#'
#' @returns A list of ggplot objects or a ggarrange object.
#'
#' @concept QC
#'
#' @examples
#' data("pbmc_small")
#' pbmc_small <- calculate_qc_metrics(pbmc_small)
#' pbmc_small <- calculate_qc_mad_outliers(pbmc_small)
#'
#' # Basic plot
#' plot_qc_metrics_outliers(pbmc_small)
#'
#' # Split violin based on a column
#' pbmc_small <- calculate_qc_mad_outliers(pbmc_small, batch = "groups")
#' plot_qc_metrics_outliers(pbmc_small, split_by = "groups")
#'
#' @export
plot_qc_metrics_outliers <- function(seurat_object,
                                     split_by = "orig.ident",
                                     fill_vector = NULL,
                                     metrics = c("nCount_RNA" = "log10",
                                                 "nFeature_RNA" = "log10",
                                                 "mitoRatio" = "identity",
                                                 "riboRatio" = "identity"),
                                     return_list = FALSE) {

  # Create empty list to store plots
  plot_list <- list()

  # Loop through metrics
  for (metric in names(metrics)) {
    plot_list[[metric]] <- plot_violin_outliers(seurat_object = seurat_object,
                                                split_by = split_by,
                                                metric = metric,
                                                title = paste0(metric, "/cell distribution"),
                                                fill_vector = fill_vector,
                                                x_trans = metrics[metric])
  }

  # Return list if to return
  if (return_list) {return(plot_list)}

  # Merge all plots
  all_plots <- ggarrange(plotlist = plot_list,
                         align = "hv" ,
                         common.legend = T)

  return(all_plots)
}


#' Plot violin outlier thresholds
#'
#' @description
#' This function creates is used to represent the distribution of the values of a desired metric (meta.data column)
#' and the corresponding thresholds (if available). It can help in the visual inspection of nmads to use when calculating
#' outliers with [calculate_qc_mad_outliers()].
#'
#' @param seurat_object Seurat object to use.
#' @param split_by Name of the column used to group_by the violins (usually should correspond to the one used as `batch` in [calculate_qc_mad_outliers()]).
#' Default to "orig.ident".
#' @param metric Name of the column to plot. Default to "nCount_RNA".
#' @param title Title to add to the plot. Default to "".
#' @param x_trans Transformation to pass to scale_x_continuous(trans). Default to "identity".
#' @param fill_vector Optional named vector used to assign fill colors to the different group of fill variable.
#' Default to NULL.
#'
#' @returns A ggplot list object with the plot.
#'
#' @concept QC
#'
#' @examples
#' data("pbmc_small")
#' pbmc_small <- calculate_qc_metrics(pbmc_small)
#' pbmc_small <- calculate_qc_mad_outliers(pbmc_small)
#'
#' # Basic plot
#' plot_violin_outliers(pbmc_small)
#'
#' # Changing metric
#' plot_violin_outliers(pbmc_small, metric = "riboRatio")
#'
#' @export
plot_violin_outliers <- function(seurat_object = seurat_object,
                                 split_by = "orig.ident",
                                 metric = "nCount_RNA",
                                 title = "",
                                 x_trans = "identity",
                                 fill_vector = NULL) {

  # Filter threshold df
  tmp_thresholds_df <- seurat_object@misc$qc_thresholds[seurat_object@misc$qc_thresholds$metric == metric,]
  tmp_thresholds_df$batch <- factor(tmp_thresholds_df$batch)

  # Create plot
  violin_plot <- ggplot() +
    geom_violin(data = seurat_object@meta.data,
                mapping = aes_string(y = split_by,
                                     x = metric,
                                     fill = split_by),
                draw_quantiles = 0.5) +
    geom_segment(data = tmp_thresholds_df,
                 aes(y = as.numeric(.data$batch) - .3, yend = as.numeric(.data$batch) + .3, x = .data$value, xend = .data$value),
                 linewidth = .8,
                 linetype = "11") +
    scale_y_discrete(expand = expansion(mult = c(0.01, 0))) +
    scale_x_continuous(trans = x_trans,
                       labels = function(x)scales::comma(x),
                       expand = expansion(mult = c(0, 0.15))) +
    theme_classic() +
    theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
          axis.title = element_text(face = "bold", size = 12),
          axis.text = element_text(size = 10),
          legend.position = "right")

  # Add fill colors if provided
  if (!is.null(fill_vector)) {
    violin_plot <- violin_plot +
      scale_fill_manual(values = fill_vector)
  }

  # 4. Return
  return(violin_plot)

}



#' Plot top expressed features
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
#' plot_top_n_features(pbmc_small)
#'
#' # Specify n
#' plot_top_n_features(pbmc_small, n = 25)
#'
#' @export
plot_top_n_features <- function(seurat_object,
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



