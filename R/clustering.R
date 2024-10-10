#' Plot reduction at different resolutions
#'
#' @description
#' This function is used to plot different clustering resolutions on a reduction.
#'
#' @param seurat_object Seurat object to use.
#' @param reduction Name of the reduction to use. Default to "pca".
#' @param column_pattern Pattern used to select meta.data column to plot. Default to "RNA_snn_res".
#' @param label_size Value to pass to DimPlot(label_size). Default to 6.
#' @param axis_prefix Prefix to add to "1" and "2" as axis dimension names. Default to "PC_".
#' @param return_list Whether to return a list of ggplot objects (TRUE) or a ggarrange object with all the plots (FALSE). Default to FALSE.
#'
#' @returns A list of ggplot objects or a ggarrange object.
#'
#' @concept Clustering
#'
#' @examples
#' data("pbmc_small")
#'
#' # With default parameters
#' plot_reduction_res(pbmc_small)
#'
#' # Selecting different reduction
#' plot_reduction_res(pbmc_small, reduction = "tsne")
#'
#' @export

plot_reduction_res <- function(seurat_object,
                               reduction = "pca",
                               column_pattern = "RNA_snn_res",
                               label_size = 6,
                               axis_prefix = "PC_",
                               return_list = FALSE) {

  # Get the columns
  column_resolutions <- grep(pattern = column_pattern, colnames(seurat_object@meta.data), value = T)

  # Create empty list
  umaps <- list()

  # Loop to create plots
  for (res in column_resolutions) {
    umaps[[res]] <- DimPlot(object = seurat_object,
                            reduction = reduction,
                            group.by = res,
                            label = T,
                            label.size = label_size,
                            repel = T,
                            shuffle = T,
                            seed = 111) +
      NoLegend() +
      coord_fixed(ratio = 1) +
      labs(x = paste0(axis_prefix, "1") , y = paste0(axis_prefix, "2"))
  }

  # Return list if to return
  if (return_list) {return(umaps)}

  # Decide ncol and nrow
  ncol <- ifelse(length(column_resolutions) <= 3,
                 yes = length(column_resolutions),
                 no = ceiling(sqrt(length(column_resolutions))))
  ncol <- ifelse(test = ncol <= 3, yes = ncol, no = 3) # max 3

  nrow <- ceiling(length(column_resolutions) / ncol)

  # Arrange plots
  sub_arrange <- ggarrange(plotlist = umaps, ncol = ncol, nrow = nrow, align = "hv")


  return(sub_arrange)
}



#' Plot splitted highlighted clusters
#'
#' @description
#' This function generate plots to highlight the position of each cluster (or meta.data categorical column) in the reduction space.
#' It create a main plot (a classic DimPlot), along with a subplot for each cluster (or category), highlighting the position of
#' the samples belonging to it in the reduction space.
#'
#' @param seurat_object Seurat object to use.
#' @param reduction Name of the reduction to use. Default to "pca".
#' @param group_by Name of the column to use to create groups of samples. Default to "orig.ident".
#' @param dims 2-values vector with the dimension to use for plotting. Default to 1:2.
#' @param cols_vector  Optional named vector that specifies the colors to use for each value of the group_by variable.
#' Default to NULL, using turbo as default palette.
#' @param axis_prefix Prefix to add to dimensions as axis dimension names. Default to "PC_".
#' @param heights 2-values vector to pass to ggarrange(height). It sets the height proportion between the main DimPlot and the
#' panel of the splitted dimplots. Default to c(1, 2).
#' @param return_list Whether to return a list of ggplot objects (TRUE) or a ggarrange object with all the plots (FALSE). Default to FALSE.
#'
#' @returns A list of ggplot objects or a ggarrange object.
#'
#' @concept Clustering
#'
#' @examples
#' data("pbmc_small")
#'
#' # With default parameters
#' plot_splitted_clusters(pbmc_small)
#'
#' # Selecting different reduction and group_by column
#' plot_splitted_clusters(pbmc_small, reduction = "tsne",
#'                        group_by = "letter.idents",
#'                        axis_prefix = "TSNE_")
#'
#' @export
plot_splitted_clusters <- function(seurat_object,
                                   reduction = "pca",
                                   group_by = "orig.ident",
                                   dims = 1:2,
                                   cols_vector = NULL,
                                   axis_prefix = "PC_",
                                   heights = c(1, 2),
                                   return_list = FALSE) {

  # Get unique values of group_by column
  group_by_levels <- sort(unique(seurat_object@meta.data[, group_by]))

  # Create cols_vector if NULL
  if (is.null(cols_vector)) {
    cols_vector <- turbo(length(group_by_levels), begin = 0.1)
    names(cols_vector) <- group_by_levels
  }

  # Create main plot
  main_plot <- suppressMessages(DimPlot(seurat_object, reduction = reduction,
                                        group.by = group_by,
                                        pt.size = 0.8,
                                        dims = dims,
                                        label = T,
                                        label.size = 6,
                                        label.box = T,
                                        shuffle = T,
                                        seed = 111,
                                        repel = T) +
                                  scale_color_manual(values = cols_vector[names(cols_vector) %in% group_by_levels], name = group_by) +
                                  scale_fill_manual(values = rep("white", length(group_by_levels))) +
                                  labs(title = group_by, x = paste0(axis_prefix, dims[1]) , y = paste0(axis_prefix, dims[2])) +
                                  coord_fixed(ratio = 1) +
                                  NoLegend())

  # Create subplots
  ## Create empty list
  sub_list <- list()

  ## Loop through group_levels
  for (group_level in group_by_levels) {

    ### Create vector of colors for this plot making all "gray" except the one active
    temp_colors <- cols_vector
    temp_colors[names(temp_colors) != group_level] <- "#D6D6D6"

    ### Create the plot
    sub_list[[group_level]] <- DimPlot(seurat_object, reduction = reduction, group.by = group_by,
                                       shuffle = T,
                                       seed = 111) +
      scale_color_manual(values = temp_colors[names(temp_colors) %in% group_by_levels]) +
      labs(title = group_level, x = paste0(axis_prefix, dims[1]) , y = paste0(axis_prefix, dims[2])) +
      coord_fixed(ratio = 1) +
      NoLegend()
  }

  # Return list if to return
  if (return_list) {return(list("main" = main_plot, "subplots" = sub_list))}

  ## Arrange the list
  #### Decide ncol and nrow
  ncol <- ifelse(length(group_by_levels) <= 4,
                 yes = length(group_by_levels),
                 no = ceiling(sqrt(length(group_by_levels))))
  ncol <- ifelse(test = ncol <= 4, yes = ncol, no = 4) # max 4

  nrow <- ceiling(length(group_by_levels) / ncol)

  sub_arrange <- ggarrange(plotlist = sub_list, ncol = ncol, nrow = nrow, align = "hv")

  #### Arrange everything
  all_arrange <- ggarrange(main_plot, sub_arrange, ncol = 1, nrow = 2, heights = heights)

  return(all_arrange)
}



#' Plot clusters metadata
#'
#' @description
#' This function creates splitted dimplots based on a meta.data categorical column, a table of the number and
#' percentages of cells per cluster (or category of a categorical meta.data columns) for each category of the meta.data column.
#'
#' @param seurat_object Seurat object to use.
#' @param reduction Name of the reduction to use. Default to "pca".
#' @param group_by Name of the column to use to create groups of samples. Default to "seurat_clusters".
#' @param split_by Name of the column to use to split the plot for. Default to "orig.ident".
#' @param ncol Number of column to pass to ggarrange(ncol) when arranging the different dimplots.
#' @param dims 2-values vector with the dimension to use for plotting. Default to 1:2.
#' @param cols_vector  Optional named vector that specifies the colors to use for each value of the group_by variable.
#' Default to NULL, using turbo as default palette.
#' @param axis_prefix Prefix to add to dimensions as axis dimension names. Default to "PC_".
#' @param transpose Whether to transpose table to better fit the image. Default to FALSE.
#' @param heights 2-values vector to pass to ggarrange(height). It sets the height proportion between the main DimPlot and the
#' panel of the splitted dimplots. Default to c(2, 1, 1).
#' @param return_list Whether to return a list of ggplot objects (TRUE) or a ggarrange object with all the plots (FALSE). Default to FALSE.
#'
#' @returns A list of ggplot objects or a ggarrange object.
#'
#' @concept Clustering
#'
#' @examples
#' data("pbmc_small")
#'
#'
#' # Selecting reduction and group_by column
#' plot_clusters_metadata(pbmc_small, reduction = "tsne",
#'                        group_by = "letter.idents", split_by = "groups",
#'                        axis_prefix = "TSNE_")
#'
#' @export
plot_clusters_metadata <- function(seurat_object,
                                   reduction = "pca",
                                   group_by = "seurat_clusters",
                                   split_by = "orig.ident",
                                   ncol = 3,
                                   dims = 1:2,
                                   cols_vector = NULL,
                                   axis_prefix = "PC_",
                                   transpose = FALSE,
                                   heights = c(2, 1, 1),
                                   return_list = FALSE) {

  # To avoid no visible binding for global variable ‘.’
  . <- NULL

  # Get unique values of group_by column
  group_by_levels <- sort(unique(seurat_object@meta.data[, group_by]))

  # Create cols_vector if NA
  if (is.null(cols_vector)) {
    cols_vector <- turbo(length(group_by_levels), begin = 0.1)
    names(cols_vector) <- group_by_levels
  }

  # Create main plot
  main_plot <- suppressMessages(DimPlot(seurat_object,
                                        reduction = reduction,
                                        group.by = group_by,
                                        split.by = split_by,
                                        dims = dims,
                                        label = T,
                                        label.size = 3,
                                        ncol = ncol,
                                        shuffle = T,
                                        seed = 111) +
                                  scale_color_manual(values = cols_vector[names(cols_vector) %in% group_by_levels], name = group_by) +
                                  labs(title = split_by, x = paste0(axis_prefix, dims[1]) , y = paste0(axis_prefix, dims[2])) +
                                  coord_fixed(ratio = 1) +
                                  NoLegend())

  # Calculate cells for each clusters/split_by
  n_cells <- SeuratObject::FetchData(seurat_object,
                       vars = c(group_by, split_by)) %>%
    dplyr::count(!!rlang::sym(group_by), !!rlang::sym(split_by)) %>%
    dplyr::group_by(!!rlang::sym(split_by)) %>%
    tibble::as_tibble() %>%
    tidyr::spread(key = !!rlang::sym(group_by), value = .data$n, fill = "-", data = .) %>%
    as.data.frame()

  # Calculate percentage of cells for each clusters/split_by
  n_table <- table(seurat_object@meta.data[,split_by], seurat_object@meta.data[,group_by])
  p_table <- round(n_table / rowSums(n_table) * 100, 1)
  p_df <- as.data.frame(p_table) %>%
    tidyr::spread(key = .data$Var2, value = .data$Freq, data = .)

  colnames(p_df)[1] <- split_by

  p_df[p_df == 0] <- "-"

  if (transpose) {
    n_cells <- t(n_cells)
    p_df <- t(p_df)
  }

  # Create tables to plot
  n_table_plot <- ggtexttable(n_cells,
                              rows = NULL,
                              theme = ttheme(base_style = "light",
                                             base_size = 12)) %>%
    tab_add_title(text = "Number of cells", face = "bold")

  p_table_plot <- ggtexttable(p_df,
                              rows = NULL,
                              theme = ttheme(base_style = "light",
                                             base_size = 12))  %>%
    tab_add_title(text = "Fraction of cells", face = "bold")

  # Return list if to return
  if (return_list) {return(list("main" = main_plot, "n_table" = n_table_plot, "p_table" = p_table_plot))}

  return(ggarrange(main_plot, n_table_plot, p_table_plot, ncol = 1, nrow = 3, heights = heights))
}
