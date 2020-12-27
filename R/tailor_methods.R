###############################
# Functions exposed to the user
###############################

#' @import Rtsne
#' @import KernSmooth
#' @import cluster
#' @import mvtnorm
#' @import ggplot2
#' @importFrom gridExtra grid.arrange
#' @importFrom flowCore flowFrame exprs
#' @importFrom mclust Mclust mclustBIC
#' @import foreach
#' @import doParallel
#' @importFrom grDevices dev.print png
#' @importFrom graphics abline lines par plot
#' @importFrom stats density dist dnorm kmeans
#' @importFrom methods is as
#' @import parallel
#' @import iterators
#' @importFrom RColorBrewer brewer.pal
#' @importFrom rlang .data
NULL

#' @title tailor_learn
#' @description This function learns a tailor model from input
#' data. It computes a preliminary binning of the data, then
#' computes a mixture model using a weighted version of the
#' expectation-maximization (EM) algorithm,
#' and finally merges mixture components which are positive/negative
#' for the same markers, using adaptive thresholds.
#' @param data A flowSet, flowFrame or a matrix containing
#' events along the rows, markers along columns.
#' @param params A list of markers to use; must be subset of
#' colnames(data).
#' @param mixture_components The number of mixture
#' components to learn. Some of these
#' are eventually merged, so it's a good idea to
#' choose a number slightly larger than
#' the number of clusters you expect to get.
#' @param min_bin_size Bins with fewer events than this threshold
#' are considered outliers, and are
#' ignored during the weighted EM algorithm. These events can
#' still be assigned to clusters during
#' the prediction phase.
#' @param max_bin_size Bins with more events than this threshold
#' are split, to ensure that the
#' weighted EM algorithm closely approximates a run of vanilla
#' EM on the entire dataset.
#' @param min_cluster_fraction Mixture components whose size is
#' smaller than this fraction are dropped from further analysis,
#' and a warning is returned.
#' @param mixtures_1D Pre-computed 1D mixture models,
#'  to be used for binning.
#' These are computed from scratch if not provided.
#' @param parallel Boolean flag; if true, uses multithreading
#' to speed up computation.
#' @param verbose If > 0, outputs milestone information.
#' If >=1, also outputs information about
#' running utilities. If >1, debugging mode.
#' @return A tailor object containing:
#' \describe{
#'   \item{fit}{The tailor model, a named list containing
#'   the mixture proportions, means and variances
#'   of all mixture components.}
#'   \item{mixtures_1D}{}
#'   \item{cat_clusters}{A named list containing information
#'   about the categorical clusters found by
#'   the model: phenotype, cluster centers, and a mapping
#'   from mixture components to categorical clusters.}
#' }
#' @examples
#' # Load data and define analytical parameters
#' fileName <- system.file("extdata", "sampled_flowset_old.rda",
#'                         package = "Tailor")
#' load(fileName)
#' tailor_params <- flowCore::colnames(fs_old)[c(7:9, 11:22)]
#'
#' # Run with default settings
#' tailor_obj <- tailor_learn(data = fs_old,
#'                           params = tailor_params,
#'                           mixture_components = 50)
#'
#' # Alternatively, customize the 1D mixtures used for binning step
#' mixtures_1D <- get_1D_mixtures(fs_old, tailor_params)
#' to_customize <- list("CD127BV421" = 2)
#' mixtures_1D <- customize_1D_mixtures(fs_old, to_customize, mixtures_1D)
#'
#' tailor_obj <- tailor_learn(data = fs_old,
#'                           params = tailor_params,
#'                           mixture_components = 50,
#'                           mixtures_1D = mixtures_1D)
#' @export
tailor_learn <- function(data, params = NULL, mixtures_1D = NULL,
                        mixture_components = 100,
                        min_bin_size = NULL, max_bin_size = NULL,
                        min_cluster_fraction = 1e-3,
                        parallel = FALSE, verbose = 0.5)
{
  data <- as_matrix(data)
  if (is.null(params)) params <- colnames(data)

  if (is.null(min_bin_size)) {
    min_bin_size <- max(5,ceiling(nrow(data) / 1e5))
  }
  if(is.null(max_bin_size)) {
    max_bin_size <- max(50, ceiling(nrow(data) / 1e4))
  }
  min_cluster_size <- ceiling(min_cluster_fraction * nrow(data))

  if (is.null(mixtures_1D)) {
    mixtures_1D <- get_1D_mixtures_default(data, params, parallel, verbose)
  } else {
    params <- names(mixtures_1D$mixtures)
  }

  if (verbose > 0) print("Binning...")
  cutoffs <- get_1D_cutoffs(mixtures_1D$mixtures, params)
  bin_summary <- map_events_to_bins(data[,params], cutoffs)

  if (verbose > 0) print("Weighted subsampling...")
  init_mixture <- get_init(data, bin_summary, params,
                           min_cluster_size, mixture_components, verbose)
  wsub <- get_weighted_subsample(data, bin_summary, params,
                                 min_bin_size, max_bin_size, verbose)


  if (verbose > 0) print("Running bulk mixture model...")
  mixture <- bulk_weighted_gmm(data = wsub$means, k = length(init_mixture$pro),
                          params = params, weights = wsub$sizes,
                          variance_correction = wsub$variances,
                          mixture = init_mixture, verbose = (verbose >= 1))

  if(verbose > 0) print("Categorical merging...")
  cat_clusters <- categorical_merging(mixture$pro, mixture$mean, mixture$variance$sigma,
                                      cutoffs, params)

  tailor_obj <- list("mixture" = mixture, "mixtures_1D" = mixtures_1D,
                     "cat_clusters" = cat_clusters)
  class(tailor_obj) <- "tailor"
  return(tailor_obj)
}



#' @title tailor_predict
#' @description Takes as input a tailor object and some
#' data (could be the data used to learn the tailor object,
#' or some new data). Computes, for each event, the mixture
#' component from which it is most likely drawn, then maps
#' this mixture component to its
#' corresponding categorical cluster.
#' @param data A flowSet, flowFrame or a matrix containing
#' events along the rows, markers along columns.
#' @param tailor_obj A tailor object containing information
#' about mixture components and categorical clusters. Can
#' be obtained as the output of tailor.learn.
#' @param n_batch A naive implementation would need
#' nrow(data)*mixture_components memory.
#' To reduce memory usage, process data in batches.
#' @param parallel Boolean flag; if true, uses multithreading
#' to process batches in parallel.
#' For optimal runtime, if parallel = TRUE, n_batch should be
#' a multiple of the number of
#' cores available, as returned by parallel::detectCores().
#' @param verbose Boolean flag; if true, outputs timing and
#' milestone information.
#' @return Two atomic vectors of integers, one giving the
#' mixture component, and the other
#' the categorical cluster, for each event.
#' @examples
#' fileName <- system.file("extdata", "sampled_flowset_old.rda",
#'                         package = "Tailor")
#' load(fileName)
#' tailor_params <- flowCore::colnames(fs_old)[c(7:9, 11:22)]
#' tailor_obj <- tailor_learn(data = fs_old,
#'                           params = tailor_params,
#'                           mixture_components = 50)
#' tailor_pred <- tailor_predict(fs_old, tailor_obj)
#' @export
tailor_predict <- function(data, tailor_obj, n_batch = 64,
                          parallel = FALSE, verbose = FALSE)
{
  params <- colnames(tailor_obj$mixture$mean)

  data <- as_matrix(data)
  data <- data[,params]

  if (parallel)
  {
    if (verbose) { cat("Analyzing ", n_batch,
                       " batches in parallel.") }
    mapping <- tailor_map_parallel(data, tailor_obj, n_batch)
  } else {
    if (verbose) { cat("Analyzing ", n_batch,
                       " batches sequentially: ") }
    mapping <- tailor_map_sequential(data, tailor_obj, n_batch, verbose)
  }

  if(verbose) {cat("\n")}

  mix_to_clust <- tailor_obj$cat_clusters$mixture_to_cluster
  return(list(mixture_mapping = mapping,
       cluster_mapping = mix_to_clust[mapping]))
}



#' @title get_1D_mixtures
#' @description Computes 1D mixture model for each marker
#' separately, for use in binning step of tailor. It is
#' difficult to find settings which work for all datasets.
#' Therefore, it is recommended to inspect the results with
#' inspect_1D_mixtures, and run get_1D_mixtures_custom
#' for problematic markers.
#' @param data A flowSet, flowFrame or a matrix containing
#' events along the rows, markers along columns.
#' @param params A list of markers to use; must be subset
#' of colnames(data).
#' @param max_mixture Will attempt to model each marker as
#' k mixture components, for
#' 1 <= k <= max_mixture. The best k is chosen based on a
#' biased version of the Integrated Complete Likelihood (ICL).
#' @param bias_ICL Bias the ICL towards more mixture components.
#' @param sample_fraction A number between 0 and 1: the
#' fraction of data points used in the
#' calculation of 1D mixture components, to improve runtime.
#' @param parallel Boolean flag; if true, uses multithreading
#' to process markers in parallel.
#' @param verbose Boolean flag; if true, outputs timing and
#' milestone information.
#' @return A named list of 1D mixture models, giving mixture
#' proportions, means and variances for each marker.
#' @examples
#' fileName <- system.file("extdata", "sampled_flowset_old.rda",
#'                         package = "Tailor")
#' load(fileName)
#' tailor_params <- flowCore::colnames(fs_old)[c(7:9, 11:22)]
#' mixtures_1D <- get_1D_mixtures(fs_old, tailor_params)
#' @export
get_1D_mixtures <- function(data, params, max_mixture = 3,
                           use_ICL = FALSE, sample_fraction = 5e-2,
                           parallel = FALSE,
                           verbose = FALSE)
{
  data <- as_matrix(data)

  # Keep all data, or sample a subset to speed up
  if (sample_fraction < 1) {
    sample_size <- ceiling(sample_fraction * nrow(data))
    sample_size <- min(max(sample_size, 1e4), nrow(data))
    sel <- sample(nrow(data), sample_size)
    data <- data[sel,params]
  } else {
    data <- data[,params]
  }

  if (parallel) {
    if (verbose) cat("Learning", length(params),
                     "1D mixtures in parallel...")
    mixtures <- learn_1D_mixtures_parallel(data,max_mixture,
                                           prior_BIC)
  }
  else {
    if(verbose) cat("Learning", length(params),
                    "1D mixtures sequentially: ")
    mixtures <- learn_1D_mixtures_sequential(data,max_mixture,
                                             use_ICL,verbose)
    # mixtures <- learn_1D_mixtures_sequential_ICL(data,max_mixture,verbose)
  }
  if(verbose) {cat("\n")}

  return(list(mixtures = mixtures))
}


#' @title customize_1D_mixtures
#' @description After visual inspection of 1D mixtures,
#' manually specify the number of
#' mixture components to learn for some of the markers.
#' @param data A flowSet, flowFrame or a matrix containing
#' events along the rows, markers along columns.
#' @param to_customize A named list, whose names are markers,
#' and values are the number of mixture
#' components to learn for each marker.
#' @param mixtures_1D 1D mixture models, obtained
#' from get_1D_mixtures.
#' @param sample_fraction A number between 0 and 1:
#' the fraction of data points used in the
#' calculation of 1D mixture components, to improve runtime.
#' @param verbose Boolean flag; if true, outputs timing
#' and milestone information.
#' @return Updated version of mixtures_1D.
#' @examples
#' fileName <- system.file("extdata", "sampled_flowset_old.rda",
#'                         package = "Tailor")
#' load(fileName)
#' tailor_params <- flowCore::colnames(fs_old)[c(7:9, 11:22)]
#' mixtures_1D <- get_1D_mixtures(fs_old, tailor_params)
#'
#' to_customize <- list("CD127BV421" = 2)
#' mixtures_1D <- customize_1D_mixtures(fs_old, to_customize,
#'                                      mixtures_1D)
#' @export
customize_1D_mixtures <- function(data, to_customize,
                                 mixtures_1D,
                                 sample_fraction = 0.2,
                                 verbose = FALSE)
{
  data <- as_matrix(data)

  for (param in names(to_customize)) {
    temp <- get_1D_mixtures_custom(data, param,
                                   k=to_customize[[param]],
                                   sample_fraction = sample_fraction,
                                   verbose = verbose)
    mixtures_1D$mixtures[[param]] <- temp[[param]]
  }

  return(mixtures_1D)
}


#' @title inspect_1D_mixtures
#' @description Plot the result of 1D mixture model
#' calculation for visual inspection.
#' Displays, for each marker, three side-by-side plots,
#' giving a kernel density estimate
#' for the data and that marker, the Gaussian mixture,
#' and the separate mixture components,
#' respectively.
#' @param data A flowSet, flowFrame or a matrix containing
#' events along the rows, markers along columns.
#' @param mixtures_1D 1D mixture models, as produced by
#' get_1D_mixtures.
#' @param params A list of markers to use; must be subset
#' of colnames(data).
#' @return Side-by-side plots of kdes and mixture components.
#' @examples
#' fileName <- system.file("extdata", "sampled_flowset_old.rda",
#'                         package = "Tailor")
#' load(fileName)
#' tailor_params <- flowCore::colnames(fs_old)[c(7:9, 11:22)]
#' mixtures_1D <- get_1D_mixtures(fs_old, tailor_params)
#'
#' inspect_1D_mixtures(fs_old, mixtures_1D, tailor_params)
#' @export
inspect_1D_mixtures <- function(data, mixtures_1D, params)
{
  data <- as_matrix(data)

  global_kdes <- make_kdes_global(data, params)

  for (param in params)
  {
    plot_kde_vs_mixture(data, global_kdes,
                            mixtures = mixtures_1D$mixtures,
                            name = param)
  }
}


#' @title categorical_labelling
#' @description If major phenotype definitions are available
#' (e.g. "CD4 Naive"), based on only a few of the markers
#' in the panel, label each categorical cluster by one of these
#' major phenotypes.
#' @param tailor_obj A tailor object, as obtained from tailor.learn.
#' @param defs A matrix or data frame giving definitions
#' of major phenotypes.
#' @return The tailor object, with updated information
#' about categorical labels.
#' @examples
#' fileName <- system.file("extdata", "sampled_flowset_old.rda",
#'                         package = "Tailor")
#' load(fileName)
#' tailor_params <- flowCore::colnames(fs_old)[c(7:9, 11:22)]
#' tailor_obj <- tailor_learn(data = fs_old,
#'                           params = tailor_params,
#'                           mixture_components = 50)
#'
#' defs <- read.csv("~/git/R/packages/Tailor/inst/extdata/pheno_definitions.csv",
#'                  row.names = 1, stringsAsFactors = FALSE)
#' tailor_obj <- categorical_labeling(tailor_obj, defs)
#' @export
categorical_labeling <- function(tailor_obj, defs)
{
  n <- nrow(tailor_obj$cat_clusters$phenotypes)
  params <- names(tailor_obj$cat_clusters$phenotypes)
  tailor_obj$cat_clusters[["func"]] <- vector(mode = "character",
                                                length = n)
  labs <- rownames(defs)
  nam <- names(defs)
  ind <- vector(mode = "integer", length = length(nam))

  for (i in seq(length(nam))) {
    ind[i] <- which(params == nam[i])
  }

  for (i in seq(n)) {
    tailor_obj$cat_clusters$func[i] <- "UNK"

    for (j in seq(nrow(defs))) {
      match <- TRUE
      for (k in seq(ncol(defs))) {
        if (defs[j,k] == "hi" &
            tailor_obj$cat_clusters$phenotypes[i,ind[k]] == 1 |
            defs[j,k] == "lo" &
            tailor_obj$cat_clusters$phenotypes[i,ind[k]] > 1) {
          match <- FALSE
        }
      }
      if (match) {
        tailor_obj$cat_clusters$func[i] <- labs[j]
        break
      }
    }
  }

  tailor_obj$cat_clusters[["labels"]] <- compute_labels(tailor_obj$cat_clusters, defs)

  return(tailor_obj)
}


#' @title plot_tailor_majpheno
#' @description Plot a t-SNE reduction to 2d of the cluster centroids.
#' @param tailor_obj A tailor object, as obtained from tailor.learn.
#' @return A 2d reduced plot of cluster centroid,
#' color-coded by major phenotype.
#' @examples
#' fileName <- system.file("extdata", "sampled_flowset_old.rda",
#'                          package = "Tailor")
#' load(fileName)
#' tailor_params <- flowCore::colnames(fs_old)[c(7:9, 11:22)]
#' tailor_obj <- tailor_learn(data = fs_old,
#'                           params = tailor_params,
#'                           mixture_components = 50)
#'
#' defs <- read.csv("~/git/R/packages/Tailor/inst/extdata/pheno_definitions.csv",
#'                  row.names = 1, stringsAsFactors = FALSE)
#' tailor_obj <- categorical_labeling(tailor_obj, defs)
#'
#' plot_tailor_majpheno(tailor_obj)
#' @export
plot_tailor_majpheno <- function(tailor_obj)
{
  pro <- tailor_obj$mixture$pro
  map <- tailor_obj$cat_clusters$mixture_to_cluster

  res <- get_tsne_clusters(tailor_obj)
  res$phenotype <- as.factor(tailor_obj$cat_clusters$func)

  cluster_ids <- seq_len(nrow(tailor_obj$cat_clusters$centers))
  res$logsize <- log(vapply(cluster_ids,
                            function(x) { sum(pro[which(map == x)]) },
                            numeric(1) ))
  res$logsize <- res$logsize * 9 / mean(res$logsize)

  g <- ggplot(res, aes(x=.data$tsne_1, y=.data$tsne_2)) +
    geom_point(aes(color = .data$phenotype, size = .data$logsize),
               alpha = 0.8) +
    scale_color_brewer(palette = "Paired") +
    guides(size = FALSE) +
    theme_bw()

  return(g)
}



#' @title plot_tailor_fluorescence
#' @description Plot a t-SNE reduction to 2d of the cluster centroids,
#' color-coded by mean fluorescence intensity for each parameter.
#' @param tailor_obj A tailor object, as obtained from tailor.learn.
#' @param midpoint Numeric, to be used as the middle of the scale for intensity plots.
#' @return A grid of 2d reduced plots of cluster centroids, one for each
#' parameter.
#' @examples
#' fileName <- system.file("extdata", "sampled_flowset_old.rda",
#'                        package = "Tailor")
#' load(fileName)
#' tailor_params <- flowCore::colnames(fs_old)[c(7:9, 11:22)]
#' tailor_obj <- tailor_learn(data = fs_old,
#'                           params = tailor_params,
#'                           mixture_components = 50)
#'
#' plot_tailor_fluorescence(tailor_obj)
#' @export
plot_tailor_fluorescence <- function(tailor_obj, midpoint = 1.5)
{
  pro <- tailor_obj$mixture$pro
  map <- tailor_obj$cat_clusters$mixture_to_cluster
  cen <- tailor_obj$cat_clusters$centers

  res <- get_tsne_clusters(tailor_obj)

  cluster_ids <- seq_len(nrow(tailor_obj$cat_clusters$centers))
  res$logsize <- log(vapply(cluster_ids,
                            function(x) { sum(pro[which(map == x)]) },
                            numeric(1) ))
  res$logsize <- res$logsize * 7 / mean(res$logsize)

  plot_list <- list()
  params = colnames(tailor_obj$mixture$mean)

  mfi_amp <- vapply(params, function(x) { max(cen[,x]) - min(cen[,x]) },
                    numeric(1)  )
  max_amplitude  <- names(which.max(mfi_amp))

  for (param in params)
  {
    res$color <- cen[,param]
    g <- ggplot(res, aes(x=.data$tsne_1, y = .data$tsne_2)) +
      geom_point(aes(color = .data$color, size = .data$logsize), alpha = 0.8) +
      scale_color_gradient2(low = "green", mid = "yellow",
                            high = "red", midpoint = midpoint,
                            guide = guide_colorbar(title = "Mean Fluorescence Intensity",
                                                   direction = "horizontal",
                                                   title.position = "top")) +
      ggtitle(param) +
      guides(size = FALSE) +
      theme_bw() +
      theme(axis.title.x=element_blank(),
            axis.text.x=element_blank(),
            axis.ticks.x=element_blank(),
            axis.title.y=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks.y=element_blank())

    if (param == max_amplitude) { color_legend <- get_legend(g) }

    plot_list[[param]] <- g + theme(legend.position = "none")
  }

  plot_list[["legend"]] <- color_legend

  ncol <- ceiling(sqrt(length(plot_list)))
  return(gridExtra::grid.arrange(grobs = plot_list, ncol = ncol))
}




#' @title plot_kdes_global
#' @description Plot 1D kdes of a dataset for visual inspection.
#' @param data A flowSet, flowFrame or a matrix containing events
#' along the rows, markers along columns.
#' @param data_overlay Some different data for comparison.
#' @param params A list of markers to use; must be subset of
#' colnames(data).
#' @return Plots of kernel density estimate for each chosen parameter.
#' @examples
#' fileName <- system.file("extdata", "sampled_flowset_old.rda",
#'                          package = "Tailor")
#' load(fileName)
#' tailor_params <- flowCore::colnames(fs_old)[c(7:9, 11:22)]
#' plot_kdes_global(fs_old, params = tailor_params)
#' @export
plot_kdes_global <- function (data, data_overlay = NULL, params = NULL)
{
  data <- as_matrix(data)
  if (is.null(params)) params <- colnames(data)

  global_kdes <- make_kdes_global(data, params)

  if (!is.null(data_overlay)) {
    data_overlay <- as_matrix(data_overlay)
    overlay_kdes <- make_kdes_global(data_overlay, params)
  }

  plot_list <- list()
  for (param in params) {
    df <- data.frame(global_kdes[[param]])
    g <- ggplot() +
      geom_line(df, mapping = aes(x = .data$x, y = .data$y)) +
      labs(title = param, x = "", y = "") +
      theme_bw()

    if (!is.null(data_overlay)) {
      df_overlay <- data.frame(overlay_kdes[[param]])
      g <- g + geom_line(df_overlay, mapping = aes(x = .data$x, y = .data$y),
                         color = "red")
    }

    plot_list[[param]] <- g
  }
  ncol <- min(5, ceiling(sqrt(length(params))))
  return(gridExtra::grid.arrange(grobs = plot_list, ncol = ncol))
}


#' @title make_cluster_phenobars
#' @description Visualize summary of cluster information, via boxplots showing means
#' and standard deviations of mean fluorescence intensity. Only markers with at least
#' one cutoff are displayed. The boxes are color-coded by mfi relative to the cutoff.
#' @param tailor_obj A tailor object, as obtained from tailor.learn.
#' @param cluster_ids A list of cluster indices to visualize. For example, those deemed
#' differentially expressed in case/controls by some statistical test.
#' @return Boxplots summarizing means and standard deviations for selected clusters.
#' @examples
#' fileName <- system.file("extdata", "sampled_flowset_old.rda",
#'                        package = "Tailor")
#' load(fileName)
#' tailor_params <- flowCore::colnames(fs_old)[c(7:9, 11:22)]
#' tailor_obj <- tailor_learn(data = fs_old,
#'                           params = tailor_params,
#'                           mixture_components = 50)
#'
#' make_cluster_phenobars(tailor_obj, cluster_ids = c(1,2))
#' @export
make_cluster_phenobars <- function(tailor_obj, cluster_ids, cluster_names = NULL) {
  plot_list <- list()
  cat <- tailor_obj$cat_clusters

  if(is.null(cluster_names))
    cluster_names <- vapply(cluster_ids, function(id) paste("Cluster", id),
                            character(1))

  params <- colnames(cat$phenotypes)
  ymin <- min(0, min(cat$centers[cluster_ids,params] - cat$sd[cluster_ids,params]))
  ymax <- max(cat$centers[cluster_ids,params] + cat$sd[cluster_ids,params])

  midpt <- c()
  cutoffs <- tailor_obj$cat_clusters$cutoffs
  for (param in params) {
    l <- length(cutoffs[[param]])
    if (l %% 2 == 1) {
      midpt <- c(midpt, cutoffs[[param]][l%/%2 + 1])
    } else {
      midpt <- c(midpt, 0.5*(cutoffs[[param]][l%/%2] + cutoffs[[param]][l%/%2 + 1]))
    }
  }

  for (i in seq_along(cluster_ids)) {
    id <- cluster_ids[i]
    df <- data.frame("name" = params,
                     "mfi" = tailor_obj$cat_clusters$centers[id,params],
                     "sd" = tailor_obj$cat_clusters$sds[id,params],
                     "midpt" = midpt)

    p <- ggplot(df) +
      geom_bar(aes(x = name, y = mfi, fill = mfi-midpt), color = "black", stat = "identity") +
      geom_linerange(aes(x = name, ymin = mfi - sd, ymax = mfi + sd)) +
      scale_fill_gradient2(low = "green", mid = "yellow", "high" = "red") +
      scale_y_continuous(limits = c(ymin, ymax), expand = expansion(mult = c(0.02,0.02))) +
      coord_flip() +
      labs(x = NULL, y = NULL, title = cluster_names[i]) +
      theme_bw() +
      theme(legend.position = "none")

    plot_list[[i]] <- p
  }

  ncol <- min(5, ceiling(sqrt(length(cluster_ids))))
  return(gridExtra::grid.arrange(grobs = plot_list, ncol = ncol))
}




