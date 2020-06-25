###############################
# Functions exposed to the user
###############################

#' @import Rtsne
#' @import KernSmooth
#' @import cluster
#' @import wadeTools
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

#' @title tailor.learn
#' @description This function learns a tailor model from input data. It
#' computes a preliminary binning of the data, then computes a mixture model using
#' a weighted version of the expectation-maximization (EM) algorithm,
#' and finally merges mixture components which are positive/negative for the same
#' markers, using adaptive thresholds.
#' @param data A flowSet, flowFrame or a matrix containing events along the rows, markers along columns.
#' @param params A list of markers to use; must be subset of colnames(data).
#' @param mixture_components The number of mixture components to learn. Some of these
#' are eventually merged, so it's a good idea to choose a number slightly larger than
#' the number of clusters you expect to get.
#' @param min_bin_size Bins with fewer events than this threshold are considered outliers, and are
#' ignored during the weighted EM algorithm. These events can still be assigned to clusters during
#' the prediction phase.
#' @param max_bin_size Bins with more events than this threshold are split, to ensure that the
#' weighted EM algorithm closely approximates a run of vanilla EM on the entire dataset.
#' @param mixtures_1D Pre-computed 1D mixture models, to be used for binning.
#' These are computed from scratch if not provided.
#' @param seed A seed for the random number generator, to be used in the initialization of the bin
#' centers.
#' @param do_tsne Boolean flag; if true, computes t-SNE reduction to 2D of the bin centers,
#' for use in visualizations.
#' @param do_variance_correction Boolean flag; if true, corrects the variances of the mixture
#' components found in weighted EM, to account for intra-bin variance. TRUE is strongly recommended.
#' @param parallel Boolean flag; if true, uses multithreading to speed up computation.
#' @param verbose If > 0, outputs milestone information. If >=1, also outputs information about
#' running utilities. If >1, debugging mode.
#' @return A tailor object containing:
#' \describe{
#'   \item{fit}{The tailor model, a named list containing the mixture proportions, means and variances
#'   of all mixture components.}
#'   \item{mixtures_1D}{}
#'   \item{cat_clusters}{A named list containing information about the categorical clusters found by
#'   the model: phenotype, cluster centers, and a mapping from mixture components to categorical clusters.}
#'   \item{tsne_centers}{Optional: a dimensional reduction to 2D of the bin centers.}
#' }
#' @export
tailor.learn = function(data, params = NULL,
                        mixture_components = 200,
                        min_bin_size = NULL, max_bin_size = NULL,
                        mixtures_1D = NULL,
                        seed = NULL,
                        do_tsne = FALSE,
                        do_variance_correction = TRUE,
                        parallel = FALSE,
                        verbose = 0.5)
{
  start.time = Sys.time()

  if(is(data, "flowSet"))
  {
    data = suppressWarnings(as(data, "flowFrame"))
    flowCore::exprs(data) = flowCore::exprs(data)[,which(flowCore::colnames(data) != "Original")]
  }

  if(is(data, "flowFrame")) data = flowCore::exprs(data)

  if (is.null(params)) params = colnames(data)
  d = length(params)

  if (!is.null(seed)) set.seed(seed)

  if (is.null(min_bin_size))
  {
    min_bin_size = max(5,ceiling(nrow(data) / 1e5))
  }
  if(is.null(max_bin_size))
  {
    max_bin_size = max(50,ceiling(nrow(data) / 1e3))
  }

  # preliminary binning

  if (is.null(mixtures_1D))
  {
    # Learn 1D mixture model for each variable
    if (verbose > 0) {print("Getting 1D mixtures...")}

    sample_fraction = 0.5
    if (nrow(data) > 5e5) sample_fraction = 0.2
    if (nrow(data) > 1e7) sample_fraction = 0.1

    mixtures_1D = get_1D_mixtures(data[,params], params, max_mixture = 3,
                                    sample_fraction = sample_fraction, seed = seed,
                                    prior_BIC = 100,
                                    parallel = parallel,
                                    verbose = (verbose >= 1))

    # Merge modes whose mean is small enough.
    # These are likely compensation artifacts.
    mixtures_1D$to_merge = find_to_merge(mixtures_1D, params,
                             negative_threshold = 0.5,
                             verbose = (verbose == 1))
  } else
  {
    # If given path to pre-computed mixtures, load them
    params = names(mixtures_1D$mixtures)
  }


  # Up-sample: assign each data point to most likely component of the model
  if (verbose > 0) { print("Assigning data to 1D mixtures...")}

  mapping = mapping_from_mixtures(data[,params], mixtures_1D$mixtures, mixtures_1D$to_merge,
                                  params,
                                  parallel = parallel, verbose = (verbose >= 1))
  phenobin = phenobin_label(mapping)


  if (verbose > 0) { print("Preparing initialization of bulk mixture model...")}

  # Parse phenobins
  phenobin_summary = get_phenobin_summary(phenobin)

  large_bins = which(phenobin_summary$bins_sorted >= min_bin_size)
  bins = as.integer(names(phenobin_summary$bins_sorted)[large_bins])
  sizes = as.vector(phenobin_summary$bins_sorted)[large_bins]
  populous = which(phenobin_summary$predictions %in% bins)

  if(verbose > 1)
  {
    cat("Total bin number before splitting: ", length(phenobin_summary$bins_sorted), "\n")
    cat("Maximum bin size: ", sizes[1], "\n")
  }



  phenobin_parameters = find_phenobin_mean(data = data,
                                           predictions = phenobin_summary$predictions,
                                           bins = bins, sizes = sizes,
                                           params = params,
                                           selected = populous,
                                           split_threshold = max_bin_size,
                                           compute_var = do_variance_correction,
                                           seed = seed,
                                           parallel = parallel,
                                           verbose = (verbose >= 1))

  nonzero = which(phenobin_parameters$sizes > 0)
  repr_sizes = phenobin_parameters$sizes[nonzero]
  repr_means = phenobin_parameters$means[nonzero,]
  repr_variances = phenobin_parameters$variances[nonzero,,]


  if (verbose > 1)
  {
    cat("Total bin number after selecting and splitting: ", nrow(repr_means), "\n")
    cat("Maximum bin size: ", max(repr_sizes), "\n")
  }



  # Prepare initialization of bulk mixture model
  candidate_bins = c(1:mixture_components)
  init_bins = which(sizes[candidate_bins] > 3 * min_bin_size) # heuristic: strive to improve
  lost = length(candidate_bins) - length(init_bins)
  if (lost > 0) cat("Warning: dropped", lost, "mixture components due to too few events at initialization.\n")

  populous = which(phenobin_summary$predictions %in% bins[init_bins])
  phenobin_parameters = find_phenobin_mean(data = data,
                                           predictions = phenobin_summary$predictions,
                                           bins = bins[init_bins], sizes = sizes[init_bins],
                                           params = params,
                                           selected = populous,
                                           split_threshold = NULL,
                                           compute_var = TRUE,
                                           seed = seed,
                                           parallel = parallel,
                                           verbose = (verbose >= 1))


  mixture = list()
  mixture$pro = sizes[init_bins]
  mixture$mean = phenobin_parameters$means
  mixture$variance$sigma = array(NaN, c(d,d,mixture_components))

  for (component in init_bins)
  {
    sigma = phenobin_parameters$variances[component,,]
    mixture$variance$sigma[,,component] = sigma

    de = det(sigma)
    cat(component, de, mixture$pro[component], "\n")
  }


  # Learn bulk mixture model on weighted subsample
  if (verbose > 0) { print("Running bulk mixture model...")}

  if (do_variance_correction)
  {
    fit = bulk_weighted_gmm(data = repr_means,
                            k = mixture_components,
                            params = params,
                            weights = repr_sizes,
                            initialize = mixture,
                            regularize_variance = TRUE,
                            variance_correction = repr_variances,
                            verbose = (verbose >= 1))
  }
  else
  {
    fit = bulk_weighted_gmm(data = repr_means,
                            k = mixture_components,
                            params = params,
                            weights = repr_sizes,
                            initialize = mixture,
                            regularize_variance = TRUE,
                            variance_correction = NULL,
                            verbose = (verbose >= 1))
  }

  if(verbose > 0) print("Categorical merging...")

  # Use the 1D mixture models to decide on +/- cutoffs
  cutoffs = get_1D_cutoffs(mixtures_1D$mixtures, mixtures_1D$to_merge, params)

  # Categorical merging
  cat_clusters = categorical_merging(fit$mixture$pro,
                             fit$mixture$mean,
                             cutoffs,
                             params)

  if (do_tsne)
  {
    if (verbose > 0) { print("Reducing phenobin centers to 2D...") }
    tsne_centers = get_tsne_centers(data = repr_means, seed = seed, probs = fit$event_probabilities)
  }


  if (do_tsne)
  {
    tailor_obj = list("mixture" = fit$mixture, "mixtures_1D" = mixtures_1D, "cat_clusters" = cat_clusters,
         "tsne_centers" = tsne_centers)
  }
  else
  {
    tailor_obj = list("mixture" = fit$mixture, "mixtures_1D" = mixtures_1D, "cat_clusters" = cat_clusters)
  }

  if (verbose > 0) print(Sys.time() - start.time)

  class(tailor_obj) = "tailor"
  tailor_obj
}



#' @title tailor.predict
#' @description Takes as input a tailor object and some data (could be the data used
#' to learn the tailor object, or some new data). Computes, for each event, the mixture
#' component from which it is most likely drawn, then maps this mixture component to its
#' corresponding categorical cluster.
#' @param data A flowSet, flowFrame or a matrix containing events along the rows, markers along columns.
#' @param tailor_obj A tailor object containing information about mixture components
#' and categorical clusters. Can be obtained as the output of tailor.learn.
#' @param n_batch A naive implementation would need nrow(data)*mixture_components memory.
#' To reduce memory usage, process data in batches.
#' @param parallel Boolean flag; if true, uses multithreading to process batches in parallel.
#' For optimal runtime, if parallel = TRUE, n_batch should be a multiple of the number of
#' cores available, as returned by parallel::detectCores().
#' @param verbose Boolean flag; if true, outputs timing and milestone information.
#' @return Two atomic vectors of integers, one giving the mixture component, and the other
#' the categorical cluster, for each event.
#' @export
tailor.predict = function(data, tailor_obj, n_batch = 64,
                          parallel = FALSE, verbose = FALSE)
{
  if(is(data, "flowSet"))
  {
    data = suppressWarnings(as(data, "flowFrame"))
    flowCore::exprs(data) = flowCore::exprs(data)[,which(flowCore::colnames(data) != "Original")]
  }

  if(is(data,"flowFrame")) data = flowCore::exprs(data)

  start.time = Sys.time()
  n = nrow(data)
  k = length(tailor_obj$mixture$pro)
  logpro = log(tailor_obj$mixture$pro)
  params = colnames(tailor_obj$mixture$mean)
  data = data[,params]

  mapping = integer(length = n)

  # n*k may be an unreasonable amount of memory to use, so break into batches
  batch_size = ceiling(n/n_batch)

  if (parallel)
  {
    if (verbose) { cat("Analyzing ", n_batch, " batches in parallel.") }
    #setup parallel backend to use many processors
    cores=detectCores()
    clust <- makeCluster(cores[1])
    registerDoParallel(clust)

    data_list = list()
    for (batch in c(1:n_batch))
    {
      start_batch = (batch - 1) * batch_size + 1
      end_batch = batch * batch_size
      if (batch == n_batch) { end_batch = n }

      data_list[[batch]] = data[c(start_batch:end_batch),params]
    }

    rm(data)

    mapping = foreach(batch_data=iter(data_list), .combine = c, .packages = c("mvtnorm")) %dopar%
      {
        posteriors = matrix(0, nrow = nrow(batch_data), ncol = k)

        for (cl in seq(k))
        {
          mean = tailor_obj$mixture$mean[cl,]
          sigma = tailor_obj$mixture$variance$sigma[,,cl]
          posteriors[,cl] = logpro[cl] + dmvnorm(batch_data, mean, sigma, log = TRUE)
        }

        # Assign each datapoint to the cluster of maximum probability
        result = apply(posteriors, 1, which.max)
        rm(posteriors)
        gc()

        result
      }

    #stop cluster
    stopCluster(clust)
  } else
  {
    if (verbose) { cat("Analyzing ", n_batch, " batches sequentially: ") }

    for (batch in c(1:n_batch))
    {
      if(verbose) {cat(batch, " ")}
      start_batch = (batch - 1) * batch_size + 1
      end_batch = batch * batch_size
      if (batch == n_batch) { end_batch = n }


      batch_data = data[c(start_batch:end_batch),params]

      # For each cluster, compute the probability that data in current batch
      # are drawn from it
      posteriors = matrix(0, nrow = end_batch - start_batch + 1, ncol = k)


      for (cl in seq(k))
      {
        weight = tailor_obj$mixture$pro[cl]
        mean = tailor_obj$mixture$mean[cl,]
        sigma = tailor_obj$mixture$variance$sigma[,,cl]

        posteriors[,cl] = weight * dmvnorm(batch_data, mean, sigma)
      }


      # Assign each datapoint to the cluster of maximum probability
      mapping[c(start_batch:end_batch)] = apply(posteriors, 1, which.max)

      gc()
    }
  }

  if(verbose) {cat("\n")}
  if(verbose) { print(Sys.time() - start.time) }

  list(mixture_mapping = mapping,
       cluster_mapping = tailor_obj$cat_clusters$mixture_to_cluster[mapping])
}



#' @title get_1D_mixtures
#' @description Computes 1D mixture model for each marker separately, for use in binning step
#' of tailor. It is difficult to find settings which work for all datasets. Therefore, it is
#' recommended to inspect the results with inspect_1D_mixtures, and run get_1D_mixtures_custom
#' for problematic markers.
#' @param data A flowSet, flowFrame or a matrix containing events along the rows, markers along columns.
#' @param params A list of markers to use; must be subset of colnames(data).
#' @param max_mixture Will attempt to model each marker as k mixture components, for
#' 1 <= k <= max_mixture. The best k is chosen based on a modified version of the Bayesian
#' Information Criterion (BIC).
#' @param prior_BIC Make this larger to favor a smaller number of mixture components.
#' @param sample_fraction A number between 0 and 1: the fraction of data points used in the
#' calculation of 1D mixture components, to improve runtime.
#' @param seed A seed for the random generator, used for subsampling and for the initialization
#' of the 1D mixture components. The results are qualitatively the same, irrespective of the
#' seed chosen.
#' @param parallel Boolean flag; if true, uses multithreading to process markers in parallel.
#' @param verbose Boolean flag; if true, outputs timing and milestone information.
#' @return A named list of 1D mixture models, giving mixture proportions, means and variances
#' for each marker.
#' @export
get_1D_mixtures = function(data, params, max_mixture = 3,
                           prior_BIC = 600, sample_fraction = 0.2, seed = NULL,
                           parallel = FALSE,
                           verbose = FALSE)
{
  if(is(data, "flowSet"))
  {
    data = suppressWarnings(as(data, "flowFrame"))
    flowCore::exprs(data) = flowCore::exprs(data)[,which(flowCore::colnames(data) != "Original")]
  }

  if(is(data,"flowFrame")) data = flowCore::exprs(data)

  start.time = Sys.time()

  data = data[,params]

  # Keep all data, or sample a subset to speed up
  if (!is.null(seed)) set.seed(seed)
  if (sample_fraction == 1)
  {
    sel = c(1:nrow(data))
  } else
  {
    sample_size = ceiling(sample_fraction * nrow(data))
    sel = sample(nrow(data), sample_size)
  }

  if (parallel)
  {
    #setup parallel backend to use many processors
    cores=detectCores()
    cl <- makeCluster(cores[1])
    registerDoParallel(cl)


    if (verbose) cat("Learning", length(params), "1D mixtures in parallel...")

    data_param = NULL # unnecessary definition, but R CMD check complains without it

    fit_list = foreach (data_param = iter(data[sel,], by = 'col'), .packages = c("mclust")) %dopar%
      {
        param = colnames(data_param)[1]

        set.seed(seed)

        # Use Bayesian information criterion to choose best k
        BIC = mclustBIC(data_param, G = 1:max_mixture, modelNames = "V", verbose = FALSE)

        # A tweak to favor smaller k
        for (k in c(1:max_mixture))
        {
          BIC[k] = BIC[k] - prior_BIC*k*log(length(data_param), base = 2)
        }

        # Fit the model with the chosen k
        fit = Mclust(data_param, x=BIC, verbose = FALSE)
        fit$parameters
      }
    names(fit_list) = params

    #stop cluster
    stopCluster(cl)
  } else
  {
    fit_list = list()

    if(verbose) cat("Learning", length(params), "1D mixtures sequentially: ")

    for (param in params)
    {
      if (verbose) {cat(param, " ")}

      dat = data[sel,param]

      # set.seed(seed)

      # Use Bayesian information criterion to choose best k
      BIC = mclustBIC(dat, G = 1:max_mixture, modelNames = "V", verbose = FALSE)

      # A tweak to favor smaller k
      for (k in c(1:max_mixture))
      {
        BIC[k] = BIC[k] - prior_BIC*k*log(length(dat), base = 2)
      }

      # Fit the model with the chosen k
      fit = Mclust(dat, x=BIC, verbose = FALSE)
      fit_list[[param]] = fit$parameters
    }
  }

  gc()

  end.time = Sys.time()
  if(verbose) {cat("\n")}
  if(verbose) {print(end.time - start.time)}

  list(mixtures = fit_list)
}


#' @title customize_1D_mixtures
#' @description After visual inspection of 1D mixtures, manually specify the number of
#' mixture components to learn for some of the markers.
#' @param data A flowSet, flowFrame or a matrix containing events along the rows, markers along columns.
#' @param to_customize A named list, whose names are markers, and values are the number of mixture
#' components to learn for each marker.
#' @param mixtures_1D 1D mixture models, obtained from get_1D_mixtures.
#' @param sample_fraction A number between 0 and 1: the fraction of data points used in the
#' calculation of 1D mixture components, to improve runtime.
#' @param seed A seed for the random generator, used for subsampling and for the initialization
#' of the 1D mixture components. The results are qualitatively the same, irrespective of the
#' seed chosen.
#' @param verbose Boolean flag; if true, outputs timing and milestone information.
#' @return Updated version of mixtures_1D.
#' @export
customize_1D_mixtures = function(data, to_customize,
                                 mixtures_1D,
                                 sample_fraction = 0.2,
                                 seed = NULL,
                                 verbose = FALSE)
{
  if(is(data, "flowSet"))
  {
    data = suppressWarnings(as(data, "flowFrame"))
    flowCore::exprs(data) = flowCore::exprs(data)[,which(flowCore::colnames(data) != "Original")]
  }

  if(is(data,"flowFrame")) data = flowCore::exprs(data)

  for (param in names(to_customize))
  {
    mixtures_1D$mixtures[[param]] = get_1D_mixtures_custom(data, param,
                                                           k=to_customize[[param]],
                                                           sample_fraction = sample_fraction,
                                                           seed = seed,
                                                           verbose = verbose)[[param]]
  }

  mixtures_1D
}


#' @title inspect_1D_mixtures
#' @description Plot the result of 1D mixture model calculation for visual inspection.
#' Displays, for each marker, three side-by-side plots, giving a kernel density estimate
#' for the data and that marker, the Gaussian mixture, and the separate mixture components,
#' respectively.
#' @param data A flowSet, flowFrame or a matrix containing events along the rows, markers along columns.
#' @param mixtures_1D 1D mixture models, as produced by get_1D_mixtures.
#' @param params A list of markers to use; must be subset of colnames(data).
#' @export
inspect_1D_mixtures = function(data, mixtures_1D, params)
{
  if(is(data, "flowSet"))
  {
    data = suppressWarnings(as(data, "flowFrame"))
    flowCore::exprs(data) = flowCore::exprs(data)[,which(flowCore::colnames(data) != "Original")]
  }

  if(is(data,"flowFrame")) data = flowCore::exprs(data)

  global_kdes = make_kdes_global(data, params)

  for (param in params)
  {
    plot_kde_vs_mixture(data, global_kdes, mixtures = mixtures_1D$mixtures, name = param)
  }
}


#' @title categorical_labelling
#' @description If major phenotype definitions are available (e.g. "CD4 Naive"),
#' based on only a few of the markers
#' in the panel, label each categorical cluster by one of these major phenotypes.
#' @param tailor_obj A tailor object, as obtained from tailor.learn.
#' @param defs A matrix or data frame giving definitions of major phenotypes.
#' @return The tailor object, with updated information about categorical labels.
#' @export
categorical_labeling = function(tailor_obj, defs)
{
  n = length(tailor_obj$cat_clusters$phenotypes)
  params = colnames(tailor_obj$mixture$mean)
  tailor_obj$cat_clusters[["labels"]] = vector(mode = "character", length = n)
  labs = rownames(defs)
  nam = names(defs)
  ind = vector(mode = "integer", length = length(nam))

  for (i in seq(length(nam)))
  {
    ind[i] = which(params == nam[i])
  }

  for (i in seq(n))
  {
    tailor_obj$cat_clusters$labels[i] = "UNK"

    for (j in seq(nrow(defs)))
    {
      match = TRUE
      for (k in seq(ncol(defs)))
      {
        if (defs[j,k] == "hi" & substr(tailor_obj$cat_clusters$phenotypes[i],ind[k],ind[k]) == "-" |
            defs[j,k] == "lo" & substr(tailor_obj$cat_clusters$phenotypes[i],ind[k],ind[k]) == "+")
        {
          match = FALSE
        }
      }
      if (match)
      {
        tailor_obj$cat_clusters$labels[i] = labs[j]
        break
      }
    }
  }
  tailor_obj
}


#' @title plot_tsne_clusters
#' @description Plot a t-SNE reduction to 2d of the cluster centroids.
#' @param tailor_obj A tailor object, as obtained from tailor.learn.
#' @param tailor_pred A list of cluster assignments, as obtained from tailor.predict.
#' @param seed A seed for the random generator, used in the initialization of t-SNE.
#' @return A 2d reduced plot of cluster centroid, color-coded by major phenotype.
#' @export
plot_tsne_clusters = function(tailor_obj, tailor_pred, seed = NULL)
{
  res = get_tsne_clusters(tailor_obj, seed = seed)
  res$phenotype = as.factor(tailor_obj$cat_clusters$labels)
  res$logsize = log(tabulate(tailor_pred$cluster_mapping))
  res$logsize = res$logsize * 9 / mean(res$logsize)

  g = ggplot(res, aes(x=.data$tsne_1, y=.data$tsne_2)) +
    geom_point(aes(color = .data$phenotype, size = .data$logsize)) +
    scale_color_brewer(palette = "Paired")

  g
}




#' @title plot_tsne_bin_centers
#' @description Plot a t-SNE reduction to 2d of the cluster centroids.
#' @param tailor_obj A tailor object, as obtained from tailor.learn.
#' @return A 2d reduced plot of bin centers, which can be seen as representatives for the cell
#' population. Color coded by major phenotype.
#' @export
plot_tsne_bin_centers = function(tailor_obj)
{
  binc = tailor_obj$tsne_centers
  binc$phenotype = tailor_obj$cat_clusters$mixture_to_cluster[binc$cluster]
  binc$phenotype = tailor_obj$cat_clusters$labels[binc$phenotype]

  g = ggplot(binc, aes(x=.data$tsne_1, y=.data$tsne_2)) +
    geom_point(aes(color = .data$phenotype)) +
    scale_color_brewer(palette = "Paired")

  g
}

#' @title plot_tsne_global
#' @description Subsample the entire dataset, and compute a 2D t-SNE embedding.
#' Plot the embedding twice, color coded by cluster membership and major phenotype, respectively.
#' @param data A flowSet, flowFrame or a matrix containing events along the rows, markers along columns.
#' @param tailor_obj A tailor object, as obtained from tailor.learn.
#' @param tailor_pred A list of cluster assignments, as obtained from tailor.predict.
#' @param defs Definitions of major phenotypes.
#' @param seed A seed for the random generator, used in the initialization of t-SNE.
#' @return A 2d embeddings, color-coded by cluster membershi and major phenotype.
#' @export
plot_tsne_global = function(data, tailor_obj, tailor_pred, defs, seed = NULL)
{
  # Preprocess input
  if(is(data, "flowSet"))
  {
    data = suppressWarnings(as(data, "flowFrame"))
    flowCore::exprs(data) = flowCore::exprs(data)[,which(flowCore::colnames(data) != "Original")]
  }

  if(is(data,"flowFrame")) data = flowCore::exprs(data)

  n_items = nrow(data)

  params = colnames(tailor_obj$mixture$mean)

  # Subsample, because t-SNE is slow
  if(!is.null(seed)) set.seed(seed)

  if (n_items > 1e4)
  {
    sel = sample(n_items, 1e4)
  } else
  {
    sel = c(1:n_items)
  }
  data = data[sel,params]
  perplexity = min((length(sel) - 1)/3, 30)

  # Learn t-SNE embedding and add column for cluster membership
  tsne = Rtsne(data, perplexity = perplexity)$Y
  colnames(tsne) = c("tsne_1", "tsne_2")
  tsne = data.frame(tsne)
  tsne$cluster = tailor_pred$cluster_mapping[sel]
  tsne$cluster_pheno = tailor_obj
  tsne$cluster = as.factor(tsne$cluster)


  # Decide on major phenotype of each subsampled data point
  cutoffs = get_1D_cutoffs(tailor_obj$mixtures_1D$mixtures, tailor_obj$mixtures_1D$to_merge,
                           params)
  data_cut = data
  for (param in params)
  {
    cutoff = cutoffs[[param]]
    data_cut[,param] = "hi"
    data_cut[which(data[,param] < cutoff), param] = "lo"
  }

  phenotype = character(length = length(sel))
  phenotype[c(1:length(sel))] = "Other"
  for (pheno in rownames(defs))
  {
    relevant = names(defs)[which(defs[pheno,] != "dc")]
    data_relevant = data_cut[,relevant]
    phenotype[which( apply(data_relevant, 1,
                           function(x) identical(as.character(x),
                                                 as.character(defs[pheno,relevant]))))] = pheno
  }

  tsne$phenotype = as.factor(phenotype)


  # Display plots
  mycolors = c(brewer.pal(name="Set1", n = 9), brewer.pal(name="Set3", n = 12))
  mycolors = mycolors[c(1:2, 10:21)]

  g = ggplot(tsne, aes(x=.data$tsne_1, y=.data$tsne_2, color = .data$phenotype)) +
      geom_point() +
      scale_color_manual(values = mycolors)

  print(g)



  # plot_list = list()
  # g = ggplot(tsne, aes(x=.data$tsne_1, y=.data$tsne_2)) +
  #   geom_point(aes(color = .data$phenotype)) +
  #   scale_color_manual(values = mycolors)
  # plot_list[[1]] = g
  #
  # xlim = ggplot_build(g)$layout$panel_scales_x[[1]]$range$range
  # ylim = ggplot_build(g)$layout$panel_scales_y[[1]]$range$range
  #
  # len = length(tailor_obj$cat_clusters$labels)
  # nplot = ceiling(len/7)
  #
  # for (plot_index in seq(1:nplot))
  # {
  #   start = 7 * (plot_index - 1) + 1
  #   end   = min(7 * plot_index, len)
  #   sel = which(tsne$cluster %in% c(start:end))
  #
  #   g = ggplot(tsne[sel,], aes(x=.data$tsne_1, y=.data$tsne_2)) +
  #     geom_point(aes(color = .data$cluster)) +
  #     scale_color_brewer(palette = "Dark2") +
  #     scale_x_continuous(limits = xlim) +
  #     scale_y_continuous(limits = ylim)
  #
  #   plot_list[[plot_index + 1]] = g
  # }
  #
  # ncol = ceiling(sqrt(nplot))
  # gridExtra::grid.arrange(grobs = plot_list, ncol = ncol)
}


#' @title plot_kdes_global
#' @description Plot 1D kdes of a dataset for visual inspection.
#' @param data A flowSet, flowFrame or a matrix containing events along the rows, markers along columns.
#' @param params A list of markers to use; must be subset of colnames(data).
#' @export
plot_kdes_global = function(data, params)
{
  # Preprocess input
  if(is(data, "flowSet"))
  {
    data = suppressWarnings(as(data, "flowFrame"))
    flowCore::exprs(data) = flowCore::exprs(data)[,which(flowCore::colnames(data) != "Original")]
  }
  if(is(data,"flowFrame")) data = flowCore::exprs(data)

  global_kdes = make_kdes_global(data, params)

  plot_list = list()
  for (param in params)
  {
    df = data.frame(global_kdes[[param]])
    g = ggplot(df, aes(x=.data$x, y=.data$y)) +
      geom_line() +
      labs(title = param, x = "", y = "") +
      theme_bw()

    plot_list[[param]] = g
  }

  ncol = min(5, ceiling(sqrt(length(params))))
  gridExtra::grid.arrange(grobs = plot_list, ncol = ncol)
}



