###############################
# Functions exposed to the user
###############################

#' @import Rtsne
#' @import KernSmooth
#' @import cluster
#' @import wadeTools
#' @import mvtnorm 
#' @import ggplot2
#' @import mclust
#' @import foreach
#' @import doParallel
NULL

#' @title tailor.learn
#' @description This function learns a tailor model from input data. It
#' computes a preliminary binning of the data, then computes a mixture model using
#' a weighted version of the expectation-maximization (EM) algorithm,
#' and finally merges mixture components which are positive/negative for the same
#' markers, using adaptive thresholds.
#' @param data A matrix containing events along the rows, markers along columns.
#' @param params A list of markers to use; must be subset of colnames(data).
#' @param mixture_components The number of mixture components to learn. Some of these 
#' are eventually merged, so it's a good idea to choose a number slightly larger than
#' the number of clusters you expect to get.
#' @param min_bin_size Bins with fewer events than this threshold are considered outliers, and are 
#' ignored during the weighted EM algorithm. These events can still be assigned to clusters during
#' the prediction phase.
#' @param max_bin_size Bins with more events than this threshold are split, to ensure that the
#' weighted EM algorithm closely approximates a run of vanilla EM on the entire dataset.
#' @param mixtures_1D_path A path for pre-computed 1D mixture models, to be used for binning.
#' These are computed from scratch if the path is not provided.
#' @param seed A seed for the random number generator, to be used in the initialization of the bin
#' centers.
#' @param do_tsne Boolean flag; if true, computes t-SNE reduction to 2D of the bin centers,
#' for use in visualizations.
#' @param do_variance_correction Boolean flag; if true, corrects the variances of the mixture
#' components found in weighted EM, to account for intra-bin variance. TRUE is strongly recommended.
#' @param parallel Boolean flag; if true, uses multithreading to speed up computation.
#' @param verbose Boolean flag; if true, outputs timing and milestone information.
#' @return A tailor object containing:
#' \describe{
#'   \item{fit}{The tailor model, a named list containing the mixture proportions, means and variances
#'   of all mixture components.}
#'   \item{mixtures_1D}{}
#'   \item{cat_clusters}{A named list containing information about the categorical clusters found by
#'   the model: phenotype, cluster centers, and a mapping from mixture components to categorical clusters.}
#'   \item{tsne_centers}{Optional: a dimensional reduction to 2D of the bin centers.}
#' }
#' @examples Write this!
#' @export
tailor.learn = function(data, params = NULL, 
                        mixture_components = 200,
                        min_bin_size = 20, max_bin_size = 200, 
                        mixtures_1D_path = NULL, 
                        seed = 135, 
                        do_tsne = TRUE,
                        do_variance_correction = TRUE,
                        parallel = FALSE,
                        verbose = TRUE) 
{ 
  # TO DO: add support for FlowFrames, FlowSets

  if (is.null(params)) params = colnames(data)
  d = length(params) 
  
  set.seed(seed)
  
  # preliminary binning 
    
  if (is.null(mixtures_1D_path)) 
  { 
    # Learn 1D mixture model for each variable 
    if (verbose) {print("Getting 1D mixtures...")} 
    mixtures_1D = list()
    
    sample_1D = floor(nrow(data)/5) # TO DO: interpolate between n/5 and n/10 based on dataset size
    
    mixtures_1D$mixtures = get_1D_mixtures(data[,params], params, max_mixture = 3,  
                                    sample_size = sample_1D, seed = seed,
                                    parallel = parallel,
                                    verbose = verbose) 
    
    # Merge modes whose mean is small enough. 
    # These are likely compensation artifacts. 
    mixtures_1D$to_merge = find_to_merge(mixtures_1D, params,  
                             negative_threshold = 0.5, 
                             verbose = FALSE) 
  } else 
  { 
    # If given path to pre-computed mixtures, load them 
    load(mixtures_1D_path) 
    tailor_params = names(mixtures_1D$mixtures)
  } 
    
    
  # Up-sample: assign each data point to most likely component of the model 
  if (verbose) { print("Assigning data to 1D mixtures...")} 
    
  mapping = mapping_from_mixtures(data[,params], mixtures_1D$mixtures, mixtures_1D$to_merge, 
                                  params, 
                                  parallel = parallel, verbose = verbose) 
  phenobin = phenobin_label(mapping) 
 
  
  if (verbose) { print("Preparing initialization of bulk mixture model...")} 
  
  # Parse phenobins 
  phenobin_summary = get_phenobin_summary(phenobin) 
  
  large_bins = which(phenobin_summary$bins_sorted >= min_bin_size) 
  bins = as.integer(names(phenobin_summary$bins_sorted)[large_bins]) 
  sizes = as.vector(phenobin_summary$bins_sorted)[large_bins] 
  populous = which(phenobin_summary$predictions %in% bins) 
  
  if(verbose) 
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
                                           verbose = FALSE) 
  repr_means = phenobin_parameters$means 
  repr_variances = phenobin_parameters$variances 
  repr_sizes = phenobin_parameters$sizes 
  
  
  if (verbose) 
  { 
    cat("Total bin number after selecting and splitting: ", nrow(repr_means), "\n") 
    cat("Maximum bin size: ", max(repr_sizes), "\n") 
  } 
  
  
  # Prepare initialization of bulk mixture model 
  init_bins = c(1:mixture_components) 
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
                                           verbose = FALSE)   
  
  
  mixture = list() 
  mixture$pro = sizes[init_bins] 
  mixture$mean = phenobin_parameters$means 
  mixture$variance$sigma = array(NaN, c(d,d,mixture_components)) 
  
  for (component in init_bins) 
  { 
    sigma = phenobin_parameters$variances[component,,]
    mixture$variance$sigma[,,component] = sigma
  } 
  
  
  # Learn bulk mixture model on weighted subsample 
  if (verbose) { print("Running bulk mixture model...")} 
  
  if (do_variance_correction)
  {
    fit = bulk_weighted_gmm(data = repr_means, 
                            k = mixture_components, 
                            params = params,  
                            weights = repr_sizes,  
                            initialize = mixture, 
                            regularize_variance = TRUE, 
                            variance_correction = repr_variances, 
                            verbose = TRUE) 
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
                            verbose = TRUE) 
  }
  
  # Use the 1D mixture models to decide on +/- cutoffs
  cutoffs = get_1D_cutoffs(mixtures_1D, to_merge, tailor_params)
  
  # Categorical merging
  cat_clusters = categorical_merging(fit$mixture$pro, 
                             fit$mixture$mean,
                             cutoffs, 
                             tailor_params)

  if (do_tsne)
  {
    if (verbose) { print("Reducing phenobin centers to 2D...") }
    tsne_centers = get_tsne_centers(data = repr_means, seed = seed, probs = fit$event_probabilities)
  }


  if (do_tsne)
  {
    tailor_obj = list("fit" = fit, "mixtures_1D" = mixtures_1D, "cat_clusters" = cat_clusters, 
         "tsne_centers" = tsne_centers) 
  }
  else
  {
    tailor_obj = list("fit" = fit, "mixtures_1D" = mixtures_1D, "cat_clusters" = cat_clusters) 
  }
  
  class(tailor_obj) = "tailor"
  tailor_obj
} 


################################################################################
# Given a computed mixture model, assign new data to components
# Input:
#     data: a matrix containing events along rows, variables along columns
#     mixture: a mixture model, e.g. the output of bulk_weighted_gmm; list with 3 components:
#             pro: a vector of length = k, the proportions of the mixture
#             mean: a matrix with dim = (k, length(params));
#                   each row is the mean of one mixture component
#             variance$sigma: an array with dim = (length(params), length(params),k)
#                   each slice is the variance of one mixture component
#     params: the variables to use during modeling; must be subset of colnames(data)
#     n_batch: integer; partition the data in batches and process separately,
#             to reduce memory usage
#     verbose: boolean flag; if true, print info about running time and current batch
# Output:
#     mapping: a vector of length = nrow(data); mapping[i] is the most likely mixture
#             component to draw event i from
# Remark about output:
#     It may be desirable to output the probabilities of assigning each event to each
#     component, rather than a hard assignment. However, this requires n*k floats,
#     which may be larger than the available memory.
################################################################################

#' @title tailor.predict
#' @description Takes as input a tailor object and some data (could be the data used
#' to learn the tailor object, or some new data). Computes, for each event, the mixture
#' component from which it is most likely drawn, then maps this mixture component to its
#' corresponding categorical cluster.
#' @param data A matrix containing events along the rows, markers along columns.
#' @param tailor_obj A tailor object containing information about mixture components
#' and categorical clusters. Can be obtained as the output of tailor.learn.
#' @param n_batch A naive implementation would need nrow(data)*mixture_components memory.
#' To reduce memory usage, process data in batches. 
#' @param parallel Boolean flag; if true, uses multithreading to process batches in parallel.
#' For optimal runtime, if parallel = TRUE, n_batch should be a multiple of the number of
#' cores available, as returned by parallel::detectCores().
#' @param verbose Boolean flag; if true, outputs timing and milestone information.
#' @return An atomic vector of integers, giving the categorical cluster for each event.
#' @export
tailor.predict = function(data, tailor_obj, n_batch = 64,
                          parallel = TRUE, verbose = FALSE)
{
  start.time = Sys.time()
  n = nrow(data)
  k = length(tailor_obj$mixture$pro)
  logpro = log(tailor_obj$mixture$pro)
  params = colnames(tailor_obj$fit$mixture$mean)
  
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
          weight = tailor_obj$mixture$pro[cl]
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
  
  tailor_obj$cat_clusters$mixture_to_cluster[mapping]
}



#' @title get_1D_mixtures
#' @description Computes 1D mixture model for each marker separately, for use in binning step
#' of tailor. It is difficult to find settings which work for all datasets. Therefore, it is
#' recommended to inspect the results with inspect_1D_mixtures, and run get_1D_mixtures_custom
#' for problematic markers.
#' @param data A matrix containing events along the rows, markers along columns.
#' @param params A list of markers to use; must be subset of colnames(data).
#' @param max_mixture Will attempt to model each marker as k mixture components, for
#' 1 <= k <= max_mixture. The best k is chosen based on a modified version of the Bayesian
#' Information Criterion (BIC).
#' @param prior_BIC Make this larger to favor a smaller number of mixture components.
#' @param sample_size An integer, at most nrow(data). Used to subsample the data in the
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
                           prior_BIC = 600, sample_size = NULL, seed = 135, 
                           parallel = FALSE,
                           verbose = FALSE) 
{ 
  start.time = Sys.time()
  
  data = data[,params]
  
  # Keep all data, or sample a subset to speed up 
  set.seed(seed) 
  if (is.null(sample_size)) 
  { 
    sel = c(1:nrow(data)) 
  } else 
  { 
    sel = sample(nrow(data), sample_size) 
  } 
  
  
  if (parallel)
  {
    #setup parallel backend to use many processors
    cores=detectCores()
    cl <- makeCluster(cores[1])
    # cl <- makeCluster(1) # For debugging only: similar to unparallelized version, with 10% overhead
    registerDoParallel(cl)
    
    if (verbose) cat("Learning", length(params), "1D mixtures in parallel...")
    
    fit_list = foreach (data_param = iter(data[sel,], by = 'col'), .packages = c("mclust")) %dopar%
      { 
        param = colnames(data_param)[1]
        
        set.seed(seed)
        
        # Use Bayesian information criterion to choose best k 
        BIC = mclustBIC(data_param, G = 1:max_mixture, model = "V", verbose = FALSE) 
        
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
      BIC = mclustBIC(dat, G = 1:max_mixture, model = "V", verbose = FALSE) 
      
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




get_1D_mixtures_custom = function(data, params, k = 3,  
                                  sample_size = NULL, seed = 135, 
                                  separate_neg = FALSE, 
                                  verbose = FALSE) 
{ 
  start.time = Sys.time() 
  fit_list = list() 
  
  set.seed(seed) 
  # Keep all data, or sample a subset to speed up 
  if (is.null(sample_size)) 
  { 
    sel = c(1:nrow(data)) 
  } else 
  { 
    sel = sample(nrow(data), sample_size) 
  } 
  
  for (param in params) 
  { 
    if (verbose) {cat(param, " ")} 
    
    dat = data[sel,param] 
    
    # If flag is true, get rid of negative values, which are compensation artifacts 
    if (separate_neg) 
    { 
      negatives = dat[which(dat < 0)] 
      dat = dat[which(dat > 0)] 
    } 
    
    # Fit the model with the chosen k 
    fit = Mclust(dat, G = k, modelNames = "V", verbose = FALSE) 
    fit_list[[param]] = fit$parameters 
    
    # If flag is true, model negative population separately, and add to mixture list 
    if (separate_neg) 
    { 
      fit = Mclust(negatives, G = 1, modelNames = "V", verbose = FALSE) 
      fit_list[[param]]$pro = c(fit$parameters$pro, fit_list[[param]]$pro) 
      fit_list[[param]]$mean = c(fit$parameters$mean, fit_list[[param]]$mean) 
      fit_list[[param]]$variance$sigmasq = c(fit$parameters$variance$sigmasq,  
                                             fit_list[[param]]$variance$sigmasq) 
      fit_list[[param]]$variance$scale = c(fit$parameters$variance$scale,  
                                           fit_list[[param]]$variance$scale) 
      
      # Rescale mixture proportions 
      ln = length(negatives) 
      lp = length(dat) 
      
      fit_list[[param]]$pro[1] = fit_list[[param]]$pro[1] * ln / (ln + lp) 
      fit_list[[param]]$pro[c(2:(k+1))] = fit_list[[param]]$pro[c(2:(k+1))] * lp / (ln + lp) 
    } 
  } 
  
  end.time = Sys.time() 
  if(verbose) {cat("\n")} 
  if(verbose) {print(end.time - start.time)} 
  
  fit_list 
} 



#' @title inspect_1D_mixtures
#' @description Plot the result of 1D mixture model calculation for visual inspection.
#' Displays, for each marker, three side-by-side plots, giving a kernel density estimate
#' for the data and that marker, the Gaussian mixture, and the separate mixture components,
#' respectively.
#' @param data A matrix containing events along the rows, markers along columns.
#' @param mixtures_1D 1D mixture models, as produced by get_1D_mixtures.
#' @param params A list of markers to use; must be subset of colnames(data).
#' @export
inspect_1D_mixtures = function(data, mixtures_1D, params)
{
  global_kdes = make_kdes_global(data, params)
  
  for (param in params)
  {
    plot_kde_vs_mixture(data, global_kdes, mixtures = mixtures_1D$mixtures, name = param)
    dev.print(png, paste("mixture_", param, ".png", sep = ""), width = 600, height = 400)
  }
}


#' @title categorical_labelling
#' @description 
#' @param categorical_clusters Information about categorical clusters obtained from a tailor
#' object.
#' @param defs A matrix or data frame giving definitions of major phenotypes (e.g. "CD4 Naive").
#' @param params A list of markers to use; must be subset of colnames(data).
#' @return WRITE THIS.
#' @export
categorical_labeling = function(categorical_clusters, defs, params)
{
  #################################################
  # Take input the list of categorical clusters
  # and definitions of major phenotypes (e.g. naive
  # CD4 cells). Label each cluster by one major
  # phenotype, or unknown.
  #################################################
  
  n = length(cats$categs)
  cats[["labels"]] = vector(mode = "character", length = n)
  labs = rownames(defs)
  nam = names(defs)
  ind = vector(mode = "integer", length = length(nam))
  for (i in seq(length(nam)))
  {
    ind[i] = which(params == nam[i])
  }
  
  for (i in seq(n))
  {
    cats$labels[i] = "UNK"
    
    for (j in seq(nrow(defs)))
    {
      match = TRUE
      for (k in seq(ncol(defs)))
      {
        if (defs[j,k] == "hi" & substr(cats$categs[i],ind[k],ind[k]) == "-" | 
            defs[j,k] == "lo" & substr(cats$categs[i],ind[k],ind[k]) == "+")
        {
          match = FALSE
        }
      }
      if (match)
      {
        cats$labels[i] = labs[j]
        break
      }
    }
  }
  
  cats
}


