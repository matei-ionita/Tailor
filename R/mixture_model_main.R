######################################################################################## 
# 
# Main wrapper for mixture modeling workflow 
# 
# Input: 
#       data: a matrix with data points along rows, variables along columns 
#       params: a list of variables to use; must be subset of colnames(data) 
#       k_bulk: number of components in the d-dimensional (bulk) mixture model 
#       n_bulk: maximum number of representatives to learn bulk model on 
#       max_steps: maximum number of weighted EM iterations to perform 
#       min_box_size: choose representative only for boxes larger than this threshold  
#       max_box_size: choose multiple representatives (split with k-means) for boxes 
#                 larger than this threshold 
#       phenobox_path: path to a list assigning each data point to a phenobox; 
#                 if provided, 1D modeling is skipped 
#       mixtures_1D_path: path to a list of 1D "mixture" objects;  
#                 if provided, 1D modeling is skipped 
#       sample_1D: number of data points to sample for 1D modeling 
#       seed: a random seed for the above; the model is v. robust w.r.t. seed 
#       kmax_1D: maximum number of mixture components for each 1D model; 
#                 BIC is used to determine 1 <= k <= kmax_1D 
#       merge_1d_negative_threshold: numeric threshold for merging 1D modes having 
#                 small mean; these are usually compensation artifacts 
#       to_merge: a list of other, user-provided, modes to merge 
#       verbose: boolean flag; if true, outputs timing and milestone information 
# Output: 
#       fit: a "mixture" object describing the d-dimensional mixture model 
#       assignment_pred: an array assigning each data point to the most likely  
#                 mixture component 
#       mixtures_1D: a list of 1D "mixture" objects; 
#       phenobox_summary: a summary of the phenobox information 
######################################################################################## 

learn_mixture_model = function(data, params, 
                               k_bulk = 200, n_bulk = 10000, max_steps = 30,
                               min_box_size = 75, max_box_size = 1000, 
                               phenobox_path = NULL, 
                               mixtures_1D_path = NULL, 
                               sample_1D = NULL, seed = 135, 
                               kmax_1D = 3, merge_1D_negative_threshold = 0.5, 
                               to_merge = NULL, 
                               n_batch = 64,
                               do_tsne = FALSE,
                               variance_correction = FALSE,
                               parallel = FALSE,
                               verbose = TRUE) 
{ 
  start.time = Sys.time() 
  d = length(params) 
  set.seed(seed)
  
  # Compute phenobox assignment. 
  if (is.null(phenobox_path)) 
  { 
    # Learn 1D mixture model for each variable 
    if (verbose) {print("Getting 1D mixtures...")} 
    
    if (is.null(mixtures_1D_path)) 
    { 
      if (is.null(sample_1D)) { sample_1D = floor(nrow(data)/10) } 
      mixtures_1D = get_1D_mixtures(data, params, max_mixture = kmax_1D,  
                                    sample_size = sample_1D, seed = seed,
                                    parallel = parallel,
                                    verbose = verbose) 
    } else 
    { 
      # If given path to pre-computed mixtures, load them 
      load(mixtures_1D_path) 
    } 
    
    
    # If merge list not already provided, merge modes whose mean is small enough. 
    # These are likely compensation artifacts. 
    if (is.null(to_merge)) 
    { 
      to_merge = find_to_merge(mixtures_1D, params,  
                               negative_threshold = merge_1D_negative_threshold, 
                               verbose = FALSE) 
    } 
    
    
    # Up-sample: assign each data point to most likely component of the model 
    if (verbose) { print("Assigning data to 1D mixtures...")} 
    
    mapping = mapping_from_mixtures(data, mixtures_1D, to_merge, params, 
                                    parallel = parallel, verbose = verbose) 
    phenobox = phenobox_label(mapping) 
    
  } else 
  { 
    # If given path to pre-computed assignment, load it 
    load(phenobox_path) 
  } 
  
  if (verbose) { print("Preparing initialization of bulk mixture model...")} 
  
  # Parse phenoboxes 
  phenobox_summary = get_phenobox_summary(phenobox) 
  
  large_boxes = which(phenobox_summary$boxes_sorted >= min_box_size) 
  boxes = as.integer(names(phenobox_summary$boxes_sorted)[large_boxes]) 
  sizes = as.vector(phenobox_summary$boxes_sorted)[large_boxes] 
  populous = which(phenobox_summary$predictions %in% boxes) 
  
  if(verbose) 
  { 
    cat("Total box number before splitting: ", length(phenobox_summary$boxes_sorted), "\n") 
    cat("Maximum box size: ", sizes[1], "\n") 
  } 
  
  phenobox_parameters = find_phenobox_mean(data = data,  
                                           predictions = phenobox_summary$predictions,  
                                           boxes = boxes, sizes = sizes,  
                                           params = params, 
                                           selected = populous, 
                                           split_threshold = max_box_size, 
                                           compute_var = TRUE,
                                           seed = seed,
                                           parallel = parallel,
                                           verbose = FALSE) 
  repr_means = phenobox_parameters$means 
  repr_variances = phenobox_parameters$variances 
  repr_sizes = phenobox_parameters$sizes 
  
  
  if (verbose) 
  { 
    cat("Total box number after selecting and splitting: ", nrow(repr_means), "\n") 
    cat("Maximum box size: ", max(repr_sizes), "\n") 
  } 
  
  
  # Prepare initialization of bulk mixture model 
  
  init_boxes = c(1:k_bulk) 
  populous = which(phenobox_summary$predictions %in% boxes[init_boxes]) 
  phenobox_parameters = find_phenobox_mean(data = data,  
                                           predictions = phenobox_summary$predictions,  
                                           boxes = boxes[init_boxes], sizes = sizes[init_boxes],  
                                           params = params, 
                                           selected = populous, 
                                           split_threshold = NULL, 
                                           compute_var = TRUE,
                                           seed = seed,
                                           parallel = parallel,
                                           verbose = FALSE)   
  
  
  mixture = list() 
  mixture$pro = sizes[init_boxes] 
  mixture$mean = phenobox_parameters$means 
  mixture$variance$sigma = array(NaN, c(d,d,k_bulk)) 
  
  for (component in init_boxes) 
  { 
    sigma = phenobox_parameters$variances[component,,]
    mixture$variance$sigma[,,component] = sigma
  } 
  
  
  # Learn bulk mixture model on weighted subsample 
  if (verbose) { print("Running bulk mixture model...")} 
  
  if (variance_correction)
  {
    fit = bulk_weighted_gmm(data = repr_means, k = k_bulk, params = params,  
                            weights = repr_sizes, steps = max_steps,  
                            initialize = mixture, regularize_variance = TRUE, 
                            variance_correction = repr_variances, 
                            verbose = TRUE) 
  }
  else
  {
    fit = bulk_weighted_gmm(data = repr_means, k = k_bulk, params = params,  
                            weights = repr_sizes, steps = max_steps,  
                            initialize = mixture, regularize_variance = TRUE, 
                            variance_correction = NULL, 
                            verbose = TRUE) 
  }
  
  # Use the learned model to predict most likely cluster for each datapoint 
  if (verbose) { print("Assigning data to bulk mixtures...")} 
  assignment_pred = predict_bulk_gmm(data = data, mixture = fit$mixture,  
                                     params = params, n_batch = n_batch, 
                                     parallel = parallel, verbose = TRUE) 
  
  if (do_tsne)
  {
    if (verbose) { print("Reducing phenobox centers to 2D...") }
    tsne_centers = get_tsne_centers(data = repr_means, seed = seed, probs = fit$event_probabilities)
  }

  
  if (verbose) { print(Sys.time() - start.time) } 
  
  if (is.null(phenobox_path)) 
  { 
    if (do_tsne)
    {
      list("fit" = fit, "assignment_pred" = assignment_pred,  
           "mixtures_1D" = mixtures_1D, "phenobox" = phenobox, "tsne_centers" = tsne_centers) 
    }
    else
    {
      list("fit" = fit, "assignment_pred" = assignment_pred,  
           "mixtures_1D" = mixtures_1D, "phenobox" = phenobox) 
    }

  } else 
  { 
    if (do_tsne)
    {
      list("fit" = fit, "assignment_pred" = assignment_pred, "tsne_centers" = tsne_centers) 
    }
    else
    {
      list("fit" = fit, "assignment_pred" = assignment_pred) 
    }
  } 
  
} 
