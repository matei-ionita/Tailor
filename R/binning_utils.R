################################################################################
################################################################################
#
# Wrappers and utilities for 1D mixture modelling in flow cytometry data
#
################################################################################
################################################################################





get_weighted_subsample <- function(data, phenobin_summary, params,
                                   min_bin_size, max_bin_size, verbose)
{
  large_bins <- which(phenobin_summary$bins_sorted >= min_bin_size)
  bins <- as.integer(names(phenobin_summary$bins_sorted)[large_bins])
  sizes <- as.vector(phenobin_summary$bins_sorted)[large_bins]
  populous <- which(phenobin_summary$predictions %in% bins)

  weighted_subsample <- find_phenobin_mean(data = data,
                                            predictions = phenobin_summary$predictions,
                                            bins = bins, sizes = sizes,
                                            params = params,
                                            selected = populous,
                                            split_threshold = max_bin_size,
                                            compute_var = TRUE,
                                            verbose = (verbose >= 1))

  return(weighted_subsample)
}



get_init <- function(data, phenobin_summary, params,
                     min_bin_size, mixture_components, verbose)
{
  large_bins <- which(phenobin_summary$bins_sorted >= min_bin_size)
  bins <- as.integer(names(phenobin_summary$bins_sorted)[large_bins])
  sizes <- as.vector(phenobin_summary$bins_sorted)[large_bins]

  candidate_bins <- seq_len(mixture_components)
  init_bins <- which(sizes[candidate_bins] > 3 * min_bin_size) # heuristic: strive to improve
  lost <- length(candidate_bins) - length(init_bins)
  if (lost > 0) cat("Warning: dropped", lost, "mixture components due to too few events at initialization.\n")

  populous <- which(phenobin_summary$predictions %in% bins[init_bins])
  init_parameters <- find_phenobin_mean(data = data,
                                        predictions = phenobin_summary$predictions,
                                        bins = bins[init_bins], sizes = sizes[init_bins],
                                        params = params,
                                        selected = populous,
                                        split_threshold = NULL,
                                        compute_var = TRUE,
                                        verbose = (verbose >= 1))


  mixture <- list()
  mixture$pro <- sizes[init_bins]
  mixture$mean <- init_parameters$means
  mixture$variance$sigma <- array(NaN, c(d,d,mixture_components))

  for (component in init_bins) {
    sigma <- init_parameters$variances[component,,]
    mixture$variance$sigma[,,component] <- sigma
  }

  return(mixture)
}

################################################################################
# Given a partition of the data into "phenobins", i.e. bins determined by 1D
# mixture modelling across each individual variable, find the mean and variance
# of each phenobin.
#
# Input:
#       data: a matrix with events along rows and variables along columns
#       predictions: a vector of length = nrow(data), phenobin assignment for each event
#       bins: those bins for which information should be computed (computing variances
#             for more than 100k or so bins may cause memory issues)
#       sizes: a vector of length = length(bins), counting the events in each box
#       params: parameters to use
#       compute_var: boolean flag; if true, computes variance in addition to mean
#       verbose: boolean flag; if true, prints information at checkpoints
# Output: list with two elements:
#       means: matrix of dim = (length(bins), length(params));
#             rows are means of events in a box
#       variance: matrix of dim = (length(bins), length(params), length(params));
#             slices are variances of events in a box
################################################################################

find_phenobin_mean <- function(data, predictions,
                              bins, sizes,
                              params, selected,
                              split_threshold = NULL,
                              compute_var = FALSE,
                              parallel = FALSE,
                              verbose = FALSE)
{
  n <- nrow(data)
  k <- length(bins)
  d <- length(params)

  if(is.null(split_threshold)) {
    n_repr <- rep(1, k)
    split_box_sizes <- sizes
  } else {
    n_repr <- ceiling(sizes/split_threshold)
    split_box_sizes <- integer(sum(n_repr))
  }

  means <- matrix(nrow = sum(n_repr), ncol = d)
  colnames(means) <- params

  if(compute_var) {
    variances <- array(0, c(sum(n_repr),d,d))
  }

  idx <- vector("list", k)
  counts <- vector("list", k)

  # Initialize lists of indices with appropriate size
  for (i in seq_len(k)) {
    box <- bins[i]
    idx[[box]] <- integer(length = sizes[i])
    counts[[box]] <- 0
  }

  if(verbose) { print("Done init!")}

  # For each event, place its index in the appropriate index list
  for (row in selected) {
    box <- predictions[row]
    counts[[box]] <- counts[[box]] + 1
    idx[[box]][counts[[box]]] <- row
  }

  if(verbose) { print("Compiled index lists.")}

  start <- 1
  for (i in seq_len(k)) {
    # Split each box until it has fewer than split_threshold events; compute their means
    box <- bins[i]
    sel <- idx[[box]]

    if (n_repr[i] == 1) {
      means[start,] <- apply(data[sel,params,drop=FALSE], 2, mean)
    } else {
      if (length(sel) > 100000) {
        # kmeans is just for binning purposes, don't care about convergence.
        # If there are many events, cap at 3 iterations to save time.
        km <- suppressWarnings(kmeans(x = round(data[sel,params,drop = FALSE],3),
                    centers = n_repr[i], iter.max = 3))
      } else {
        km <- suppressWarnings(kmeans(x = round(data[sel,params,drop = FALSE],3),
                    centers = n_repr[i]))
      }
      means[c(start:(start + n_repr[i]-1)),] <- km$centers
    }

    if (compute_var) {
      # Compute the in-box variance

      if (length(sel) > 1) {
        # If box wasn't split, compute variance for entire box
        if (n_repr[i] == 1) {
          split_box_sizes[start] <- sizes[i]
          var <- var(data[sel,params,drop = FALSE])
          variances[start,,] <- var
        } else {
          # Otherwise, compute for each piece
          data_this_box <- data[sel,params]

          for (cl in seq_len(n_repr[i])) {
            ind_cl <- which(km$cluster ==cl)
            split_box_sizes[start + cl -1] <- length(ind_cl)

            var <- var(data_this_box[ind_cl,,drop = FALSE])
            if (sum(is.na(var)) != 0) {
              var <- 0
            }

            variances[start+cl-1,,] <- var
          }
        }

      } else {
        var <- 0
        variances[start,,] <- var
      }

    }
    start <- start + n_repr[i]
  }

  if (compute_var) {
    l <- list("means" = means, "variances" = variances,
              "n_repr" = n_repr, "sizes" = split_box_sizes)
  } else {
    l <- list("means" = means, "n_repr" = n_repr, "sizes" = split_box_sizes)
  }
  return(l)
}




# Find pairs of mixture components that represent the same biological population
# So far, just merging components with mean < 0.5; any differences below this value
# are likely artifacts of compensation.
find_to_merge <- function(mixtures, params, negative_threshold = 0.5, verbose = FALSE)
{
  d <- length(mixtures)
  to_merge <- list()

  for (param in params) {
    mixture <- mixtures[[param]]
    k <- length(mixture$pro)

    small <- which(mixture$mean < negative_threshold)

    if (length(small) >1) {
      to_merge[[param]] <- as.vector(small)
    }
  }

  return(to_merge)
}


# Map each data point to the most probable mixture component for all parameters.
# to_merge is a list specifying components that should be merged; probabilities
# are added in this case.
mapping_from_mixtures <- function(data, mixtures, to_merge, params,
                                 parallel = FALSE, verbose = FALSE)
{
  d <- length(mixtures)
  n <- nrow(data)

  col_list <- NULL

  if (parallel) {
    if (verbose) { cat("Mapping", n, "events to", d, "univariate mixtures in parallel...") }
    cl <- start_parallel_cluster()

    data_param <- NULL # unnecessary definition, but R CMD check complains without it

    mapping <- foreach(data_param = iter(data, by='col'),
                      .combine = cbind, .packages = c("mvtnorm")) %dopar%
    {
      param <- colnames(data_param)[1]
      mixture <- mixtures[[param]]
      k <- length(mixture$pro)
      n_small <- length(to_merge[[param]])

      # If there is only one component left after merging, skip this parameter
      if (k < 2 | k == n_small) {
        rep(0L,n)
      } else {
        to_sum <- NULL
        if (n_small > 1) {
          sum_base <- to_merge[[param]][1]
          to_sum <- to_merge[[param]][c(2:n_small)]
        }

        # Compute posterior probabilities for each class (up to a normalization factor)
        posteriors <- matrix(0, nrow = n, ncol = k - length(to_sum))
        summed <- 0

        for (i in seq(k)) {
          # If class i is to be merged with class sum_base, add the probabilities
          if (i %in% to_sum) {
            summed <- summed + 1
            posteriors[,sum_base] <- posteriors[,sum_base] + mixture$pro[i] * dnorm(data_param, mean = mixture$mean[i],
                                                                                   sd = sqrt(mixture$variance$sigmasq[i]))
          } else {
            posteriors[,i-summed] <- mixture$pro[i] * dnorm(data_param, mean = mixture$mean[i],
                                                           sd = sqrt(mixture$variance$sigmasq[i]))
          }
        }

        # The class assignment is given by the maximum posterior probability
        apply(posteriors, 1, which.max)
      }
    }

    stopCluster(cl)

    colnames(mapping) <- params
    for (param in params) {
      if (mapping[1,param] != 0) {col_list = c(col_list, param)}
    }
    mapping <- mapping[,col_list]
  } else {
    mapping <- matrix(nrow = n, ncol = 0)

    if (verbose) { cat("Mapping", n, "events to", d, "univariate mixtures sequentially: ") }

    for (ind in seq(d)) {
      param <- params[ind]
      mixture <- mixtures[[param]]
      k <- length(mixture$pro)
      n_small <- length(to_merge[[param]])

      if (verbose) { cat(param, " ") }

      # If there is only one component left after merging, skip this parameter
      if (k < 2 | k == n_small) { next }

      to_sum <- NULL
      if (n_small > 1) {
        sum_base <- to_merge[[param]][1]
        to_sum <- to_merge[[param]][c(2:n_small)]
      }

      # Compute posterior probabilities for each class (up to a normalization factor)
      posteriors <- matrix(0, nrow = n, ncol = k - length(to_sum))
      summed <- 0

      for (i in seq(k)) {
        # If class i is to be merged with class sum_base, add the probabilities
        if (i %in% to_sum) {
          summed <- summed + 1
          posteriors[,sum_base] <- posteriors[,sum_base] + mixture$pro[i] * dnorm(data[,param], mean = mixture$mean[i],
                                                                                 sd = sqrt(mixture$variance$sigmasq[i]))
        } else {
          posteriors[,i-summed] <- mixture$pro[i] * dnorm(data[,param], mean = mixture$mean[i],
                                                         sd = sqrt(mixture$variance$sigmasq[i]))
        }
      }
      # The class assignment is given by the maximum posterior probability
      class <- apply(posteriors, 1, which.max)

      # Put together the class assignments for all parameters
      mapping <- cbind(mapping, class, deparse.level = 0)
      col_list <- c(col_list, param)

      rm(posteriors, class)
      gc()
    }
    colnames(mapping) <- col_list
  }
  if(verbose) {cat("\n")}

  return(mapping)
}


# Collapse all class assignments into a string; this is the phenobin label
phenobin_label <- function(mapping)
{
  apply(mapping, 1, function(x) {paste(x, collapse = "")})
}


# Process phenobin assignment, and return list of bins sorted by size
get_phenobin_summary <- function(phenobin)
{
  t <- table(phenobin)
  t_sorted <- rev(order(t))
  box_idx <- t
  box_idx[names(box_idx)] <- seq_len(length(box_idx))
  predictions <- as.vector(box_idx[phenobin])
  t_numeric <- table(predictions)
  t_numeric_sorted <- rev(order(t_numeric))

  l <- list("predictions" = predictions, "bins_sorted" = t_numeric[t_numeric_sorted])
  return(l)
}


get_1D_mixtures_default <- function(data, params, parallel, verbose)
{
  if (verbose > 0) {print("Getting 1D mixtures...")}

  sample_fraction <- 0.5
  if (nrow(data) > 5e5) sample_fraction <- 0.2
  if (nrow(data) > 1e7) sample_fraction <- 0.1

  mixtures_1D <- get_1D_mixtures(data[,params], params, max_mixture = 3,
                                 sample_fraction = sample_fraction,
                                 parallel = parallel,
                                 verbose = (verbose >= 1))

  # Merge modes whose mean is small enough.
  # These are likely compensation artifacts.
  mixtures_1D$to_merge <- find_to_merge(mixtures_1D, params,
                                        negative_threshold = 0.5,
                                        verbose = (verbose == 1))
}


get_1D_mixtures_custom <- function(data, params, k = 3,
                                  sample_fraction = 0.2,
                                  separate_neg = FALSE,
                                  verbose = FALSE)
{
  fit_list <- list()

  # Keep all data, or sample a subset to speed up
  if (sample_fraction == 1) {
    sel <- seq_len(nrow(data))
  } else
  {
    sample_size <- ceiling(sample_fraction * nrow(data))
    sel <- sample(nrow(data), sample_size)
  }

  for (param in params) {
    if (verbose) {cat(param, " ")}

    dat <- data[sel,param]
    # If flag is true, get rid of negative values, which are compensation artifacts
    if (separate_neg) {
      negatives <- dat[which(dat < 0)]
      dat <- dat[which(dat > 0)]
    }

    # Fit the model with the chosen k
    fit <- Mclust(dat, G = k, modelNames = "V", verbose = FALSE)
    fit_list[[param]] <- fit$parameters

    # If flag is true, model negative population separately, and add to mixture list
    if (separate_neg) {
      fit <- Mclust(negatives, G = 1, modelNames = "V", verbose = FALSE)
      fit_list[[param]]$pro <- c(fit$parameters$pro, fit_list[[param]]$pro)
      fit_list[[param]]$mean <- c(fit$parameters$mean, fit_list[[param]]$mean)
      fit_list[[param]]$variance$sigmasq <- c(fit$parameters$variance$sigmasq,
                                             fit_list[[param]]$variance$sigmasq)
      fit_list[[param]]$variance$scale <- c(fit$parameters$variance$scale,
                                           fit_list[[param]]$variance$scale)

      # Rescale mixture proportions
      ln <- length(negatives)
      lp <- length(dat)

      fit_list[[param]]$pro[1] <- fit_list[[param]]$pro[1] * ln / (ln + lp)
      fit_list[[param]]$pro[c(2:(k+1))] <- fit_list[[param]]$pro[c(2:(k+1))] * lp / (ln + lp)
    }
  }
  if(verbose) {cat("\n")}

  return(fit_list)
}




#########################
# Some plotting utilities
#########################

plot_distribution_1d <- function(dat, mixtures, name,
                                min = -2, max = 4,
                                separate = FALSE)
{
  p <- mixtures[[name]]
  pts <- seq(from = min(dat[,name]), to = max(dat[,name]), length.out = 500)
  k <- length(p$pro)

  # Vector of colors - assuming not more than 4 components
  colors <- c("blue", "red", "green", "orange")

  # Select the first mixture component
  m <- p$mean[1]
  sd <- sqrt(p$variance$sigmasq[1])
  g <- p$pro[1] * dnorm(pts, mean = m, sd = sd)

  # Plot the first mixture component
  if (separate) {
    plot(pts, g, type = "l", xlim = c(min,max), ylim = c(0,1),
         xlab = name, ylab = "density",
         col = "blue", main = "components", cex.main = 2)
  }

  # Same for other components
  if (k > 1) {
    for (comp in c(2:k)) {
      m <- p$mean[comp]
      sd <- sqrt(p$variance$sigmasq[comp])

      g1 <-  p$pro[comp] * dnorm(pts, mean = m, sd = sd)

      if (separate) {
        lines(pts, g1, col = colors[comp])
      } else
      {
        g <- g + g1
      }
    }
  }

  if (!separate) {
    plot(pts, g, type = "l", xlim = c(min,max), ylim = c(0,1),
         xlab = name, ylab = "density", main = "Gaussian mixture", cex.main = 2)
  }
}








