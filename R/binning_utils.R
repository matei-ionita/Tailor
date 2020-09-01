#####################################
#####################################
#
# Wrappers and utilities for binning
#
#####################################
#####################################





get_weighted_subsample <- function(data, bin_summary, params,
                                   min_bin_size, max_bin_size,
                                   verbose)
{
  large_bins <- which(bin_summary$bins_sorted >= min_bin_size)
  bins <- as.integer(names(bin_summary$bins_sorted)[large_bins])
  sizes <- as.vector(bin_summary$bins_sorted)[large_bins]
  populous <- which(bin_summary$predictions %in% bins)

  weighted_subsample <- find_bin_mean(data = data,
                                      predictions = bin_summary$predictions,
                                      bins = bins, sizes = sizes,
                                      params = params,
                                      selected = populous,
                                      split_threshold = max_bin_size,
                                      verbose = (verbose >= 1))

  return(weighted_subsample)
}



get_init <- function(data, bin_summary, params,
                     min_bin_size, mixture_components,
                     verbose)
{
  large_bins <- which(bin_summary$bins_sorted >= min_bin_size)
  bins <- as.integer(names(bin_summary$bins_sorted)[large_bins])
  sizes <- as.vector(bin_summary$bins_sorted)[large_bins]

  candidate_bins <- seq_len(mixture_components)
  init_bins <- which(sizes[candidate_bins] > 3 * min_bin_size)
  lost <- length(candidate_bins) - length(init_bins)
  if (lost > 0) cat("Warning: dropped", lost,
                    "mixture components due to too few events.\n")

  populous <- which(bin_summary$predictions %in% bins[init_bins])
  init_parameters <- find_bin_mean(data = data,
                                  predictions = bin_summary$predictions,
                                  bins = bins[init_bins], sizes = sizes[init_bins],
                                  params = params,
                                  selected = populous,
                                  split_threshold = NULL,
                                  verbose = (verbose >= 1))


  d = length(params)
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



find_bin_mean <- function(data, predictions, bins, sizes, params,
                               selected, split_threshold = NULL,
                               parallel = FALSE, verbose = FALSE)
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
  variances <- array(0, c(sum(n_repr),d,d))

  idx <- get_bin_list(bins, sizes, selected, predictions, k)

  start <- 1
  for (i in seq_len(k)) {
    box <- bins[i]
    sel <- idx[[box]]

    if (n_repr[i] == 1) {
      means[start,] <- apply(data[sel,params,drop=FALSE], 2, mean)
    } else {
      km <- split_bin_kmeans(data, params, sel, n_repr[i])
      means[c(start:(start + n_repr[i]-1)),] <- km$centers
    }

    bin_info <- compute_bin_variance(data, params, sel, n_repr[i], start,
                                     km, sizes[i], variances, split_box_sizes)
    variances <- bin_info$variances
    split_box_sizes <- bin_info$sizes

    start <- start + n_repr[i]
  }

  l <- list("means" = means, "variances" = variances,
            "n_repr" = n_repr, "sizes" = split_box_sizes)
  return(l)
}


get_bin_list <- function(bins, sizes, selected, predictions, k)
{
  idx <- vector("list", k)
  counts <- vector("list", k)

  # Initialize lists of indices with appropriate size
  for (i in seq_len(k)) {
    box <- bins[i]
    idx[[box]] <- integer(length = sizes[i])
    counts[[box]] <- 0
  }

  # For each event, place its index in the appropriate index list
  for (row in selected) {
    box <- predictions[row]
    counts[[box]] <- counts[[box]] + 1
    idx[[box]][counts[[box]]] <- row
  }

  return(idx)
}


split_bin_kmeans <- function(data, params, sel, k)
{
  if (length(sel) > 100000) {
    # kmeans is just for binning purposes, don't care about convergence.
    # If there are many events, cap at 3 iterations to save time.
    km <- suppressWarnings(kmeans(x = round(data[sel,params,drop = FALSE],3),
                                  centers = k, iter.max = 3))
  }
  else {
    km <- suppressWarnings(kmeans(x = round(data[sel,params,drop = FALSE],3),
                                  centers = k))
  }

  return(km)
}


compute_bin_variance <- function(data, params, sel, n_repr, start, km,
                                 size, variances, split_box_sizes)
{
  if (length(sel) > 1) {

    # If box wasn't split, compute variance for entire box
    if (n_repr == 1) {
      split_box_sizes[start] <- size
      var <- var(data[sel,params,drop = FALSE])
      variances[start,,] <- var
    } else {
      # Otherwise, compute for each piece
      for (cl in seq_len(n_repr)) {
        ind_cl <- which(km$cluster ==cl)
        split_box_sizes[start + cl -1] <- length(ind_cl)

        var <- var(data[sel[ind_cl],params,drop = FALSE])
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

  l = list(variances = variances, sizes = split_box_sizes)
  return(l)
}



# Find pairs of mixture components that represent the
# same biological population. So far, just merging components
# with mean < 0.5; any differences below this value
# are likely artifacts of compensation.
find_to_merge <- function(mixtures, params, negative_threshold = 0.5,
                          verbose = FALSE)
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


map_events_to_bins <- function(data, cutoffs)
{
  params  <- names(cutoffs)
  cutoffs <- as.numeric(cutoffs)

  mapping <- t(t(data[,params]) - cutoffs) # inelegant, but faster than sweep
  mapping[mapping >  0] <- 2L
  mapping[mapping <= 0] <- 1L
  mode(mapping) <- "integer"

  bins <- apply(mapping, 1, function(x) {paste(x, collapse = "")})
  bin_summary <- get_bin_summary(bins)
  return(bin_summary)
}


# Process bin assignment, and return list of bins sorted by size
get_bin_summary <- function(bins)
{
  t <- table(bins)
  t_sorted <- rev(order(t))
  box_idx <- t
  box_idx[names(box_idx)] <- seq_len(length(box_idx))
  predictions <- as.vector(box_idx[bins])
  t_numeric <- table(predictions)
  t_numeric_sorted <- rev(order(t_numeric))

  l <- list("predictions" = predictions,
            "bins_sorted" = t_numeric[t_numeric_sorted])
  return(l)
}


get_1D_mixtures_default <- function(data, params, parallel, verbose)
{
  if (verbose > 0) {print("Getting 1D mixtures...")}

  sample_fraction <- 1
  if (nrow(data) > 1e4) sample_fraction <- 0.5
  if (nrow(data) > 1e5) sample_fraction <- 0.1
  if (nrow(data) > 5e5) sample_fraction <- 0.06
  if (nrow(data) > 1e6) sample_fraction <- 0.03

  mixtures_1D <- get_1D_mixtures(data[,params], params, max_mixture = 3,
                                 sample_fraction = sample_fraction,
                                 parallel = parallel,
                                 verbose = (verbose >= 1))

  # Merge modes whose mean is small enough.
  # These are likely compensation artifacts.
  mixtures_1D$to_merge <- find_to_merge(mixtures_1D, params,
                                        negative_threshold = 0.5,
                                        verbose = (verbose == 1))

  return(mixtures_1D)
}


get_1D_mixtures_custom <- function(data, params, k = 3, sample_fraction = 0.2,
                                  separate_neg = FALSE, verbose = FALSE)
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


learn_1D_mixtures_parallel <- function(data, max_mixture, prior_BIC)
{
  cl <- start_parallel_cluster()

  data_param <- NULL
  fit_list <- foreach (data_param = iter(data, by = 'col'),
                       .packages = c("mclust")) %dopar%
    {
      param <- colnames(data_param)[1]

      # Use Bayesian information criterion to choose best k
      BIC <- mclustBIC(data_param, G = seq_len(max_mixture),
                       modelNames = "V", verbose = FALSE)

      # A tweak to favor smaller k
      for (k in seq_len(max_mixture))
      {
        BIC[k] <- BIC[k] - prior_BIC*k*log(length(data_param), base = 2)
      }

      # Fit the model with the chosen k
      fit <- Mclust(data_param, x=BIC, verbose = FALSE)
      fit$parameters
    }
  names(fit_list) <- colnames(data)

  #stop cluster
  stopCluster(cl)

  return(fit_list)
}


learn_1D_mixtures_sequential <- function(data, max_mixture, prior_BIC, verbose)
{
  fit_list <- list()

  for (param in colnames(data)) {
    if (verbose) {cat(param, " ")}

    dat <- data[,param]

    # Use Bayesian information criterion to choose best k
    BIC <- mclustBIC(dat, G = seq_len(max_mixture),
                     modelNames = "V", verbose = FALSE)

    # A tweak to favor smaller k
    for (k in seq_len(max_mixture)) {
      BIC[k] <- BIC[k] - prior_BIC*k*log(length(dat), base = 2)
    }

    # Fit the model with the chosen k
    fit <- Mclust(dat, x=BIC, verbose = FALSE)
    fit_list[[param]] <- fit$parameters
  }

  return(fit_list)
}




#########################
# Some plotting utilities
#########################

plot_distribution_1d <- function(dat, mixtures, name,
                                min = -2, max = 4,
                                separate = FALSE)
{
  mix <- mixtures[[name]]
  pts <- seq(from = min(dat[,name]), to = max(dat[,name]),
             length.out = 500)
  k <- length(mix$pro)

  # Vector of colors - assuming not more than 4 components
  colors <- c("blue", "red", "green", "orange")

  # Select the first mixture component
  m <- mix$mean[1]
  sd <- sqrt(mix$variance$sigmasq[1])
  g <- mix$pro[1] * dnorm(pts, mean = m, sd = sd)

  # Plot the first mixture component
  if (separate) {
    df <- data.frame(cbind(pts,g))
    names(df) <- c("x", "y")

    p <- ggplot() +
      geom_line(df, mapping = aes(x=.data$x, y=.data$y), color = "blue") +
      labs(title = "Mixture components", x = "", y = "") +
      theme_bw()

    # plot(pts, g, type = "l", xlim = c(min,max), ylim = c(0,1),
    #      xlab = name, ylab = "density",
    #      col = "blue", main = "components", cex.main = 2)
  }

  # Same for other components
  if (k > 1) {
    for (comp in c(2:k)) {
      m <- mix$mean[comp]
      sd <- sqrt(mix$variance$sigmasq[comp])

      g1 <-  mix$pro[comp] * dnorm(pts, mean = m, sd = sd)

      if (separate) {
        df$y <- g1
        p <- p + geom_line(df, mapping = aes(x=.data$x, y=.data$y),
                           color = colors[comp])
        # lines(pts, g1, col = colors[comp])
      } else
      {
        g <- g + g1
      }
    }
  }

  if (!separate) {
    df <- data.frame(cbind(pts,g))
    names(df) <- c("x", "y")

    p <- ggplot(df, aes(x=.data$x, y=.data$y)) +
      geom_line() +
      labs(title = "Gaussian mixture", x = "", y = "") +
      theme_bw()

  }
  return(p)
}



