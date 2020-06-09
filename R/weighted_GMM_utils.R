################################################################################
################################################################################
#
# Wrappers and utilities for mixture modelling in flow cytometry data
#
################################################################################
################################################################################


################################################################################
# Wrapper for weighted GMM algorithm
# Input:
#     data: a matrix containing events along the rows, variables along columns
#     k: the number of mixture components (clusters)
#     params: the variables to use during modeling; must be subset of colnames(data)
#     weights: a vector of length = nrow(data), specifying the weight of each event
#     initialize: a list of three elements, used to initialize the components:
#             pro: a vector of length = k, the proportions of the mixture
#             mean: a matrix with dim = (k, length(params));
#                   each row is the mean of one mixture component
#             variance$sigma: an array with dim = (length(params), length(params),k)
#                   each slice is the variance of one mixture component
#     seed: a seed for the random generator, only used in random initialization
#     max_steps: hard limit on number of steps for the weighted EM algorithm;
#             from my experience, more than the default of 10 leads to overfitting
#     relative_loglik_threshold: if the relative improvement of the log likelihood
#             objective function drops below threshold, terminate the algorithm
#             before max_steps are reached
#     regularize_variance: a boolean flag; if true, adds a small multiple of the
#             d x d identity matrix to each variance
#     variance_correction: variance lost in downsampling, added back at the end
#     verbose: boolean flag; if true, prints info about running time
# Output:
#     list with two elements, resulting from the EM algorithm:
#             mixture: a list describing the mixture; same format as `initialize`
#             event_probabilities: a matrix with dim = (nrow(data),k);
#                   the entry (ij) is the probability that event i is drawn from
#                   component j
################################################################################

bulk_weighted_gmm = function(data, k, params, weights = NULL,
                             initialize = NULL, seed = 91,
                             regularize_variance = TRUE, variance_correction = NULL,
                             max_steps = 10,
                             relative_loglik_threshold = 5e-3,
                             verbose = FALSE)
{
  start.time = Sys.time()

  n = nrow(data)

  if (is.null(initialize))
  {
    # Random initialization
    mixture = initialize_mixture(data,k,params,seed,verbose)
  }
  else
  {
    # Initialization given in function argument
    mixture = initialize
    mixture$pro = mixture$pro / sum(mixture$pro)
  }

  # If no weights are given, assume they are equal
  if (is.null(weights))
  {
    weights = rep(1, n)
  }

  # Normalize the weight vector
  sum_weights = sum(weights)
  weights = weights / sum_weights

  if(verbose) { cat("Total weight: ", sum_weights, "\n")}

  # Begin with one Expectation step; compute probabilities of assigning events
  # to the mixture components given in initialization
  e_step_result = e_step(data, mixture, params, verbose = FALSE)
  event_probabilities = e_step_result$probs
  probsum = e_step_result$probsum

  loglik = sum( log(probsum) * weights) * sum_weights / n

  # Variable that will keep track of log likelihood for each iteration
  loglik_history = loglik


  for (step in c(1:steps))
  {
    if(steps == 0) {break}

    if(verbose) {cat("Starting step ", step, "\n")}

    # Maximization step, with update rule modified to account for weights
    prev_mixture = mixture
    mixture = m_step(data, mixture, event_probabilities,
                     params, weights, variance_correction,
                     regularize_variance, verbose = FALSE)

    # Expectation step, identical to un-weighted EM algorithm
    e_step_result = e_step(data, mixture, params, verbose = FALSE)
    probsum = e_step_result$probsum

    loglik = sum( log(probsum) * weights) * sum_weights / n
    loglik_history = c(loglik_history, loglik)

    # If log likelihood decreased, stop the algorithm
    # This can only happen due to variance correction in m_step
    loglik_rel_change = abs((loglik - loglik_history[step]) / loglik_history[step])
    if ((step > 1) & (loglik_rel_change < 10*relative_loglik_threshold))
    {
      mixture = prev_mixture

      if (!is.null(variance_correction))
      {
        mixture = correct_variances(mixture, variance_correction, event_probabilities, weights, verbose)
      }

      break
    } else
    {
      event_probabilities = e_step_result$probs
    }
  }

  if(verbose) {print(Sys.time() - start.time)}

  list("mixture" = mixture, "event_probabilities" = event_probabilities,
       "loglik" = loglik_history)
}




################################################################################
# Utilities for weighted EM algorithm
################################################################################


# Quick and dirty random initialization, if no other is provided by user
# Works horribly for more than 2-3 dimensions
initialize_mixture = function(data, k, params, seed, verbose = FALSE)
{
  # Random initialization, with identity variance: modify later
  mixture = list()
  d = length(params)

  mixture$pro = rep(1/k, k)

  set.seed(seed)
  s = sample(nrow(data), k)
  mixture$mean = data[s,params,drop = FALSE]

  mixture$variance$sigma = array(NaN, c(d,d,k))
  for (component in c(1:k))
  {
    mixture$variance$sigma[,,component] = diag(d)
  }

  mixture
}


# Expectation step: given mixture parameters, update probabilities of event assignment
# Identical to unweighted version
e_step = function(data, mixture, params,
                  verbose = FALSE)
{
  start.time = Sys.time()
  n = nrow(data)
  k = length(mixture$pro)
  d = ncol(data)

  probs = matrix(0, nrow = n, ncol = k)

  # For each cluster, compute the probability that the datapoints belong to it
  for (cl in seq(k))
  {
    proportion = mixture$pro[cl]
    mean = mixture$mean[cl,]
    variance = mixture$variance$sigma[,,cl]

    if (d == 1)
    {
      # dnorm works with standard deviation, not variance
      probs[,cl] = proportion * dnorm(data, mean, sqrt(variance))
    }
    else
    {
      probs[,cl] = proportion * dmvnorm(data, mean, variance)
    }
  }


  # Normalize the probability along each row

  probsum = apply(probs, 1, sum)
  min = min(which(probsum > 0))
  probsum[probsum == 0] = min # To avoid division by 0
  probs = probs/probsum


  # Just for debugging, plot cluster assignments
  #assignment_2d = apply(probs, 1, which.max)
  #g = ggplot(data.frame(data),
  #           aes(x=CD3, y=CD19, color = as.factor(assignment_2d) ) ) + geom_point()
  #print(g)

  if(verbose) { print(Sys.time() - start.time) }

  list("probs" = probs, "probsum" = probsum)
}



# Maximization step: given event probabilities, update mixture parameters
# Modified from vanilla EM, to take weights into account, and correct for missing variance
m_step = function(data, mixture, event_probabilities, params,
                  weights, variance_correction, regularize_variance, verbose = FALSE)
{
  k = ncol(event_probabilities)
  d = ncol(data)
  n = nrow(data)
  weighted_probs = event_probabilities * weights

  sqrt_probs = sqrt(weighted_probs)


  mixture$pro = apply(weighted_probs, 2, sum)

  mixture$mean = t(weighted_probs)%*%data / mixture$pro

  if (verbose) { cat("Updating mixture parameters ") }

  for (slice in c(1:k))
  {

    if (d==1)
    {
      # TO DO: rewrite d == 1 case, this is junk
      out = shifted_data * shifted_data
      out = out * weighted_probs[,slice]

      summed = sum(out) / mixture$pro[slice]
      dim(summed) = c(1,1)
      mixture$variance$sigma[,,slice] = summed

    }
    else
    {
      scaled_data = sqrt_probs[,slice] * data
      scaled_mean = apply(scaled_data, 2, mean)

      scaled_var = var(scaled_data)
      var = n/mixture$pro[slice] * (scaled_var + outer(scaled_mean, scaled_mean) )
      var = var - outer(mixture$mean[slice,], mixture$mean[slice,])

      if(verbose & slice%%100 == 0)
      {
        cat(slice, " ")
      }

      # If flag is true, don't let variance get too close to 0
      if (regularize_variance)
      {
        var = var + 0.01 * diag(d)
      }


      mixture$variance$sigma[,,slice] = var

      rm(var, scaled_var, scaled_mean, scaled_data)
    }
  }

  mixture
}


# Variance correction: account for missing variance in each mixture component,
# arising from replacing a bunch of events with the center of their phenobin
correct_variances = function(mixture, variance_correction, event_probabilities, weights, verbose = FALSE)
{
  if (verbose)
  {
    print("Performing variance correction...")
  }
  k = ncol(event_probabilities)
  weighted_probs = event_probabilities * weights

  for (slice in c(1:k))
  {
    this_correction = apply(weighted_probs[,slice] * variance_correction, c(2:3), sum)
    this_correction = this_correction / mixture$pro[slice]

    mixture$variance$sigma[,,slice] = mixture$variance$sigma[,,slice] + this_correction
  }

  mixture
}



##########################
# Plotting utils
##########################


plot_mixture_assignments = function(data, params, k, assignment, title,
                                    outlier_threshold = 100,
                                    verbose = FALSE)
{
  start.time = Sys.time()

  # Some auxiliary objects
  tab_pred = tabulate(assignment)
  global_kdes = make_kdes_global(data, params)

  for (cl in c(1:k))
  {
    # If the cluster has too few events, they're probably outliers; discard
    if (tab_pred[cl] < outlier_threshold)
    {
      cat("Skipped cluster ", cl, ": ", tab_pred[cl], " events.\n")
      next
    }

    # Select the events from current cluster, and overlay their kde with the global kde
    sel = which(assignment == cl)
    l = length(sel)

    if (l < 1000)
    {
      events_cl = data[sel,]
    } else if (l < 5000)
    {
      events_cl = data[sample(sel,floor(l/2)),]
    } else
    {
      events_cl = data[sample(sel,floor(l/4)),]
    }

    plot_cluster_histograms(global_kdes, events_cl, params)

    dev.print(png, tight(title, "_cl", cl, ".png"), width = 600, height = 400)
  }

  if(verbose) { print(Sys.time() - start.time) }
}

