#####################################################################
#####################################################################
#
# Wrappers and utilities for mixture modelling in flow cytometry data
#
#####################################################################
#####################################################################


bulk_weighted_gmm <- function(data, k, params, weights = NULL, mixture = NULL,
                             regularize_variance = TRUE, variance_correction = NULL,
                             max_steps = 10, relative_loglik_threshold = 5e-3,
                             verbose = FALSE)
{
  n <- nrow(data)

  if (is.null(mixture)) {
    mixture <- initialize_mixture(data,k,params,verbose)
  } else {
    mixture$pro <- mixture$pro / sum(mixture$pro)
  }

  if (is.null(weights)) { weights <- rep(1, n) }
  sum_weights = sum(weights)
  weights <- weights / sum_weights

  e_step_result <- e_step(data, mixture, params, verbose = FALSE)
  event_probabilities <- e_step_result$probs
  probsum <- e_step_result$probsum

  loglik <- sum( log(probsum) * weights) * sum_weights / n
  loglik_history <- loglik

  for (step in seq_len(max_steps)) {
    if(verbose) {cat("Starting step ", step, "\n")}
    prev_mixture <- mixture
    mixture <- m_step(data, mixture, event_probabilities,
                     params, weights, variance_correction,
                     regularize_variance, verbose = FALSE)

    e_step_result <- e_step(data, mixture, params, verbose = FALSE)
    probsum <- e_step_result$probsum

    loglik <- sum( log(probsum) * weights) * sum_weights / n
    loglik_history <- c(loglik_history, loglik)

    loglik_rel_change <- abs((loglik - loglik_history[step]) / loglik_history[step])
    if ((step > 1) & (loglik_rel_change < 10*relative_loglik_threshold)) {
      mixture <- prev_mixture

      if (!is.null(variance_correction)) {
        mixture <- correct_variances(mixture, variance_correction, event_probabilities, weights, verbose)
      }
      break
    } else { event_probabilities <- e_step_result$probs }
  }

  return(mixture)
}




################################################################
# Utilities for weighted EM algorithm
################################################################


# Quick and dirty random initialization, if not provided by user
# Works horribly for more than 2-3 dimensions
initialize_mixture <- function(data, k, params, verbose = FALSE)
{
  # Random initialization, with identity variance: modify later
  mixture <- list()
  d <- length(params)

  mixture$pro <- rep(1/k, k)

  s <- sample(nrow(data), k)
  mixture$mean <- data[s,params,drop = FALSE]

  mixture$variance$sigma <- array(NaN, c(d,d,k))
  for (component in seq_len(k)) {
    mixture$variance$sigma[,,component] <- diag(d)
  }

  return(mixture)
}


# Expectation step: given mixture parameters,
# update probabilities of event assignment.
# Identical to unweighted version
e_step <- function(data, mixture, params,
                  verbose = FALSE)
{
  n <- nrow(data)
  k <- length(mixture$pro)
  d <- ncol(data)

  probs <- matrix(0, nrow = n, ncol = k)

  # For each cluster, compute the probability that the
  # datapoints belong to it
  for (cl in seq(k)) {
    proportion <- mixture$pro[cl]
    mean <- mixture$mean[cl,]
    variance <- mixture$variance$sigma[,,cl]

    if (d == 1) {
      # dnorm works with standard deviation, not variance
      probs[,cl] <- proportion * dnorm(data, mean, sqrt(variance))
    }
    else {
      probs[,cl] <- proportion * dmvnorm(data, mean, variance)
    }
  }

  # Normalize the probability along each row
  probsum <- apply(probs, 1, sum)
  min <- min(which(probsum > 0))
  probsum[probsum == 0] <- min # To avoid division by 0
  probs <- probs/probsum

  l <- list("probs" = probs, "probsum" = probsum)
  return(l)
}



# Maximization step: given event probabilities, update mixture parameters
# Modified from vanilla EM, to take weights into account,
# and correct for missing variance
m_step <- function(data, mixture, event_probabilities, params,
                  weights, variance_correction,
                  regularize_variance, verbose = FALSE)
{
  k <- ncol(event_probabilities)
  d <- ncol(data)
  n <- nrow(data)

  weighted_probs <- event_probabilities * weights
  sqrt_probs <- sqrt(weighted_probs)

  mixture$pro <- apply(weighted_probs, 2, sum)
  mixture$mean <- t(weighted_probs)%*%data / mixture$pro

  for (slice in seq_len(k)) {
    if (d==1)
    {
      # TO DO: rewrite d == 1 case
    }
    else {
      scaled_data <- sqrt_probs[,slice] * data
      scaled_mean <- apply(scaled_data, 2, mean)

      scaled_var <- var(scaled_data)
      var <- n/mixture$pro[slice] * (scaled_var +
                                    outer(scaled_mean, scaled_mean) )
      var <- var - outer(mixture$mean[slice,], mixture$mean[slice,])

      if (regularize_variance) {
        var <- var + 0.01 * diag(d)
      }

      mixture$variance$sigma[,,slice] <- var
      rm(var, scaled_var, scaled_mean, scaled_data)
    }
  }

  return(mixture)
}


# Variance correction: account for missing variance
# in each mixture component, arising from replacing a
# bunch of events with the center of their bin
correct_variances <- function(mixture, variance_correction,
                              event_probabilities, weights,
                              verbose = FALSE)
{
  if (verbose) {
    print("Performing variance correction...")
  }

  k <- ncol(event_probabilities)
  weighted_probs <- event_probabilities * weights

  for (slice in seq_len(k)) {
    this_correction <- apply(weighted_probs[,slice] * variance_correction,
                             c(2:3), sum)
    this_correction <- this_correction / mixture$pro[slice]

    mixture$variance$sigma[,,slice] <- mixture$variance$sigma[,,slice] +
      this_correction
  }

  return(mixture)
}


