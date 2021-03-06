

#################################################
# Utilities for categorical merging:
# merge mixture components which are unlikely to
# represent biologically distinct populations
#################################################


##############################################
# Take as input a list of 1D mixture models
# Output a list of cutoffs for each parameter:
# a cutoff is the threshold between distinct
# Gaussian mixture components
##############################################

get_1D_cutoffs <- function(mixtures, params)
{
  cutoffs <- list()

  for (param in params) {
    mix <- mixtures[[param]]
    n_cutoffs <- length(mix$pro) - 1

    if (n_cutoffs == 0) {
      cutoffs[[param]] <- NULL
      next
    }
    cutoffs[[param]] <- integer(0)

    for (i1 in seq(n_cutoffs)) {
      i2 <- i1 + 1
      pos <- get_single_cutoff(mix, i1, i2)
      if (pos < 0.5)
        next
      l <- length(cutoffs[[param]])
      if (l == 0)
        cutoffs[[param]] <- c(cutoffs[[param]],pos)
      else if (cutoffs[[param]][l] < pos)
        cutoffs[[param]] <- c(cutoffs[[param]],pos)
    }
  }

  return(cutoffs)
}

get_single_cutoff <- function(mix, i1, i2, er = 1e-4) {
  pro1 <- mix$pro[i1]
  pro2 <- mix$pro[i2]
  mean1 <- mix$mean[i1]
  mean2 <- mix$mean[i2]
  sd1 <- sqrt(mix$variance$sigmasq[i1])
  sd2 <- sqrt(mix$variance$sigmasq[i2])

  f <- function(x) pro1 * dnorm(x,mean1,sd1) - pro2 * dnorm(x,mean2,sd2)
  pos <- binary_search(f, range = c(mean1, mean2), val = 0, error = er)

  if (abs(pos - mean1) < er) {
    pos <- binary_search(f, range = c(mean1-sd1, mean1), val = 0, error = er)
  }
  if (abs(pos - mean2) < er) {
    pos <- binary_search(f, range = c(mean2, mean2+sd2), val = 0, error = er)
  }

  return(unname(pos))
}


print_kdes_with_cutoffs <- function(data, cutoffs, parameters)
{
  data = as_matrix(data)
  global_kdes <- make_kdes_global(data, parameters)

  # determine size of the grid of plots
  i <- 1
  while (i * (i+1) < length(parameters))
  {
    i <- i+1
  }
  par(mfrow = c(i,i+1), mar = c(1,0,3,0) + 0.1)

  for (parameter in parameters) {
    kde <- global_kdes[[parameter]]
    plot(kde, type = "l", main = parameter)

    # We model some markers as unimodal: no cutoff
    if(!is.null(cutoffs[[parameter]])) {
      abline(v = cutoffs[[parameter]])
    }
  }
}


###################################################
# Take input the means and mixture proportions
# of mixture components, as well as cutoffs between
# biologically distinct modes for each marker.
# Merge mixture components if they lie on the same
# side of the cutoff, for all markers.
# Output the centroid and putative phenotype of
# each categorical cluster, as well as mapping from
# mixture components to categorical clusters.
###################################################
categorical_merging <- function(pros, means, sigmas, cutoffs, params)
{
  n <- nrow(means)
  categs <- vector(mode = "character", length = n)
  cat_mapping <- vector(mode = "integer", length = n)

  for (param in params) {
    if (is.null(cutoffs[[param]])) {
      for (i in seq(n)) {
        categs[i] <- paste(categs[i], 1, sep = "")
      }
      next
    }

    for (i in seq(n)) {
      j <- 1
      while(j <= length(cutoffs[[param]]) && means[i,param] > cutoffs[[param]][j])
        j <- j+1
      categs[i] <- paste(categs[i], j, sep = "")
    }
  }

  un <- unique(categs)

  for (i in seq(n)) {
    cat_mapping[i] <- which(un == categs[i])
  }

  centers <- get_cluster_means(un, params, cat_mapping, pros, means)
  sds <- get_cluster_sds(un, params, cat_mapping, pros, means, sigmas, centers)
  phenotypes <- get_unique_phenotypes(params, un)
  # cutoffs <- vapply(cutoffs, unname, numeric(1))

  l <- list(phenotypes = phenotypes, mixture_to_cluster = cat_mapping,
            centers = centers, sds = sds, cutoffs = cutoffs)
  return(l)
}

get_cluster_means <- function(unique, params, cat_mapping, pros, means) {
  centers <- matrix(0, nrow = length(unique), ncol = length(params))
  colnames(centers) <- params
  for (i in seq(length(unique))) {
    weight <- 0
    components <- which(cat_mapping == i)
    for (comp in components) {
      centers[i,] <- centers[i,] +  pros[comp] * means[comp,]
      weight <- weight + pros[comp]
    }
    centers[i,] <- centers[i,] / weight
  }

  return(centers)
}


get_cluster_sds <- function(unique, params, cat_mapping, pros, means, sigmas, centers) {
  sds <- matrix(0, nrow = length(unique), ncol = length(params))
  colnames(sds) <- params
  for (i in seq(length(unique))) {
    weight <- 0
    var <- matrix(0, nrow = length(params), ncol = length(params))
    components <- which(cat_mapping == i)
    for (comp in components) {
      diff <- means[comp,] - centers[i,]
      var <- var + pros[comp] * sigmas[,,comp] + pros[comp] * outer(diff, diff)
      weight <- weight + pros[comp]
    }
    sds[i,] <- sqrt(diag(var / weight))
  }

  return(sds)
}


get_unique_phenotypes <- function(names, unique) {
  split <- strsplit(unique, "")
  phenotypes <- do.call(rbind, split)
  phenotypes <- data.frame(phenotypes)
  names(phenotypes) <- names

  n_levels <- vapply(phenotypes, function(x) length(unique(x)), integer(1))
  phenotypes <- phenotypes[,which(n_levels > 1)]

  return(phenotypes)
}





