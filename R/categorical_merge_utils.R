

#################################################
# Utilities for categorical merging:
# merge mixture components which are unlikely to
# represent biologically distinct populations
#################################################


get_1D_cutoffs <- function(mixtures, to_merge, params)
{
  ##############################################
  # Take as input a list of 1D mixture models
  # Output a list of cutoffs for each parameter:
  # a cutoff is the threshold between distinct
  # Gaussian mixture components
  ##############################################

  cutoffs <- list()

  # Ignore distinction between mixture components
  # specified in to_merge; we don't expect these
  # to be biologically different
  for (param in params)
  {
    mix <- mixtures[[param]]
    if (length(mix$pro) == 1)
    {
      cutoffs[[param]] <- NULL
      next
    }

    i1 <- 1
    i2 <- 2

    if (length(mix$pro) == 3)
    {
      i1 <- 2
      i2 <- 3
    }

    if (!is.null(to_merge[[param]]))
    {
      if (to_merge[[param]][1] == 2)
      {
        i1 <- 1
        i2 <- 2
      }
    }

    pro1 <- mix$pro[i1]
    pro2 <- mix$pro[i2]
    mean1 <- mix$mean[i1]
    mean2 <- mix$mean[i2]
    sd1 <- sqrt(mix$variance$sigmasq[i1])
    sd2 <- sqrt(mix$variance$sigmasq[i2])

    # look for point where second gaussian becomes larger than the first
    pos <- mean1
    while (pro1 * dnorm(pos, mean = mean1, sd = sd1)
           > pro2 * dnorm(pos, mean = mean2, sd = sd2))
    {
      pos <- pos + 0.05
    }

    cutoffs[[param]] <- pos
  }

  cutoffs
}


print_kdes_with_cutoffs <- function(data, cutoffs, parameters)
{
  ###################################################
  # print kdes with vertical lines displaying cutoffs
  ###################################################

  if(is(data, "flowSet"))
  {
    data <- suppressWarnings(as(data, "flowFrame"))
    flowCore::exprs(data) <- flowCore::exprs(data)[,which(flowCore::colnames(data) != "Original")]
  }

  if(is(data, "flowFrame")) data <- flowCore::exprs(data)


  global_kdes <- make_kdes_global(data, parameters)

  # determine size of the grid of plots
  i <- 1
  while (i * (i+1) < length(parameters))
  {
    i <- i+1
  }

  par(mfrow = c(i,i+1), mar = c(1,0,3,0) + 0.1)

  for (parameter in parameters)
  {
    kde <- global_kdes[[parameter]]
    plot(kde, type = "l", main = parameter)

    # We model some markers as unimodal: no cutoff
    if(!is.null(cutoffs[[parameter]]))
    {
      abline(v = cutoffs[[parameter]])
    }
  }
}



categorical_merging <- function(pros, means, cutoffs, params)
{
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

  n <- nrow(means)
  categs <- vector(mode = "character", length = n)
  cat_mapping <- vector(mode = "integer", length = n)

  for (param in params)
  {

    if (is.null(cutoffs[[param]]))
    {
      for (i in seq(n))
      {
        categs[i] <- paste(categs[i], "-", sep = "")
      }
      next
    }

    v <- means[,param] - cutoffs[[param]]
    for (i in seq(n))
    {
      if (v[i]<0)
      {
        categs[i] <- paste(categs[i], "-", sep = "")
      }
      else
      {
        categs[i] <- paste(categs[i], "+", sep = "")
      }
    }
  }

  un <- unique(categs)
  for (i in seq(n))
  {
    cat_mapping[i] <- which(un == categs[i])
  }

  centers <- matrix(0, nrow = length(un), ncol = length(params))
  colnames(centers) <- params
  for (i in seq(length(un)))
  {
    weight <- 0
    components <- which(cat_mapping == i)
    for (comp in components)
    {
      centers[i,] <- centers[i,] +  pros[comp] * means[comp,]
      weight <- weight + pros[comp]
    }
    centers[i,] <- centers[i,] / weight
  }

  list(phenotypes = un, mixture_to_cluster = cat_mapping, centers = centers)
}








