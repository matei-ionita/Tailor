#################################################
# Utilities for visualization and plotting: kdes,
# t-SNE reduction to 2D of cluster centroids
#################################################


make_kdes_global = function(data, parameters)
{
  # Make kde of global distribution for each parameter
  kdes = list()
  for (param in parameters)
  {
    kdes[[param]] = bkde(data[,param])
  }

  kdes
}


plot_kde_vs_mixture = function(data, global_kdes, mixtures, name,
                              xmin = -2, xmax = 5)
{
  # Set up the stage for 3 plots side by side
  par(mfrow = c(1,3), mar = c(5, 4, 4, 1))

  # Plot kde of the desired variable
  plot(global_kdes[[name]], type = "l", xlim = c(xmin,xmax), ylim = c(0,1),
       xlab = name, ylab = "density", main = "kde", cex.main = 2)

  # Plot gaussian mixture
  plot_distribution_1d(data, mixtures, name = name,
                       min = xmin, max = xmax)

  # Plot mixture components separately
  plot_distribution_1d(data, mixtures, name = name,
                       min = xmin, max = xmax, separate = TRUE)
}


plot_cluster_histograms = function(global_kdes, cluster = NULL,
                                   cluster2 = NULL, parameters,
                                   weights = NULL,
                                   overlay_hist = TRUE)
{
  # Overlay global kde with cluster kde, for each parameter

  i = 1
  while (i * (i+1) < length(parameters))
  {
    i = i+1
  }

  par(mfrow = c(i,i+1), mar = c(1,0,3,0) + 0.1)

  for (parameter in parameters)
  {
    # Plot kernel density estimates for each of the parameters
    if (is.null(cluster2))
    {
      kde = global_kdes[[parameter]]
    }
    else
    {
      kde = bkde(cluster2[,parameter])
    }
    plot(kde, type = "l", main = parameter)

    if (overlay_hist)
    {
      if (is.null(weights))
      {
        kde1 = bkde(cluster[,parameter])
      }
      else
      {
        kde1 = density(cluster[,parameter], weights = weights)
      }
      scal = 0.8 * max(kde$y) / max(kde1$y)
      kde1$y = kde1$y * scal
      lines(kde1, col = "red")
    }

    gc()
  }
}




get_tsne_centers = function(data, probs, seed = 137)
{
  # Get t-SNE reduction to 2D of bin centers,
  # and map each bin to its most likely cluster

  n_items = nrow(data)
  perplexity = min((n_items - 1)/3, 30)
  set.seed(seed)

  res = Rtsne(data, perplexity = perplexity)$Y
  cluster = apply(probs, 1, which.max)
  res = cbind(res, cluster)

  colnames(res) = c("tsne_1", "tsne_2", "cluster")

  data.frame(res)
}
