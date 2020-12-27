#################################################
# Utilities for visualization and plotting: kdes,
# t-SNE reduction to 2D of cluster centroids
#################################################


make_kdes_global <- function(data, parameters)
{
  # Make kde of global distribution for each parameter
  kdes <- list()
  for (param in parameters) {
    kdes[[param]] <- bkde(data[,param])
  }

  return(kdes)
}


plot_kde_vs_mixture <- function(data, global_kdes, mixtures, name,
                              xmin = -2, xmax = 5)
{
  df <- data.frame(global_kdes[[name]])
  p1 <- ggplot(df, aes(x=.data$x, y=.data$y)) +
    geom_line() +
    labs(title = paste(name, "kde"), x = "", y = "") +
    theme_bw()

  # Plot gaussian mixture
  p2 <- plot_distribution_1d(data, mixtures, name = name,
                             min = xmin, max = xmax)

  # # Plot mixture components separately
  p3 <- plot_distribution_1d(data, mixtures, name = name,
                            min = xmin, max = xmax, separate = TRUE)

  return(gridExtra::grid.arrange(p1, p2, p3, ncol = 3))
}

# Overlay global kde with cluster kde, for each parameter
plot_cluster_histograms <- function(global_kdes, cluster = NULL,
                                   cluster2 = NULL, parameters,
                                   weights = NULL,
                                   overlay_hist = TRUE)
{
  i <- 1
  while (i * (i+1) < length(parameters)) {
    i <- i+1
  }
  par(mfrow = c(i,i+1), mar = c(1,0,3,0) + 0.1)

  for (parameter in parameters) {
    # Plot kernel density estimates for each of the parameters
    if (is.null(cluster2)) {
      kde <- global_kdes[[parameter]]
    } else {
      kde <- bkde(cluster2[,parameter])
    }
    plot(kde, type = "l", main = parameter)

    if (overlay_hist) {
      if (is.null(weights)) {
        kde1 <- bkde(cluster[,parameter])
      } else {
        kde1 <- density(cluster[,parameter], weights = weights)
      }
      scal <- 0.8 * max(kde$y) / max(kde1$y)
      kde1$y <- kde1$y * scal
      lines(kde1, col = "red")
    }
  }
}


get_tsne_centers <- function(data, probs)
{
  # Get t-SNE reduction to 2D of bin centers,
  # and map each bin to its most likely cluster

  n_items <- nrow(data)
  perplexity <- min((n_items - 1)/3, 30)

  res <- Rtsne(data, perplexity = perplexity)$Y
  cluster <- apply(probs, 1, which.max)
  res <- cbind(res, cluster)
  colnames(res) <- c("tsne_1", "tsne_2", "cluster")

  return(data.frame(res))
}


get_tsne_clusters <- function(tailor_obj)
{
  # Get t-SNE reduction to 2D of cluster centroids
  centers <- tailor_obj$cat_clusters$centers

  n_items <- nrow(centers)
  perplexity <- min((n_items - 1)/3, 30)

  res <- Rtsne(dist(centers), perplexity = perplexity)$Y
  colnames(res) <- c("tsne_1", "tsne_2")

  return(data.frame(res))
}


get_legend <- function(my_plot)
{
  tmp <- ggplot_gtable(ggplot_build(my_plot))
  leg <- which(vapply(tmp$grobs, function(x) x$name, character(1)) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}



compute_labels <- function(obj, defs) {
  k <- length(obj$func)
  labels <- sapply(seq(1,k), function(i) cluster_phenotype_label(obj,i, defs))
  return(labels)
}

cluster_phenotype_label = function(obj, cluster, defs = NULL, ignore_markers = c()) {
  phen = obj$phenotypes[cluster,]
  hi = names(phen)[which(phen == 3)]
  med = names(phen)[which(phen == 2)]
  lo = names(phen)[which(phen == 1)]
  func = obj$func[cluster]


  if (!is.null(defs)) {
    ex = colnames(defs)
    idx = which(defs[func, ] == "dc")

    if (length(idx > 0)) {
      ex = ex[-idx]
    }

    ex = c(ex, ignore_markers)
    incl_hi = hi[!(hi %in% ex)]
    incl_med = med[!(med %in% ex)]
    incl_lo = lo[!(lo %in% ex)]
  } else {
    incl_hi = hi
    incl_med = med
    incl_lo = lo
  }

  if(length(incl_hi) > 0) {
    hit <- paste(sep = "", incl_hi, "++")
  } else {
    hit <- numeric(0)
  }
  if(length(incl_med) > 0) {
    medt <- paste(sep = "", incl_med, "+")
  } else {
    medt <- numeric(0)
  }
  if(length(incl_lo) > 0) {
    lot <- paste(sep = "", incl_lo, "-")
  } else {
    lot <- numeric(0)
  }

  res <- func
  if (length(hit) > 0) {
    res <- paste(c(func, hit), collapse = " ")
  }
  if (length(medt) > 0) {
    res <- paste(c(func, medt), collapse = " ")
  }
  # if (length(lot) > 0) {
  #   res <- paste(c(res , lot), collapse = " ")
  # }

  return(res)
}





