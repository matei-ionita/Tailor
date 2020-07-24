as_matrix <- function(data) {
  if (is(data, "flowSet")) {
    data <- flowCore::exprs(suppressWarnings(as(data, "flowFrame")))
    data <- data[,which(colnames(data) != "Original")]
  }

  if (is(data, "flowFrame")) {
    data <- flowCore::exprs(data)
  }
  return(data)
}


start_parallel_cluster <- function() {
  cores <- detectCores()
  cl <- makeCluster(cores[1])
  registerDoParallel(cl)

  return(cl)
}


tailor_map_parallel <- function(data, tailor_obj, n_batch)
{
  k <- length(tailor_obj$mixture$pro)
  logpro <- log(tailor_obj$mixture$pro)

  n <- nrow(data)
  params = colnames(data)
  mapping <- integer(length = n)
  batch_size <- ceiling(n/n_batch)

  clust <- start_parallel_cluster()

  data_list <- list()
  for (batch in seq_len(n_batch))
  {
    start_batch <- (batch - 1) * batch_size + 1
    end_batch <- batch * batch_size
    if (batch == n_batch) { end_batch <- n }

    data_list[[batch]] <- data[c(start_batch:end_batch),params]
  }
  rm(data)

  batch_data <- NULL
  mapping <- foreach(batch_data = iter(data_list), .combine = c, .packages = c("mvtnorm")) %dopar%
    {
      probs <- matrix(0, nrow = nrow(batch_data), ncol = k)

      for (cl in seq(k))
      {
        mean <- tailor_obj$mixture$mean[cl,]
        sigma <- tailor_obj$mixture$variance$sigma[,,cl]
        probs[,cl] <- logpro[cl] + dmvnorm(batch_data, mean, sigma, log = TRUE)
      }

      # Assign each datapoint to the cluster of maximum probability
      result <- apply(probs, 1, which.max)
      rm(probs)
      gc()

      result
    }

  stopCluster(clust)
  return(mapping)
}



tailor_map_sequential <- function(data, tailor_obj, n_batch, verbose)
{
  k <- length(tailor_obj$mixture$pro)
  logpro <- log(tailor_obj$mixture$pro)

  n <- nrow(data)
  params = colnames(data)
  mapping <- integer(length = n)
  batch_size <- ceiling(n/n_batch)

  for (batch in seq_len(n_batch)) {
    if(verbose) {cat(batch, " ")}
    start_batch <- (batch - 1) * batch_size + 1
    end_batch <- batch * batch_size
    if (batch == n_batch) { end_batch = n }

    batch_data <- data[c(start_batch:end_batch),params]

    # For each cluster, compute the probability that data in current batch
    # are drawn from it
    probs <- matrix(0, nrow = end_batch - start_batch + 1, ncol = k)

    for (cl in seq(k)) {
      weight <- tailor_obj$mixture$pro[cl]
      mean <- tailor_obj$mixture$mean[cl,]
      sigma <- tailor_obj$mixture$variance$sigma[,,cl]

      probs[,cl] <- weight * dmvnorm(batch_data, mean, sigma)
    }

    # Assign each datapoint to the cluster of maximum probability
    mapping[c(start_batch:end_batch)] <- apply(probs, 1, which.max)
  }

  return(mapping)
}


binary_search <- function(f, range, val = 0, error = 1e-4) {
  lower <- range[1]
  upper <- range[2]
  mid   <- (lower + upper)/2
  trend <- ifelse(f(upper) >= f(lower), 1, -1)

  if (abs(f(mid) - val) < error) return(mid)
  if (upper - mid < error) return(upper)
  if (mid - lower < error) return(lower)

  diff <- trend * (f(mid) - val)
  if (diff < 0) {
    return(binary_search(f, c(mid, upper)))
  }
  else {
    return(binary_search(f, c(lower, mid)))
  }
}


