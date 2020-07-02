as_matrix <- function(data) {
  if(is(data, "flowSet")) {
    data <- suppressWarnings(as(data, "flowFrame"))
    flowCore::exprs(data) <- flowCore::exprs(data)[,which(flowCore::colnames(data) != "Original")]
  }

  if(is(data,"flowFrame")) data <- flowCore::exprs(data)

  return(data)
}


start_parallel_cluster <- function() {
  cores <- detectCores()
  cl <- makeCluster(cores[1])
  registerDoParallel(cl)

  return(cl)
}
