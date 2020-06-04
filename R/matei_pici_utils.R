######################### 
# Metacluster percentages 
######################### 

# This is Wade's function, with minor modifications: 
# Accepting any clustering algorithm, not necessarily the one in FlowSOM 
# Speeding up by removing stuff from inner loop 


cluster_percentages = function(counts, mapping, verbose = FALSE) 
{ 
  n_instances = length(counts) 
  n_clust = max(mapping) 
  n_data = length(mapping) 
  
  pctg = matrix(NA, nrow = n_instances, ncol = n_clust) 
  colnames(pctg) = paste("c", 1:n_clust, sep = "") 
  
  start = 1 
  for (i in 1:n_instances) 
  { 
    if(verbose) { cat("Instance", i, "\n") } 
    size = counts[i] 
    
    sample_mapping = mapping[start : (start + size - 1)] 
    hst = hist(sample_mapping, breaks = (0:n_clust) + 0.5, plot = FALSE)$counts 
    
    pctg[i,] = hst / size 
    
    start = start + size 
  } 
  
  pctg 
} 


# for each instance in the fsom data, calculate percentage of events in each 
# metacluster 
# Input: 
# fsom, a FlowSOM object 
# metaClustering, a factor containing cluster assignments for SOM nodes 
# Value:  matrix, rows are instances, cols are metaclusters (rowsums = 1) 
metacluster_percentages_fast = function(fsom, metaClustering, verbose = FALSE)  
{ 
  n_instances = length(fsom$FlowSOM$metaData) 
  n_clust = length(levels(metaClustering)) 
  n_nodes = length(metaClustering) 
  pctg = matrix(NA, nrow = n_instances, ncol = n_clust) 
  colnames(pctg) = paste("c", 1:n_clust, sep = "") 
  
  for (i in 1:n_instances)  
  { 
    if (verbose) {cat("Instance", i)} 
    strt = fsom$FlowSOM$metaData[[i]][1] 
    fini = fsom$FlowSOM$metaData[[i]][2] 
    n_events = 1 + fini - strt 
    
    # events per cluster 
    dat = fsom$FlowSOM$map$mapping[strt:fini,1] 
    dis = fsom$FlowSOM$map$mapping[strt:fini,2] 
    dat_no_outliers = dat[which(dis<4)] 
    hst = hist(dat_no_outliers, breaks = (0:n_nodes) + .5, plot = FALSE)$counts 
    
    for (j in 1:n_clust)  
    { 
      #if (verbose) {cat(j, " ")} 
      # resolve nodes to metaclusters 
      pctg[i, j] = sum(hst[which(as.numeric(metaClustering) == j)]) / n_events 
    } 
    if (verbose) {cat("\n")} 
  } 
  
  pctg 
} 



######################################################################### 
# Hierarchical clustering, based on distance and dip test for unimodality 
######################################################################### 

# Compute distance matrix 
# Input: centroids, a matrix whose rows are the parameter values for the cluster centroids 
# Value: distMatrix 

compute_distance_matrix = function(centroids)  
{ 
  n_clust = nrow(centroids) 
  distMatrix = matrix(Inf, nrow = n_clust, ncol = n_clust) 
  
  for (i in seq(n_clust-1))  
  { 
    for (j in c((i+1):n_clust)) 
    { 
      distMatrix[i,j] = dist(rbind(centroids[i,], centroids[j,])) 
    } 
  } 
  
  distMatrix 
} 

# remark: for the precomputed SOM nodes, distances go up to 10 


# Update distance matrix, after the merger of two clusters 
# Assuming that centroids has been updated already, and new cluster is the last row 

update_distance_matrix = function(centroids, distMatrix, c1, c2) 
{ 
  n_clust = nrow(centroids) 
  newDist = numeric(length = n_clust-1) 
  newCentroid = centroids[n_clust,] 
  
  # Compute distances to newly created cluster 
  for (i in seq(n_clust-1)) 
  { 
    newdist = dist(rbind(centroids[i,], newCentroid)) 
    
    if (i < c1) 
    { 
      newDist[i] = max(newdist, distMatrix[i,c1], distMatrix[i,c2]) 
      next 
    } 
    
    if (i + 1 < c2) 
    { 
      newDist[i] = max(newdist, distMatrix[c1,i+1], distMatrix[i+1,c2]) 
      next 
    } 
    
    newDist[i] = max(newdist, distMatrix[c1,i+2], distMatrix[c2,i+2]) 
  } 
  
  # Drop rows, cols corresponding to clusters that were merged 
  distMatrix = distMatrix[-c(c1,c2), -c(c1,c2)] 
  
  # Add row, col corresponding to new cluster 
  distMatrix = cbind(distMatrix, newDist, deparse.level = 0) 
  distMatrix = rbind(distMatrix, rep(Inf, n_clust)) 
  
  distMatrix 
} 



# Compute Hartigans' dip test, do determine whether the distribution 
# resulting from merger of two clusters is unimodal. 
# Small values indicate unimodal distribution; 0.02 is a reasonable threshold 

# If plot = FALSE, returns the numeric value of the dip test 
# If plot = TRUE, prints the numeric value and returns the plot of the distribution  

dip_test = function(events1, events2, centroid1, centroid2, l1, l2, 
                    lda = FALSE, plot = FALSE, p_value = FALSE) 
{ 
  set.seed(999) 
  #l = min(l1,l2, max(floor(l1/10),floor(l2/10))) 
  l = min(l1,l2) 
  if (l1 > 10000 & l2 > 10000) { l = min(l, max(floor(l1/10), floor(l2/10)))} 
  events1 = events1[sample(seq(l1),l),] 
  events2 = events2[sample(seq(l2),l),] 
  # Compute direction of maximal variation between clusters 
  if (lda) 
  { 
    # Use the first LDA component: more accurate, but very slow 
    events = rbind(data.frame(val = events1, clust = "1"),  
                   data.frame(val = events2, clust = "2") ) 
    diff = lda(clust ~ ., events, prior = c(1,1)/2)$scaling[,1] 
    remove(events) 
  } 
  else 
  { 
    # Poor man's LDA: difference between centroids 
    diff = centroid2 - centroid1 
  } 
  # Project data onto the difference vector 
  dat1 = data.frame(val = events1%*%diff, clust = "1") 
  dat2 = data.frame(val = events2%*%diff, clust = "2") 
  dat = rbind(dat1, dat2) 
  
  # Compute dip test 
  if (p_value) { dip_value = dip.test(dat[,"val"]) } 
  else { dip_value = dip(dat[,"val"]) } 
  
  if (plot) 
  { 
    #cat("Dip value:", dip_value, "\n") 
    print(dip_value) 
    # Sanity check: plot the distribution 
    ggplot(dat, aes(x=val, fill=clust)) + 
      geom_histogram() 
  } 
  else { dip_value } 
} 


# A tool to evaluate the performance of the dip test 
evaluate_dip = function(fsom, clara, c1, c2, p_value = FALSE) 
{ 
  if(!is.null(fsom)) 
  { 
    cols = fsom$FlowSOM$map$colsUsed 
    centroids = fsom$FlowSOM$map$codes 
    mapping = fsom$FlowSOM$map$mapping[,1] 
    
    events1 = fsom$FlowSOM$data[which(mapping == c1),cols] 
    events2 = fsom$FlowSOM$data[which(mapping == c2),cols] 
  } 
  else if(!is.null(clara)) 
  { 
    centroids <- clara$medoids 
    mapping <- clara$clustering 
    
    events1 = clara$data[which(mapping == c1),] 
    events2 = clara$data[which(mapping == c2),] 
  } 
  
  dip_test(events1, events2, centroids[c1,], centroids[c2,],  
           nrow(events1), nrow(events2), 
           lda = TRUE, plot = TRUE, p_value = p_value) 
} 



# main function of hierarchical algorithm 

# Input: 
# fsom, a FlowSOM object, with nodes that need to be clustered 
# kmin, an integer giving the minimum allowed number of clusters 
# dip_threshold, a merger is prevented unless the dip test is below the threshold 
# dist_absolute_threshold, exit the algorithm if the minimum distance between two clusters 
# becomes greater than this threshold 

# Value: node_mapping, a list with cluster membership of the SOM nodes 
hierarchical_clustering_dip = function(fsom, 
                                       clara = NULL, 
                                       kmin = 1, 
                                       dip_threshold = 0.03, 
                                       dist_absolute_threshold = 2.4, 
                                       perform_dip = TRUE, 
                                       dip_lda = FALSE, 
                                       verbose = FALSE) 
{ 
  # Various initializations: 
  if (!is.null(fsom)) 
  { 
    cols <- fsom$FlowSOM$map$colsUsed 
    centroids <- fsom$FlowSOM$map$codes 
    mapping <- fsom$FlowSOM$map$mapping 
    
    n_clust = nrow(centroids) 
    
    if (verbose) { cat("Precomputing lists of events:", "\n")} 
    events_in_nodes = vector(mode = "list", length = n_clust) 
    for (node in seq(n_clust)) 
    { 
      if (verbose & node %% 10 == 0) { cat(node, " ") } 
      events_in_nodes[[node]] = fsom$FlowSOM$data[which(mapping[,1] == node & mapping[,2] < 4),cols] 
    } 
    if (verbose) cat("\n", "Done precomputing.", "\n") 
  } 
  else if (!is.null(clara)) 
  { 
    cols <- colnames(clara$data) 
    centroids <- clara$medoids 
    mapping <- clara$clustering 
    
    n_clust = nrow(centroids) 
    
    if (verbose) { cat("Precomputing lists of events:", "\n")} 
    events_in_nodes = vector(mode = "list", length = n_clust) 
    for (node in seq(n_clust)) 
    { 
      if (verbose & node %% 10 == 0) { cat(node, " ") } 
      events_in_nodes[[node]] = clara$data[which(mapping == node),] 
    } 
    if (verbose) cat("\n", "Done precomputing.", "\n") 
  } 
  
  
  distMatrix = compute_distance_matrix(centroids) 
  
  n_clust_original = n_clust 
  node_mapping = seq(n_clust) 
  
  
  # A flag to keep track of algorithm exit conditions 
  done_algorithm = FALSE 
  
  while (n_clust > kmin) 
  { 
    # Find cluster pair to merge 
    found_pair = FALSE 
    
    while (!found_pair) 
    { 
      # Candidate: pair with minimum distance 
      minDist = Inf 
      
      for (i in seq(n_clust-1)) { 
        for (j in c((i+1):n_clust)) 
        { 
          if (distMatrix[i,j] < minDist) 
          { 
            minDist = distMatrix[i,j] 
            c1 = i 
            c2 = j 
          } 
        } 
      } 
      
      start.time = Sys.time() 
      # Compute which events belong to the candidate clusters 
      nodes1 = which(node_mapping == c1) 
      nodes2 = which(node_mapping == c2) 
      events1 = matrix(nrow = 0, ncol = length(cols)) 
      events2 = matrix(nrow = 0, ncol = length(cols)) 
      for (node in nodes1) 
      { 
        events1 = rbind(events1, events_in_nodes[[node]]) 
      } 
      for (node in nodes2) 
      { 
        events2 = rbind(events2, events_in_nodes[[node]]) 
      } 
      l1 = nrow(events1) 
      l2 = nrow(events2) 
      end.time = Sys.time() 
      
      
      if (verbose) { cat("Currently", n_clust, "clusters, min distance is", minDist, "\n") } 
      #if (verbose) { cat("Event membership took", end.time - start.time, "\n")} 
      
      # If the minimum distance is too large (e.g. infinite), exit algorithm 
      if (minDist >= dist_absolute_threshold)  
      {  
        done_algorithm = TRUE 
        break 
      } 
      
      # Perform dip test 
      if (l1 == 0 | l1 == 1 | l2 == 0 | l2 == 1 | perform_dip == FALSE) 
      { 
        # If one cluster is empty or a singleton, merge it without performing dip test 
        found_pair = TRUE 
      } 
      else 
      { 
        # If dip is smaller than threshold, accept the pair. 
        # Otherwise, set distance to infinity, so it's not considered again. 
        start.time = Sys.time() 
        dip_value = dip_test(events1, events2, centroids[c1,], centroids[c2,], l1, l2, 
                             lda = dip_lda, plot = FALSE) 
        end.time= Sys.time() 
        #if (verbose) { cat("Dip test took", end.time - start.time, "\n")} 
        
        if (dip_value < dip_threshold) 
        { 
          found_pair = TRUE 
          # print(which(node_mapping == c1 | node_mapping == c2)) 
        } 
        else 
        { 
          distMatrix[c1,c2] = Inf 
          if (verbose) { cat("Dip test failed!\n") } 
          print(which(node_mapping == c1 | node_mapping == c2)) 
        } 
      } 
      
      remove(events1) 
      remove(events2) 
      gc() 
    } 
    
    if (done_algorithm) { break } 
    
    # Perform the merger 
    
    # Merge the centroids and update distance matrix 
    newCentroid = l1 * centroids[c1,] + l2 * centroids[c2,] 
    newCentroid = newCentroid / (l1+l2) 
    centroids = centroids[-c(c1,c2),] 
    centroids = rbind(centroids, newCentroid, deparse.level = 0) 
    
    distMatrix = update_distance_matrix(centroids, distMatrix, c1, c2) 
    
    # Update the cluster membership of SOM nodes 
    for (node in seq(n_clust_original)) 
    { 
      clust = node_mapping[node] 
      
      if (clust < c1) { next } 
      
      if (clust == c1 | clust == c2)  
      {  
        node_mapping[node] = n_clust - 1  
        next 
      } 
      
      if (clust > c1 & clust < c2)  
      {  
        node_mapping[node] = clust-1  
        next 
      } 
      
      if (clust > c2) 
      { 
        node_mapping[node] = clust -2 
      } 
    } 
    
    n_clust = n_clust - 1 
    
  } 
  
  if (verbose) { cat("Done, obtained ", n_clust, "clusters.\n")} 
  
  s = 0 
  for (node in seq(n_clust_original)) 
  { 
    s = s+nrow(events_in_nodes[[node]]) 
  } 
  if(verbose) { cat("Fraction of events discarded:", 1-s/nrow(fsom$FlowSOM$data), "\n")} 
  
  factor(node_mapping) 
} 



############################### 
# Analyze the clusters obtained 
############################### 


# Make kde of global distribution for each parameter 
make_kdes_global = function(data, parameters) 
{ 
  kdes = list() 
  for (param in parameters) 
  { 
    kdes[[param]] = bkde(data[,param]) 
  } 
  
  kdes 
} 

# Overlay global kde with cluster kde, for each parameter 
plot_cluster_histograms = function(global_kdes, cluster = NULL, cluster2 = NULL,
                                   parameters,  
                                   weights = NULL, overlay_hist = TRUE) 
{ 
  i = 1
  while (i * (i+1) < length(parameters))
  {
    i = i+1
  }
  
  par(mfrow = c(i,i+1), mar = c(1,0,3,0) + 0.1) 
  
  for (parameter in parameters) 
  { 
    # Plot kernel density estimates for each of the parameters 
    #kde = bkde(data[,parameter]) 
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



################################################################################ 
# Modify Wade's get_samples function to read and prepare the lymphocytes dataset 
################################################################################ 


# retrieve and prepare files as specified by a list of files and a dbase 
get_samples_lymph = function(dbase, flist, compensate=TRUE, transform=TRUE,  
                             derail=TRUE, nice.names = TRUE, verbose=FALSE)  
{ 
  
  n = length(flist) 
  ff_list = list() 
  for (i in 1:n) { 
    
    fn = flist[i] 
    ff = read.FCS(fn) 
    
    if (verbose) { 
      cat("processing", i, "of", n, "...") 
      n_orig = nrow(ff) 
      cat("n_orig =", n_orig, "-") 
    } 
    
    if (compensate) {ff = autocomp(ff)} 
    
    if (derail) { 
      ff = Subset(ff, rectangleGate("FSC-A"=c(-Inf, 262142))) 
    } 
    if (transform) { 
      ff = doTransform(ff, cols = 2:7, method = 'linear') 
      ff = doTransform(ff, cols = 8:18, method = 'biexp') 
    } 
    if (nice.names) { 
      dnames = colnames(ff)    # detector names 
      
      names = parameters(ff)$desc 
      names[11] = "UNK" 
      names[c(1:7)] = c("Time", "FSC-A", "FSC-H", "FSC-W", "SSC-A", "SSC-H", "SSC-W") 
      #repair any ugly names... 
      #names[7] <- 'CD279 BB515' 
      #names[15] = 'DumpLD BV510' 
      scat_names = c("FSC-A", "FSC-H", "FSC-W", "SSC-A", "SSC-H", "SSC-W") 
      
      parameters(ff)$desc[2:7] = scat_names 
      
      # separate out the antibody from the fluorochrome, save the detector (original name) 
      # re-constitute, with just the antibody as the name, and the detector and fluorochrome 
      # as the desc 
      fidx = 8:18 
      #tmp = strsplit(x = names[fidx], split = "<", fixed = TRUE) 
      #ab = sapply(tmp, function(x){x[1]}) 
      #fl = sapply(tmp, function(x){x[2]}) 
      
      #colnames(ff) = c("Time", scat_names, ab) 
      #parameters(ff)$desc = c("Time", scat_names, paste(dnames[fidx], fl)) 
      colnames(ff) = c("Time", scat_names, names[fidx]) 
      parameters(ff)$desc = c("Time", scat_names, dnames[fidx]) 
    } 
    ff_list[[i]] = ff 
    if (verbose) { 
      n_derailed = nrow(ff) 
      cat(" n_derailed =", n_derailed) 
      cat(" ...done.\n") 
    } 
  } 
  fs = flowSet(ff_list) 
  
  #ff_list[[1]] 
} 




############################################### 
# Modifications to the FlowSOM algorithm 
############################################### 

# Skip the time-consuming and ineffective ensemble clustering step; 
# Modify the neighborhood function for the SOM grid,  
# to control the elasticity of the map. 
# ... 

# Wrapper which constructs the FlowSOM object 
FlowSOM_nometa = function (input, pattern = ".fcs", compensate = FALSE, spillover = NULL,  
                           transform = FALSE, toTransform = NULL,  
                           transformFunction = flowCore::logicleTransform(),  
                           scale = TRUE, scaled.center = TRUE, scaled.scale = TRUE,  
                           silent = TRUE, colsToUse, importance = NULL,  
                           seed = NULL, ...)  
{ 
  if (!is.null(seed)) { 
    set.seed(seed) 
  } 
  
  # Read FCS files; if necessary, compensate, transform, scale. 
  t <- system.time(fsom <- ReadInput(input, pattern = pattern,  
                                     compensate = compensate, spillover = spillover,  
                                     transform = transform, toTransform = toTransform,  
                                     transformFunction = transformFunction,  
                                     scale = scale, scaled.center = scaled.center,  
                                     scaled.scale = scaled.scale,  
                                     silent = silent)) 
  gc() 
  if (!silent)  
    message(t[3], "\n") 
  
  # Build the map object: learn the position of the nodes, assign datapoints to nodes 
  t <- system.time(fsom <- BuildSOM_elastic(fsom, colsToUse, silent = silent,  
                                            importance = importance, ...)) 
  gc() 
  if (!silent)  
    message(t[3], "\n") 
  
  # Build the minimal spanning tree 
  t <- system.time(fsom <- BuildMST(fsom, silent = silent)) 
  gc() 
  if (!silent)  
    message(t[3], "\n") 
  
  # Erased metaclustering step from here 
  
  list(FlowSOM = fsom) 
} 


# Wrapper for the map object 
BuildSOM_elastic = function (fsom, colsToUse = NULL, silent = FALSE, ...)  
{ 
  if (!"data" %in% names(fsom)) { 
    stop("Please run the ReadInput function first!") 
  } 
  if (!silent)  
    message("Building SOM\n") 
  if (is.null(colsToUse)) { 
    colsToUse <- seq_len(ncol(fsom$data)) 
  } 
  
  # Call the SOM function, which actually constructs the map 
  fsom$map <- SOM_elastic(fsom$data[, colsToUse], silent = silent,  
                          ...) 
  gc() 
  fsom$map$colsUsed <- colsToUse 
  #fsom$map$medianValues <- t(sapply(seq_len(fsom$map$nNodes),  
  #                                  function(i) { 
  #                                    apply(subset(fsom$data, fsom$map$mapping[, 1] ==  
  #                                                   i), 2, stats::median) 
  #                                  })) 
  #fsom$map$medianValues[is.nan(fsom$map$medianValues)] <- 0 
  #colnames(fsom$map$medianValues) <- colnames(fsom$data) 
  #fsom$map$sdValues <- t(sapply(seq_len(fsom$map$nNodes),  
  #                              function(i) { 
  #                                apply(subset(fsom$data, fsom$map$mapping[, 1] ==  
  #                                               i), 2, stats::sd) 
  #                              })) 
  #fsom$map$sdValues[is.nan(fsom$map$sdValues)] <- 0 
  #colnames(fsom$map$sdValues) <- colnames(fsom$data) 
  fsom 
} 



# Build map: modified neighborhood function. 
# Decreased the interaction between SOM nodes: no reason to expect that 
# data lies along a 2-dimensional submanifold. 

# Original: radius = stats::quantile(nhbrdist, 0.67) * c(1, 0), rlen= 10 
SOM_elastic = function (data, xdim = 10, ydim = 10, rlen = 10, mst = 1,  
                        alpha = c(0.05, 0.01),  
                        radius = stats::quantile(nhbrdist, 0.67) * c(1, 0), 
                        init = FALSE, distf = 2,  
                        silent = FALSE, codes = NULL, importance = NULL)  
{ 
  if (!is.null(importance)) { 
    data <- data * rep(importance, each = nrow(data)) 
  } 
  grid <- expand.grid(seq_len(xdim), seq_len(ydim)) 
  nCodes <- nrow(grid) 
  if (is.null(codes)) { 
    if (init) { 
      starters <- FlowSOM:::Initialize(data, nCodes) 
      message("Initialization ready\n") 
    } 
    else { 
      starters <- sample(1:nrow(data), nCodes, replace = FALSE) 
    } 
    codes <- data[starters, , drop = FALSE] 
  } 
  nhbrdist <- as.matrix(stats::dist(grid, method = "maximum")) 
  if (mst == 1) { 
    radius <- list(radius) 
    alpha <- list(alpha) 
  } 
  else { 
    radius <- seq(radius[1], radius[2], length.out = mst +  
                    1) 
    radius <- lapply(1:mst, function(i) { 
      c(radius[i], radius[i + 1]) 
    }) 
    alpha <- seq(alpha[1], alpha[2], length.out = mst +  
                   1) 
    alpha <- lapply(1:mst, function(i) { 
      c(alpha[i], alpha[i + 1]) 
    }) 
  } 
  
  
  # Changed radius[[i]] to 0 
  
  for (i in seq_len(mst)) { 
    res <- .C("C_SOM", data = as.double(data), codes = as.double(codes),  
              nhbrdist = as.double(nhbrdist), alpha = as.double(alpha[[i]]),  
              radius = as.double(radius[[i]]), xdists = double(nCodes),  
              n = as.integer(nrow(data)), px = as.integer(ncol(data)),  
              ncodes = as.integer(nCodes), rlen = as.integer(rlen),  
              distf = as.integer(distf)) 
    codes <- matrix(res$codes, nrow(codes), ncol(codes)) 
    colnames(codes) <- colnames(data) 
    nhbrdist <- FlowSOM:::Dist.MST(codes) 
  } 
  if (!silent)  
    message("Mapping data to SOM\n") 
  mapping <- FlowSOM:::MapDataToCodes(codes, data) 
  list(xdim = xdim, ydim = ydim, rlen = rlen, mst = mst, alpha = alpha,  
       radius = radius, init = init, distf = distf, grid = grid,  
       codes = codes, mapping = mapping, nNodes = nCodes) 
} 



####################################### 
# Map events to clusters 
####################################### 
# E.g. if downsampled before clustering 
# Complexity (O(n*k*d^2)), too slow 

events_to_clusters = function(centroids, events, verbose = FALSE) 
{ 
  n = nrow(events) 
  k = nrow(centroids) 
  mapping = integer(length = n) 
  distances = numeric(length = n) 
  
  for (i in seq(n)) 
  { 
    dmin = Inf 
    jmin = 0 
    
    for (j in seq(k)) 
    { 
      d = dist(rbind(events[i,],centroids[j,])) 
      if (d < dmin) 
      { 
        dmin = d 
        jmin = j 
      } 
    } 
    
    mapping[i] = jmin 
    distances[i] = dmin 
    
    if(verbose & i %%10000 == 0) { cat("Done", i , "\n")} 
  } 
  
  cbind(mapping, distances) 
} 




compute_somnode_info = function(fsom) 
{ 
  n_nodes = nrow(fsom$FlowSOM$map$codes) 
  
  somnode_info = data.frame(matrix(nrow = n_nodes, ncol = 4) ) 
  names(somnode_info) = c("size", "mean", "median", "fraction>4") 
  
  for (node in seq(n_nodes)) 
  { 
    if (node %% 50 == 0) {print(node)} 
    events_node = fsom$FlowSOM$map$mapping[which(fsom$FlowSOM$map$mapping[,1] == node),2] 
    l = length(events_node) 
    somnode_info[node,] = c(l, mean(events_node), median(events_node), 
                            length(which(events_node > 4))/l) 
  } 
  
  somnode_info 
} 


# Assuming that dat contains a boolean column which distinguishes the two arms 
kde_by_arm = function(dat, parameter) 
{ 
  Ipi = bkde(dat[which(dat[,"Nivo"] == FALSE),parameter]) 
  Nivo = bkde(dat[which(dat[,"Nivo"] == TRUE),parameter]) 
  
  plot(Ipi, type = "l", col = "red",  
       xlab = "expression", ylab = "density",  
       main = parameter) 
  lines(Nivo, col = "green") 
  legend("topright", legend = c("Ipi", "Ipi + Nivo"), col = c("red", "green"),lty=1:2, cex=0.8) 
} 


# From Winston Chang's R cookbook 
multiplot <- function(..., plotlist=NULL, cols) { 
  require(grid) 
  
  # Make a list from the ... arguments and plotlist 
  plots <- c(list(...), plotlist) 
  
  numPlots = length(plots) 
  
  # Make the panel 
  plotCols = cols                          # Number of columns of plots 
  plotRows = ceiling(numPlots/plotCols) # Number of rows needed, calculated from # of cols 
  
  # Set up the page 
  grid.newpage() 
  pushViewport(viewport(layout = grid.layout(plotRows, plotCols))) 
  vplayout <- function(x, y) 
    viewport(layout.pos.row = x, layout.pos.col = y) 
  
  # Make each plot, in the correct location 
  for (i in 1:numPlots) { 
    curRow = ceiling(i/plotCols) 
    curCol = (i-1) %% plotCols + 1 
    print(plots[[i]], vp = vplayout(curRow, curCol )) 
  } 
  
} 





######################################################################################### 
# Modified version of flowFP, because the original drops the last column, for some reason 
######################################################################################### 

flowFP_mod = function (fcs, model = NULL, sampleClasses = NULL, ...)  
{ 
  flowFP:::checkType(fcs, c("flowSet", "flowFrame"), "flowFP") 
  if (is.null(model)) { 
    model = flowFPModel_mod(fcs, ...) 
  } 
  fp = as(model, "flowFP") 
  fp@tags = list() 
  numFeatures = 2^fp@.cRecursions 
  if (is(fcs, "flowSet")) { 
    fp@counts = matrix(as.integer(0), nrow = length(fcs),  
                       ncol = numFeatures) 
    fp@sampleNames = sampleNames(fcs) 
    for (i in 1:length(fcs)) { 
      if (model@dequantize) { 
        fcs[[i]] = flowFP:::dequantize(fcs[[i]]) 
      } 
      fp@tags[[i]] = flowFP:::tag_events(fcs[[i]]@exprs, model) 
      fp@counts[i, ] = flowFP:::count_events(fp@tags[[i]], numFeatures) 
    } 
    if (!is.null(sampleClasses)) { 
      sampleClasses(fp) <- sampleClasses 
    } 
  } 
  else { 
    if (model@dequantize) { 
      fcs = flowFP:::dequantize(fcs) 
    } 
    fp@sampleNames = identifier(fcs) 
    fp@counts = matrix(as.integer(0), nrow = 1, ncol = numFeatures) 
    fp@tags[[1]] = flowFP:::tag_events(fcs@exprs, model) 
    fp@counts[1, ] = flowFP:::count_events(fp@tags[[1]], numFeatures) 
  } 
  return(fp) 
} 



flowFPModel_mod = function (fcs, name = "Default Model", parameters = NULL, nRecursions = "auto",  
                            dequantize = TRUE, sampleSize = NULL, excludeTime = TRUE)  
{ 
  flowFP:::checkType(fcs, c("flowSet", "flowFrame"), "flowFPModel") 
  if (dequantize)  
    fcs = flowFP:::dequantize(fcs) 
  if (is(fcs, "flowSet")) { 
    trainingSet = sampleNames(fcs) 
    if (is.null(sampleSize)) { 
      training_events = median(fsApply(fcs, nrow)) 
      fcs <- as(fcs, "flowFrame") 
      exprs(fcs) <- subset(exprs(fcs), subset = T, select = colnames(fcs)[1:(ncol(exprs(fcs)))]) 
    } 
    else { 
      nevents = vector(mode = "numeric", length = length(fcs)) 
      nevents[1] = min(sampleSize, nrow(fcs[[1]])) 
      tmp = fcs[[1]] 
      ind = sample(1:nrow(fcs[[1]]), nevents[1]) 
      exprs(tmp) = exprs(tmp)[ind, ] 
      if (length(fcs) > 1) { 
        for (i in 2:length(fcs)) { 
          nevents[i] = min(sampleSize, nrow(fcs[[i]])) 
          ind = sample(1:nrow(fcs[[i]]), nevents[i]) 
          exprs(tmp) = rbind(exprs(tmp), exprs(fcs[[i]])[ind,  
                                                         ]) 
        } 
      } 
      training_events = median(nevents) 
      fcs = tmp 
    } 
  } 
  else { 
    if (length(fcs) == 0) { 
      stop("Can not create a model from a flowFrame with zero events\n") 
    } 
    trainingSet = identifier(fcs) 
    if (is.null(sampleSize)) { 
      training_events = nrow(fcs) 
    } 
    else { 
      training_events = min(sampleSize, nrow(fcs)) 
      ind = sample(1:nrow(fcs), training_events) 
      exprs(fcs) = exprs(fcs)[ind, ] 
    } 
  } 
  parameters = flowFP:::parse_parameters(colnames(fcs), parameters,  
                                         excludeTime) 
  if (nRecursions == "auto") { 
    nRecursions = (log2(training_events)%/%1) - 3 
  } 
  flowFP:::validate_params(fcs, parameters, nRecursions) 
  model = new("flowFPModel", name = name, parameters = parameters,  
              nRecursions = nRecursions, .cRecursions = nRecursions,  
              trainingSet = trainingSet, dequantize = dequantize,  
              trainingSetParams = colnames(fcs)) 
  model@.tmp_tags = vector(mode = "integer", length = nrow(fcs)) 
  model@.tmp_tags[] = as.integer(1) 
  for (i in 1:model@nRecursions) { 
    model <- flowFP:::bin_level(fcs, model, i) 
  } 
  model@binBoundary = flowFP:::createBinBoundaries(model, fcs) 
  model@.tmp_tags = vector(mode = "integer", length = 0) 
  return(model) 
} 






############################ 
# Box information for flowFP 
############################ 



calculate_bin_phenotypes_fast = function(fp, fs, method=c("median", "mean"), 
                                         do_range = TRUE, verbose = FALSE) { 
  start.time = Sys.time() 
  parameters = parameters(fp) 
  n_bins = 2 ^ nRecursions(fp) 
  n_parameters = length(parameters) 
  center = matrix(NA, nrow = n_bins, ncol = n_parameters) 
  colnames(center) = parameters 
  range = matrix(NA, nrow = n_bins, ncol = n_parameters) 
  colnames(range) = parameters 
  
  max_count = max(fp@counts) 
  n_events = length(fp@tags[[1]]) 
  
  # for convenience, lump the frames in fs together and recalculate the fingerprint 
  ff = as(fs, "flowFrame") 
  #fp = flowFP(fcs = ff, model = as(fp, "flowFPModel")) 
  
  counts = integer(length = n_bins) 
  idx_list = vector("list", n_bins) 
  for (clas in c(1:n_bins)) 
  { 
    idx_list[[clas]] = integer(length= max_count) 
  } 
  
  for (event in c(1:n_events)) 
  { 
    clas = fp@tags[[1]][event] 
    new_count = counts[clas] + 1 
    counts[clas] = new_count 
    idx_list[[clas]][new_count] = event 
  } 
  
  for (i in 1:n_bins) { 
    idx = idx_list[[i]] 
    if (length(idx) == 0) { 
      next 
    } 
    
    for (j in 1:n_parameters) { 
      p = parameters[j] 
      vals = exprs(ff)[idx, p] 
      if (method == "median") { 
        center[i, j] = median(vals, na.rm = TRUE) 
        if (do_range) 
        { 
          range[i, j] = (quantile(x = vals, probs = 0.75, na.rm = TRUE) - 
                           quantile(x = vals, probs = 0.25, na.rm = TRUE)) / 2 
        } 
        
      } else { 
        
      } 
    } 
    
  } 
  if(verbose) {print(Sys.time() - start.time)} 
  return(list(center = center, range = range)) 
} 






####################################################### 
# Utils for nd mixture models with weighted subsampling 
####################################################### 


get_weighted_subsample_1d = function(data, mi, ma, max_count = 20000, nodes = 10, 
                                     step = 0.02, seed = 91, verbose = FALSE) 
{ 
  s = matrix(nrow = 0, ncol = 1) 
  colnames(s) = colnames(data) 
  w = c() 
  v = c() 
  pts = seq(from = mi, to = ma, by = step) 
  n = nrow(dat) 
  l = length(pts) 
  
  counts = integer(l) 
  idx_list = vector("list", l) 
  for (clas in c(1:l)) 
  { 
    idx_list[[clas]] = integer(length = max_count) 
  } 
  
  for (i in c(1:n)) 
  { 
    clas = ceiling((data[i] - mi) / step) 
    if(clas < 1) {clas = 1} 
    if(clas > l) {clas = l} 
    
    new_count = counts[clas] + 1 
    counts[clas] = new_count 
    idx_list[[clas]][new_count] = i 
  } 
  
  for (clas in c(1:l)) 
  { 
    count = counts[clas] 
    if (count == 0) { next } 
    sel = idx_list[[clas]][c(1:count)] 
    
    if (count <= nodes) 
    { 
      # If the bin is very sparse, take all events; there's no correction to variance 
      samp = data[sel] 
      weights = rep(1,count) 
      var = numeric(count) 
    } 
    else 
    { 
      # If the bin is dense, find 5-10 good representatives 
      set.seed(seed) 
      km = kmeans(data[sel], nodes, iter.max = 10, nstart = 1, algorithm = "MacQueen") 
      samp = km$centers 
      weights = km$size 
      
      # Find the correction to variance 
      var = numeric(nodes) 
      for (i in c(1:count)) 
      { 
        var[km$cluster[i]] = var[km$cluster[i]] + (data[sel[i]] - km$cluster[i])^2 
      } 
    } 
    
    s = rbind(s,samp) 
    w = c(w, weights) 
    v = c(v, var) 
  } 
  
  list(samp = s, weight = w, variance = v) 
} 







# Automatized sequential gating for basic functional phenotypes 
umap_add_phenotype = function(layout, dat) 
{ 
  # Find the manual gating strategy for Tregs 
  #cd4 = which(dat[,"CD4"] > 1.9) 
  #gate_treg = data.frame(dat[cd4,c("CD25", "CD127")]) 
  #ggplot(gate_treg, aes(x=CD25, y=CD127) ) + 
  #  geom_point(shape = 20, size = 1, alpha = 0.5) 
  
  # Gating 
  proj = data.frame(layout) 
  proj = cbind(proj, dat[,"CD4"] > 2) 
  proj = cbind(proj, dat[,"CD8"] > 2) 
  proj = cbind(proj, dat[,"CD45RA"] > 1.8) 
  proj = cbind(proj, dat[,"CD197"] > 1.8) 
  proj = cbind(proj, dat[,"CD25"] > 1.75) 
  proj = cbind(proj, dat[,"CD127"] > 1.1)  
  proj = cbind(proj, dat[,"Ki67"] > 1) 
  names(proj) = c("umap1", "umap2", "CD4", "CD8", "CD45RA","CD197", 
                  "CD25","CD127","Ki67") 
  
  # Phenotyping 
  phenotype = character(nrow(proj)) 
  phenotype[c(1:nrow(proj))] = "" 
  cd4 = which(proj[,"CD4"]) 
  cd8 = which(proj[,"CD8"]) 
  naive = which(proj[,"CD45RA"] & proj[,"CD197"]) 
  emra = which(proj[,"CD45RA"] & !proj[,"CD197"]) 
  cm = which(!proj[,"CD45RA"] & proj[,"CD197"]) 
  em = which(!proj[,"CD45RA"] & !proj[,"CD197"]) 
  treg = which(proj[,"CD4"] & proj[,"CD25"] & !proj[,"CD127"]) 
  
  phenotype[cd4] = sapply(phenotype[cd4], FUN = function(x) {paste(x,"CD4",sep = "")}, simplify = TRUE) 
  phenotype[cd8] = sapply(phenotype[cd8], FUN = function(x) {paste(x,"CD8",sep = "")}, simplify = TRUE) 
  phenotype[naive] = sapply(phenotype[naive], FUN = function(x) {paste(x,"naive",sep = " ")}, simplify = TRUE) 
  phenotype[emra] = sapply(phenotype[emra], FUN = function(x) {paste(x,"emra",sep = " ")}, simplify = TRUE) 
  phenotype[cm] = sapply(phenotype[cm], FUN = function(x) {paste(x,"cm",sep = " ")}, simplify = TRUE) 
  phenotype[em] = sapply(phenotype[em], FUN = function(x) {paste(x,"em",sep = " ")}, simplify = TRUE) 
  phenotype[which(proj[,"CD4"] & proj[,"CD8"])] = "DP" 
  phenotype[which(!proj[,"CD4"] & !proj[,"CD8"])] = "DN" 
  phenotype[treg] = "Treg" 
  #phenotype[ which(proj[,"Ki67"])] = "Ki67+" 
  
  cbind(proj, phenotype) 
} 





###################################### 
# Pregate CD4 / CD8 
###################################### 

gate_cd4_cd8 = function(ff, location, show = FALSE) { 
  params = c("CD4", "CD8") 
  # pregate 
  gate = rectangleGate("CD8"=c(1.5, Inf), "CD4"=c(-Inf, 2.2)) 
  tmp = Subset(ff, gate) 
  # find the blob 
  bb = blob.boundary(ff = tmp, parameters = params,  
                     location = location, height = .04)  # upped height from .2 to .3 5/1/2019 
  inflate_dist = .5 
  bb_infl = inflate.contour(get.hull(bb), dist = inflate_dist) 
  ff_gated = Subset(ff, polygonGate(.gate = bb_infl)) 
  
  # Gate the inflated contour again 
  gate = rectangleGate("CD8"=c(1.3, Inf), "CD4"=c(-Inf, 2.4)) 
  ff_gated = Subset(ff_gated, gate) 
  
  if(show) { 
    pplot(ff, params, xlim = c(-1, 5.4), ylim = c(-1, 5.4)) 
    lines(bb) 
    lines(bb_infl, col = 'red', lwd = 2) 
    xline(2.4) 
    yline(1.3) 
    annotate_in_box(xll = 10000, yll = 10000, xur = 200000, nrows = 1, label = sprintf("N = %d", nrow(ff_gated)), cex = 1.5) 
  } 
  
  ff_gated 
} 




################################ 
# Utils for statistical analysis 
################################ 

make_cluster_percentage_df = function(rbase, fs_name, outlier_threshold) 
{ 
  # Load the relevant data frame 
  load(tight(rbase, "fs_", fs_name, ".rda")) 
  counts = get(tight("fs_", fs_name))@phenoData@data[[tight("n_", fs_name)]] 
  
  # Load the cluster assignment, and compute cluster percentages 
  assignment_path = tight(rbase, "mapping_", fs_name, ".rda") 
  load(assignment_path) # mapping 
  pctg = cluster_percentages(counts, mapping, verbose = TRUE) 
  
  # Remove outlier clusters 
  outliers= which(tabulate(mapping) < outlier_threshold) 
  pctg = pctg[,-outliers] 
  
  # combine with clinical metadata for ease of analysis 
  pdat = get(tight("fs_", fs_name))@phenoData@data 
  data.frame(pctg, pid = pdat$pid, visit = pdat$visit, arm = pdat$arm) 
} 


# Fold change over all visits, to track just persistent changes, but with more stat power 
fold_change_all_visits = function(df, pv_threshold = 0.00005) 
{ 
  n_clust = ncol(df) - 3 
  clusters_list = colnames(df)[c(1:n_clust)] 
  
  pvs = data.frame(matrix(nrow = 5, ncol = n_clust)) 
  row.names(pvs) = levels(df$visit)[-c(1,7)] 
  names(pvs) = clusters_list 
  small_pv = list() 
  
  bl = "BL" 
  df_bl = df[which(df$visit == bl),] 
  current_visit = "all" 
  
  
  df_current = df[which(df$visit != "BL"),] 
  df_merged = merge(df_bl, df_current, by.x = c("pid", "arm"), by.y = c("pid","arm"), sort = FALSE) 
  
  fold_change = data.frame(matrix(nrow = nrow(df_merged), ncol = n_clust)) 
  names(fold_change) = clusters_list 
  
  for (cluster in clusters_list) 
  { 
    cy = paste(cluster,".y", sep="") 
    cx = paste(cluster,".x", sep="") 
    fold_change[,cluster] = df_merged[,cy] / df_merged[,cx] - 1 
  } 
  fold_change = cbind(df_merged[,"arm", drop = FALSE], fold_change) 
  
  fold_change[is.na(fold_change)] = 0 
  fold_change[,clusters_list][fold_change[,clusters_list] > 10] = 10  
  
  
  
  
  # No plots, just p-values 
  for (cluster in clusters_list) 
  { 
    fmla = formula(paste(cluster, "~ arm")) 
    pv = wilcox.test(fmla, data = fold_change, exact = FALSE)$p.value 
    
    pvs[current_visit, cluster] = pv 
  } 
  
  small_pv[[current_visit]] = which(pvs[current_visit,] < pv_threshold) 
  #small_pv[[current_visit]] = small_pv[[current_visit]][!is.na(small_pv[[current_visit]])] 
  
  #fold_change[,names(small_pv[[current_visit]])] 
  
  
  
  #signif_cl = c() 
  #for (visit in levels(df$visit)) 
  #{ 
  #  if (visit == "BL") {next} 
  #   
  #  signif_cl = c(signif_cl, small_pv[[visit]]) 
  #} 
  #signif_cl = unique(signif_cl) 
  
  #signif_changes = data.frame(matrix(nrow = length(small_pv[[current_visit]]), ncol = 2)) 
  #names(signif_changes) = c("cluster", "pv") 
  
  #signif_changes$cluster = small_pv[[current_visit]] 
  #signif_changes$pv = pvs[current_visit, small_pv[[current_visit]] ] 
  
  signif_changes = which(pvs[current_visit,] < pv_threshold) 
  pvs[current_visit, signif_changes] 
  
  
  #for (visit in names(signif_changes)) 
  #{ 
  #  for (cl in row.names(signif_changes)) 
  #  { 
  #    signif_changes[cl,visit] = pvs[visit,cl] 
  #  } 
  #} 
  
} 


# Given a dataframe of cluster percentages and patient metadata, compute fold change over BL 
# (More accurately, log of relative size: log(size_visit / size_BL) ) 
compute_fold_change = function(df) 
{ 
  n_clust = ncol(df) - 3 
  clusters_list = colnames(df)[c(1:n_clust)] 
  
  df_bl = df[which(df$visit == "BL"),] 
  df_visits = df[which(df$visit != "BL"),] 
  
  # Merge baseline df with visit df on patient ID and arm 
  df_merged = merge(df_bl, df_visits,  
                    by.x = c("pid", "arm"),  
                    by.y = c("pid","arm"),  
                    sort = FALSE) 
  n_data = nrow(df_merged) 
  
  fold_change = data.frame(matrix(nrow = n_data, ncol = n_clust)) 
  names(fold_change) = clusters_list 
  
  # Compute raw fold change (more accurately, relative size) 
  for (cluster in clusters_list) 
  { 
    cy = paste(cluster,".y", sep="") 
    cx = paste(cluster,".x", sep="") 
    fold_change[,cluster] = df_merged[,cy] / df_merged[,cx] 
  } 
  
  # Clean NaN, Inf, 0 and apply log 
  fold_change[is.na(fold_change)] = 1 # NaN only occurs from 0/0, which means no change 
  fold_change[fold_change < 1e-3] = 1e-3 
  fold_change[fold_change > 1e+3] = 1e+3 
  fold_change = log(fold_change, base = 2) 
  
  # Merge with metadata 
  fold_change = cbind(df_merged[,c("arm","pid","visit.y")], fold_change) 
  
  fold_change 
} 


# Model increase or decrease over baseline as binomial trial 
fold_change_binomial = function(fold_change, visit = "all", pv_threshold = 5e-5) 
{ 
  if (visit != "all") 
  { 
    fold_change = fold_change[fold_change$visit.y == visit ,] 
  } 
  
  n_clust = ncol(fold_change) - 3 
  n_data  = nrow(fold_change) 
  
  pvs = numeric(n_clust) 
  n_increased = integer(n_clust) 
  
  # For each cluster, measure increases vs decreases over baseline, and apply binomial test 
  for (i in c(1:n_clust)) 
  { 
    cluster = clusters_list[i] 
    n_increased[i] = length(which(fold_change[,cluster] > 0)) 
    pvs[i] = binom.test(x = n_increased[i], n = n_data, alternative = "two.sided")$p.value 
  } 
  
  # For which clusters is the test significant? 
  signif = which(pvs < pv_threshold) 
  
  if (length(signif) == 0) 
  { 
    print("No cluster is statistically significant; try increasing threshold.") 
    return(NULL) 
  } 
  
  signif_data = data.frame(matrix(nrow = length(signif), ncol = 3)) 
  names(signif_data) = c("cluster", "pv", "change") 
  signif_data$cluster = clusters_list[signif] 
  signif_data$pv = signif(pvs[signif], 3) 
  signif_data$change = 1 
  signif_data$change[n_increased[signif] < n_data / 2] = -1 
  
  signif_data 
} 



# TO DO: finish this 
plot_fold_change_binomial = function(fold_change, visit = "all", clusters) 
{ 
  if (visit != "all") 
  { 
    fold_change = fold_change[fold_change$visit.y == visit ,] 
  } 
  
  n_data  = nrow(fold_change) 
  
  for (cl in clusters) 
  { 
    ggplot(data = fold_change, aes(x = fold_change$arm)) 
  } 
} 



get_1D_cutoffs = function(mixtures, to_merge, params)
{
  cutoffs = list()
  
  for (param in params)
  {
    mix = mixtures[[param]]
    if (length(mix$pro) == 1)
    {
      cutoffs[[param]] = NULL
      next
    }
    
    i1 = 1
    i2 = 2
    
    if (length(mix$pro) == 3)
    {
      i1 = 2
      i2 = 3
    }
    
    if (!is.null(to_merge[[param]]))
    {
      if (to_merge[[param]][1] == 2)
      {
        i1 = 1
        i2 = 2
      }
    }
    
    
    
    pro1 = mix$pro[i1]
    pro2 = mix$pro[i2]
    mean1 = mix$mean[i1]
    mean2 = mix$mean[i2]
    sd1 = sqrt(mix$variance$sigmasq[i1])
    sd2 = sqrt(mix$variance$sigmasq[i2])
    
    # look for point where second gaussian becomes larger than the first
    pos = mean1
    while (pro1 * dnorm(pos, mean = mean1, sd = sd1) > pro2 * dnorm(pos, mean = mean2, sd = sd2))
    {
      pos = pos + 0.05
    }
    
    cutoffs[[param]] = pos
  }
  
  cutoffs
}

print_kdes_with_cutoffs = function(kdes, cutoffs, parameters)
{
  i = 1
  while (i * (i+1) < length(parameters))
  {
    i = i+1
  }
  
  par(mfrow = c(i,i+1), mar = c(1,0,3,0) + 0.1) 
  
  for (parameter in parameters) 
  { 
    kde = global_kdes[[parameter]]
    plot(kde, type = "l", main = parameter)
    if(!is.null(cutoffs[[parameter]]))
    {
      abline(v = cutoffs[[parameter]])
    }
  }
}



categorical_merging = function(pros, means, cutoffs, params)
{
  n = nrow(means)
  categs = vector(mode = "character", length = n)
  cat_mapping = vector(mode = "integer", length = n)
  
  for (param in params)
  {
    if (is.null(cutoffs[[param]]))
    {
      for (i in seq(n))
      {
        categs[i] = paste(categs[i], "-", sep = "")
      }
      next
    }
    
    v = means[,param] - cutoffs[[param]]
    for (i in seq(n))
    {
      if (v[i]<0)
      {
        categs[i] = paste(categs[i], "-", sep = "")
      }
      else
      {
        categs[i] = paste(categs[i], "+", sep = "")
      }
    }
  }
  
  un = unique(categs)
  for (i in seq(n))
  {
    cat_mapping[i] = which(un == categs[i])
  }
  
  centers = matrix(0, nrow = length(un), ncol = length(params))
  colnames(centers) = params
  for (i in seq(length(un)))
  {
    weight = 0
    components = which(cat_mapping == i)
    for (comp in components)
    {
      centers[i,] = centers[i,] +  pros[comp] * means[comp,]
      weight = weight + pros[comp]
    }
    centers[i,] = centers[i,] / weight
  }
  
  list(categs = un, cat_mapping = cat_mapping, centers = centers)
}





categorical_labeling = function(cats, defs, params)
{
  n = length(cats$categs)
  cats[["labels"]] = vector(mode = "character", length = n)
  labs = rownames(defs)
  nam = names(defs)
  ind = vector(mode = "integer", length = length(nam))
  for (i in seq(length(nam)))
  {
    ind[i] = which(params == nam[i])
  }
  
  for (i in seq(n))
  {
    cats$labels[i] = "UNK"
    
    for (j in seq(nrow(defs)))
    {
      match = TRUE
      for (k in seq(ncol(defs)))
      {
        if (defs[j,k] == "hi" & substr(cats$categs[i],ind[k],ind[k]) == "-" | 
            defs[j,k] == "lo" & substr(cats$categs[i],ind[k],ind[k]) == "+")
        {
          match = FALSE
        }
      }
      if (match)
      {
        cats$labels[i] = labs[j]
        break
      }
    }
  }
  
  cats
}


tsne_categorical_clusters = function(centers, seed = 137)
{
  n_items = nrow(centers)
  perplexity = min((n_items - 1)/3, 30)
  set.seed(seed)
  res = Rtsne(dist(centers), perplexity = perplexity)$Y
  colnames(res) = c("tsne_1", "tsne_2")
  
  data.frame(res)
}



get_tsne_centers = function(data, seed, probs)
{
  n_items = nrow(data)
  perplexity = min((n_items - 1)/3, 30)
  set.seed(seed)
  
  res = Rtsne(data, perplexity = perplexity)$Y
  cluster = apply(probs, 1, which.max)
  res = cbind(res, cluster)
  
  colnames(res) = c("tsne_1", "tsne_2", "cluster")
  
  data.frame(res)
}



