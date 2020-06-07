

####################################### 
# Map events to clusters 
####################################### 
# E.g. if downsampled before clustering 
# Complexity (O(n*k*d^2)), too slow 

# events_to_clusters = function(centroids, events, verbose = FALSE) 
# { 
#   n = nrow(events) 
#   k = nrow(centroids) 
#   mapping = integer(length = n) 
#   distances = numeric(length = n) 
#   
#   for (i in seq(n)) 
#   { 
#     dmin = Inf 
#     jmin = 0 
#     
#     for (j in seq(k)) 
#     { 
#       d = dist(rbind(events[i,],centroids[j,])) 
#       if (d < dmin) 
#       { 
#         dmin = d 
#         jmin = j 
#       } 
#     } 
#     
#     mapping[i] = jmin 
#     distances[i] = dmin 
#     
#     if(verbose & i %%10000 == 0) { cat("Done", i , "\n")} 
#   } 
#   
#   cbind(mapping, distances) 
# } 
# 






# ####################################################### 
# # Utils for nd mixture models with weighted subsampling 
# ####################################################### 
# 
# 
# get_weighted_subsample_1d = function(data, mi, ma, max_count = 20000, nodes = 10, 
#                                      step = 0.02, seed = 91, verbose = FALSE) 
# { 
#   s = matrix(nrow = 0, ncol = 1) 
#   colnames(s) = colnames(data) 
#   w = c() 
#   v = c() 
#   pts = seq(from = mi, to = ma, by = step) 
#   n = nrow(dat) 
#   l = length(pts) 
#   
#   counts = integer(l) 
#   idx_list = vector("list", l) 
#   for (clas in c(1:l)) 
#   { 
#     idx_list[[clas]] = integer(length = max_count) 
#   } 
#   
#   for (i in c(1:n)) 
#   { 
#     clas = ceiling((data[i] - mi) / step) 
#     if(clas < 1) {clas = 1} 
#     if(clas > l) {clas = l} 
#     
#     new_count = counts[clas] + 1 
#     counts[clas] = new_count 
#     idx_list[[clas]][new_count] = i 
#   } 
#   
#   for (clas in c(1:l)) 
#   { 
#     count = counts[clas] 
#     if (count == 0) { next } 
#     sel = idx_list[[clas]][c(1:count)] 
#     
#     if (count <= nodes) 
#     { 
#       # If the bin is very sparse, take all events; there's no correction to variance 
#       samp = data[sel] 
#       weights = rep(1,count) 
#       var = numeric(count) 
#     } 
#     else 
#     { 
#       # If the bin is dense, find 5-10 good representatives 
#       set.seed(seed) 
#       km = kmeans(data[sel], nodes, iter.max = 10, nstart = 1, algorithm = "MacQueen") 
#       samp = km$centers 
#       weights = km$size 
#       
#       # Find the correction to variance 
#       var = numeric(nodes) 
#       for (i in c(1:count)) 
#       { 
#         var[km$cluster[i]] = var[km$cluster[i]] + (data[sel[i]] - km$cluster[i])^2 
#       } 
#     } 
#     
#     s = rbind(s,samp) 
#     w = c(w, weights) 
#     v = c(v, var) 
#   } 
#   
#   list(samp = s, weight = w, variance = v) 
# } 
# 
# 


#################################################
# Utilities for categorical merging:
# merge mixture components which are unlikely to
# represent biologically distinct populations
#################################################


get_1D_cutoffs = function(mixtures, to_merge, params)
{
  ##############################################
  # Take as input a list of 1D mixture models
  # Output a list of cutoffs for each parameter:
  # a cutoff is the threshold between distinct 
  # Gaussian mixture components
  ##############################################
  
  cutoffs = list()
  
  # Ignore distinction between mixture components
  # specified in to_merge; we don't expect these
  # to be biologically different
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
    while (pro1 * dnorm(pos, mean = mean1, sd = sd1) 
           > pro2 * dnorm(pos, mean = mean2, sd = sd2))
    {
      pos = pos + 0.05
    }
    
    cutoffs[[param]] = pos
  }
  
  cutoffs
}


print_kdes_with_cutoffs = function(kdes, cutoffs, parameters)
{
  ###################################################
  # print kdes with vertical lines displaying cutoffs
  ###################################################
  
  # determine size of the grid of plots
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
    
    # We model some markers as unimodal: no cutoff
    if(!is.null(cutoffs[[parameter]]))
    {
      abline(v = cutoffs[[parameter]])
    }
  }
}



categorical_merging = function(pros, means, cutoffs, params)
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
  #################################################
  # Take input the list of categorical clusters
  # and definitions of major phenotypes (e.g. naive
  # CD4 cells). Label each cluster by one major
  # phenotype, or unknown.
  #################################################
  
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




