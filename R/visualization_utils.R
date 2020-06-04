parallel_pheno_mod = function(mfi, idx_bin = 1:length(mfi[[1]]), parameters = names(mfi), col = 'red', mfi_colors = FALSE, show_contours = TRUE, bars = TRUE, ...) {
  
  if (is.numeric(parameters)) {
    pnames = colnames(ff)[parameters]
  } else {
    pnames = parameters
  }
  if (show_contours) {
    cont_col = 'darkgray'
  } else {
    cont_col = 'white'
  }
  #mfi = data.frame(mfi)
  # hyphens in parameter names get turned into dots here.  Change back
  colnames(mfi) = sub(pattern = ".", replacement = "-", x = colnames(mfi), fixed = TRUE)
  
  
  parallel_contours_mod(mfi, parameters = parameters, col = cont_col, ...)
  mfi = mfi[idx_bin, parameters]
  
  med_vec = vector(mode = 'numeric')
  q1_vec = vector(mode = 'numeric')
  q3_vec = vector(mode = 'numeric')
  for (i in 1:length(parameters)) {
    tmp = fivenum(mfi[, i])
    med_vec[i] = tmp[3]
    q1_vec[i] = tmp[2]
    q3_vec[i] = tmp[4]
  }
  # draw the median
  if (bars) {
    if (mfi_colors) {
      col = pcolor(med_vec, min_value = 0, max_value = 5)
    }
    add_bars(vals = med_vec, yvals = 1:length(parameters), col = col)
  } else {
    x = med_vec
    y = 1:length(parameters)
    lines(x, y, col = col, lwd = 3)
  }
  
  # draw the flags
  for (i in 1:length(parameters)) {
    if (bars) {
      draw_flag(y = i, q1 = q1_vec[i], q3 = q3_vec[i], med = NA, cex = 2, lwd = 2)
    } else {
      draw_flag(y = i, q1 = q1_vec[i], q3 = q3_vec[i], med = med_vec[i], cex = 2, lwd = 2)
    }
  }
}

parallel_contours_mod = function(mfi, parameters = names(mfi), col = 'blue', ...) {
  if (is.numeric(parameters)) {
    pnames = colnames(ff)[parameters]
  } else {
    pnames = parameters
  }
  plot(0, 0, pch = '', xlim = c(0, bx(262143)), ylim = c(1 - .3, length(parameters) + .3),
       xaxt = 'n', yaxt = 'n',
       xlab = '', ylab = '')
  ax(1, instrument = 'diva', type = 'biexp')
  axis(side = 2, labels = pnames, at = 1:length(pnames), las = 1, ...)

  x1 = matrix(t(mfi), nrow(mfi) * length(parameters))
  x2 = rep(c(1:length(parameters)), nrow(mfi))
  x = cbind(x1, x2)
  #x = matrix(NA, nrow = nrow(mfi) * length(parameters), ncol = 2)
  #k = 1
  #for (i in 1:nrow(mfi)) {
  #  if (i %% 10000 == 0) { print (i)}
  #  for (p in 1:length(parameters)) {
  #    x[k, 1] = mfi[i, p]
  #    x[k, 2] = p
  #    k = k + 1
  #  }
  #}

  kde = bkde2D(x = x, bandwidth = c(.1, 1), gridsize = c(501, 501))
  kde$fhat = kde$fhat / max(kde$fhat)
  contour(x = kde$x1, y = kde$x2, z = kde$fhat, col = col,
          xaxt = 'n', yaxt = 'n',
          drawlabels = FALSE,
          # levels = seq(.01, .2, length.out = 20),
          add = TRUE)
  
}


dotplot_by_arm = function(dat1, dat2, lab1, lab2, title)
{
  df = data.frame(  rbind( cbind(dat1, rep(lab1, length(dat1))) , 
                           cbind(dat2, rep(lab2, length(dat2))) ) )
  names(df) = c("cluster_pctg", "arm")
  
  ggplot(df, aes(x=arm, y=cluster_pctg)) + 
    geom_dotplot(binaxis='y', stackdir='center') +
    labs(title = title)
}

