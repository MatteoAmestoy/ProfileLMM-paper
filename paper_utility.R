library(patchwork)



sim_lif = function(nInd,rInd,paramLat){
  n = nInd*rInd
  df = data.frame(Zind = rep(1:nInd,rInd))
  df$Zind = df$Zind-1
  # Time invariant individual caracteristics
  df[c('FEgender','FE1')] = t(as.matrix(as(factor(df$Zind),Class = "sparseMatrix")))%*%
    cbind(sample(c(0,1),nInd,replace = T),rmvn(nInd,rep(0,1),diag(1)))

  df[c('FE2','FE3')] = rmvn(n,rep(0,2),diag(2))


  theta0 = {} # parameters to estimate
  if(paramLat$name == 'Correct'){
    rho = paramLat$rho
    spread = paramLat$spread
    theta0$Lat_split = cbind(sample(c(-1,0,1),n,replace = T),
                             sample(c(-1,0,1),n,replace = T))

    theta0$Lat = as.factor(paste(theta0$Lat_split[,1],paste(theta0$Lat_split[,2],sep=''),sep=''))
    theta0$label=    levels(theta0$Lat)
    levels(theta0$Lat) = 1:9

    df[c('Pol1','Pol2','Soc1','Soc2')] = cbind(rmvn(n,theta0$Lat_split[,1]%*%t(rep(1,2)),spread*((1-rho)*diag(2)+rho)),
                                               rmvn(n,theta0$Lat_split[,2]%*%t(rep(1,2)),spread*((1-rho)*diag(2)+rho)))

  }else {
    if(paramLat$name == 'Flat'){
      df[c('Pol1','Soc1')] = matrix(runif(2*n,min = -2,max=2),nrow=n,ncol=2)
      df[c('Pol2','Soc2')] = df[c('Pol1','Soc1')] +rmvn(n,c(0,0),0.1*diag(2))
    }else if(paramLat$name == 'Wrong'){
      rho = paramLat$rho
      spread = paramLat$spread
      samples = cbind(c(0.65,0.65),
                      c(-0.65,-0.65),
                      c(0.65,-0.65),
                      c(-0.65,0.65),
                      c(0,0))
      latSamples = t(samples[,sample(1:5,n,replace = T)])

      df[c('Pol1','Pol2','Soc1','Soc2')] = cbind(rmvn(n,latSamples[,1]%*%t(rep(1,2)),spread*((1-rho)*diag(2)+rho)),
                                                 rmvn(n,latSamples[,2]%*%t(rep(1,2)),spread*((1-rho)*diag(2)+rho)))
    }

    score = matrix(0,nrow=n,ncol = 9)
    idx_ = 1
    latName = c()
    for (pol in c(-1,0,1)){
      for (soc in c(-1,0,1)){
        latName = c(latName,paste(pol,soc))
        score[,idx_] = exp(2*(-(df$Pol1-pol)**2-(df$Soc1-soc)**2))
        idx_ = idx_ +1
      }
    }
    print(latName)
    # score = score/rowSums(score)
    theta0$Lat = c()
    for (idx_ in 1:n){
      theta0$Lat = c(theta0$Lat,sample(latName,1,prob = score[idx_,]))
    }
    # print(score[1:5,])
    # print(theta0$Lat[1:5])
    theta0$Lat = as.factor(theta0$Lat)
    theta0$label=    levels(theta0$Lat)
    levels(theta0$Lat) = 1:9


  }

  # print(theta0$Lat[1:5])
  #   print(theta0$label[1:5])
  #

  # Time
  Ti = c()
  for (t in 1:rInd){
    Ti= c(Ti,runif(nInd,t-1,t))
  }
  df$Ti = Ti

  Bt = bs(df$Ti,df = 2,degree = 2)
  attributes(Bt) <- attributes(Bt)["dim"]
  df[c('RE1','RE2')] = Bt
  rm(t,Ti,Bt)


  REnames = c('RE1','RE2')
  FEnames = c('FEgender','FE1','FE2','FE3')
  Latnames =  c('FE1','FE2')
  Assignnames = c('Pol1','Pol2','Soc1','Soc2')

  theta0$betaFE = rnorm(length(FEnames)+1) #adding intercept

  theta0$WLat = 2*(0.9*diag(length(Latnames)+1)+0.1)

  for(cidx in 1:9){
    theta0$alphaLat[(1:(length(Latnames)+1))+(length(Latnames)+1)*(cidx-1)] = rmvn(1,rep(0,length(Latnames)+1),theta0$WLat)
  }
  bbb = (0:8)*(length(Latnames)+1)
  for(cidx in 1:(length(Latnames)+1)){
    theta0$alphaLat[bbb+cidx] = theta0$alphaLat[bbb+cidx]-mean(theta0$alphaLat[bbb+cidx])
  }


  theta0$SigRE = diag(length(REnames))+0.2
  theta0$alphaRE = c()
  for(ind in 1:nInd){
    theta0$alphaRE = c(theta0$alphaRE,rmvn(1,rep(0,length(REnames)),theta0$SigRE))
  }

  theta0$sigma = 1

  df$Y = drop(as.matrix(cbind(1,df[FEnames]))%*%theta0$betaFE +
                t(KhatriRao(t(t(as(factor(theta0$Lat),Class = "sparseMatrix"))),t(cbind(1,df[Latnames])),make.dimnames = TRUE))%*%theta0$alphaLat+
                t(KhatriRao(t(t(as(factor(df$Zind),Class = "sparseMatrix"))),t(df[REnames]),make.dimnames = TRUE))%*%theta0$alphaRE+
                rnorm(n,0,theta0$sigma))

  df$FEgender = factor(df$FEgender)


  covList={}
  covList$FE = FEnames
  covList$RE = REnames
  covList$Lat = Latnames
  covList$Assign = Assignnames
  covList$REunit = 'Zind'
  covList$Y = 'Y'

  return(list(df = df,covList = covList,theta0=theta0))

}

sim_lif_v2 = function(nInd,rInd,paramLat,theta0={}){
  n = nInd*rInd
  df = data.frame(Zind = rep(1:nInd,rInd))
  df$Zind = df$Zind-1
  # Time invariant individual caracteristics
  df[c('FEgender','FE1')] = t(as.matrix(as(factor(df$Zind),Class = "sparseMatrix")))%*%
    cbind(sample(c(0,1),nInd,replace = T),rmvn(nInd,rep(0,1),diag(1)))

  df[c('FE2','FE3')] = rmvn(n,rep(0,2),diag(2))

  REnames = c('RE1','RE2')
  FEnames = c('FEgender','FE1','FE2','FE3')
  Latnames =  c('FE2')
  Assignnames = c('Exp1','Exp2')

  if (length(names(theta0))==0){
    theta0 = {} # parameters to estimate
    theta0$betaFE = 2*rnorm(length(FEnames)+1) #adding intercept
    theta0$WLat = 4*diag(length(Latnames)+1)
    theta0$SigRE = diag(length(REnames))+0.2
    theta0$sigma = 0.5

    for(cidx in 1:9){
      theta0$alphaLat[(1:(length(Latnames)+1))+(length(Latnames)+1)*(cidx-1)] = rmvn(1,rep(0,length(Latnames)+1),theta0$WLat)
    }
    bbb = (0:8)*(length(Latnames)+1)
    for(cidx in 1:(length(Latnames)+1)){
      theta0$alphaLat[bbb+cidx] = theta0$alphaLat[bbb+cidx]-mean(theta0$alphaLat[bbb+cidx])
    }
  }

  theta0$alphaRE = c()
  for(ind in 1:nInd){
    theta0$alphaRE = c(theta0$alphaRE,rmvn(1,rep(0,length(REnames)),theta0$SigRE))
  }

  if(paramLat$name == 'Correct'|paramLat$name == 'Simple'){
    Lat_split = cbind(sample(c(-1,0,1),n,replace = T),
                      sample(c(-1,0,1),n,replace = T))

    df[c('Exp1','Exp2')] = Lat_split + rmvn(n,c(0,0),paramLat$spread*((1-paramLat$rho)*diag(2)+paramLat$rho))

  }else if(paramLat$name == 'Flat'){
    df[c('Exp1','Exp2')] = matrix(runif(2*n,min = -2,max=2),nrow=n,ncol=2)

  }else if(paramLat$name == 'Wrong'){
    df[c('Exp1','Exp2')] = rmvn(n,c(0,0),paramLat$spread*((1-paramLat$rho)*diag(2)+paramLat$rho))
  }

  score = matrix(0,nrow=n,ncol = 9)
  idx_ = 1
  latName = c()
  precMatClus  =  inv(paramLat$spreadClus*((1-paramLat$rhoClus)*diag(2)+paramLat$rhoClus))
  for (pol in c(-1,0,1)){
    for (soc in c(-1,0,1)){
      latName = c(latName,paste(pol,soc))
      mu = as.matrix(c(pol,soc))
      score[,idx_] = -diag(as.matrix(df[,c('Exp1','Exp2')])%*%
                            precMatClus%*%
                            t(as.matrix(df[,c('Exp1','Exp2')])))+2*as.matrix(df[,c('Exp1','Exp2')])%*%
        precMatClus%*%mu-rep(t(mu)%*%precMatClus%*%mu,n)
      idx_ = idx_ +1
    }
  }

  # score = score/rowSums(score)
  theta0$Lat = c()
  for (idx_ in 1:n){
    p_ = exp(score[idx_,]-max(score[idx_,]))/sum(exp(score[idx_,]-max(score[idx_,])))
    theta0$Lat = c(theta0$Lat,sample(latName,1,prob = p_))
  }
  theta0$Lat = as.factor(theta0$Lat)
  theta0$label=    levels(theta0$Lat)
  levels(theta0$Lat) = 1:9

  # Time
  Ti = c()
  for (t in 1:rInd){
    Ti= c(Ti,runif(nInd,t-1,t))
  }
  df$Ti = Ti

  Bt = bs(df$Ti,df = 2,degree = 2)
  attributes(Bt) <- attributes(Bt)["dim"]
  df[c('RE1','RE2')] = Bt
  rm(t,Ti,Bt)




  df$Y = drop(as.matrix(cbind(1,df[FEnames]))%*%theta0$betaFE +
                t(KhatriRao(t(t(as(factor(theta0$Lat),Class = "sparseMatrix"))),t(cbind(1,df[Latnames])),make.dimnames = TRUE))%*%theta0$alphaLat+
                t(KhatriRao(t(t(as(factor(df$Zind),Class = "sparseMatrix"))),t(df[REnames]),make.dimnames = TRUE))%*%theta0$alphaRE+
                rnorm(n,0,theta0$sigma))

  df$FEgender = factor(df$FEgender)


  covList={}
  covList$FE = FEnames
  covList$RE = REnames
  covList$Lat = Latnames
  covList$Assign = Assignnames
  covList$REunit = 'Zind'
  covList$Y = 'Y'

  return(list(df = df,covList = covList,theta0=theta0))

}



bins = function(vec,n){seq(min(vec),max(vec), (max(vec)-min(vec))/n) }

ellipse = function(theta,Sigma,r,n){
  n_points <- n
  xy <- cbind(sin(seq(0, 2 * pi, length.out = n_points)),
              cos(seq(0, 2 * pi, length.out = n_points)))

  # then we scale the dimensions
  ev <- eigen(Sigma)
  xy[, 1] <- xy[, 1] * 1
  xy[, 2] <- xy[, 2] * sqrt(min(ev$values) / max(ev$values))

  # then rotate
  phi <- atan(ev$vectors[2, 1] / ev$vectors[1, 1])
  R <- matrix(c(cos(phi), sin(phi), -sin(phi), cos(phi)), 2)
  xy <- tcrossprod(R, xy)


  return(list(xx = sqrt(r) * xy[1, ] + theta[1], yy=sqrt(r) * xy[2, ] + theta[2]))
}


ClusterPurity <- function(clusters, classes) {
  sum(apply(table(classes, clusters), 2, max)) / length(clusters)
}

create_co_occurrence_matrix_r <- function(clustering_labels_matrix) {
 
  num_observations <- nrow(clustering_labels_matrix)
  num_clustering_runs <- ncol(clustering_labels_matrix)

  # Initialize an N x N matrix with zeros to store co-occurrence counts
  co_occurrence_matrix <- matrix(0L, nrow = num_observations, ncol = num_observations)

  # Iterate through each clustering run (column)
  for (k in 1:num_clustering_runs) {
    # Get the cluster labels for the current run
    current_labels <- clustering_labels_matrix[, k]

    # Find unique cluster IDs in this run
    unique_cluster_ids <- unique(current_labels)

    # For each unique cluster, identify its members and update the co-occurrence matrix
    for (cluster_id in unique_cluster_ids) {
      # Find all observation indices that belong to the current cluster_id
      obs_indices <- which(current_labels == cluster_id)

      # If there is more than one observation in this cluster,
      # increment the co-occurrence counts for all pairs within this cluster.
      # This is the key vectorized operation: it adds 1 to all cells
      # corresponding to pairs within the 'obs_indices' subset of the matrix.
      if (length(obs_indices) > 1) {
        co_occurrence_matrix[obs_indices, obs_indices] <- co_occurrence_matrix[obs_indices, obs_indices] + 1
      } else {
        # Even if a cluster has only one member, that member is still clustered with itself.
        # This ensures diagonals are correctly incremented for all runs.
        co_occurrence_matrix[obs_indices, obs_indices] <- co_occurrence_matrix[obs_indices, obs_indices] + 1
      }
    }
  }

  return(co_occurrence_matrix)
}

plot_single_gaussian_density_zones <- function(means, covariances, x_range, y_range,
                                               grid_res, contour_levels,
                                               title, show_points,
                                               data_points, point_color,component_colors = NULL,
                                               x_axis_label = "X-axis", y_axis_label = "Y-axis") { # Added x_axis_label, y_axis_label
  # --- Input Validation (simplified for internal use, main function handles broad validation) ---
  num_components <- length(means)
  if (num_components == 0) {
    return(ggplot() + labs(title = "No Components to Plot"))
  }

  # --- 1. Create Grid ---
  x_coords <- seq(x_range[1], x_range[2], length.out = grid_res)
  y_coords <- seq(y_range[1], y_range[2], length.out = grid_res)
  grid_data <- expand_grid(x = x_coords, y = y_coords)

  # --- 2. Calculate Density for Each Gaussian ---
  density_list <- vector("list", num_components)
  for (i in 1:num_components) {
    density_list[[i]] <- dmvnorm(as.matrix(grid_data), mean = means[[i]], sigma = covariances[[i]])
  }

  # Combine densities into a matrix where columns are components
  density_matrix <- do.call(cbind, density_list)
  colnames(density_matrix) <- paste0("density_comp", 1:num_components)

  # --- 3. Determine Dominant Gaussian ---
  grid_data$dominant_comp <- apply(density_matrix, 1, which.max)
  grid_data$dominant_comp <- factor(grid_data$dominant_comp) # Convert to factor for coloring

  # --- 4. Prepare data for contours (if requested) ---
  contour_data <- NULL
  if (!is.null(contour_levels)) {
    contour_data_list <- vector("list", num_components)
    for (i in 1:num_components) {
      z_matrix <- matrix(density_list[[i]], nrow = grid_res, ncol = grid_res, byrow = FALSE)
      contour_data_list[[i]] <- data.frame(
        x = rep(x_coords, each = grid_res),
        y = rep(y_coords, times = grid_res),
        z = as.vector(z_matrix),
        component = factor(i)
      )
    }
    contour_data <- do.call(rbind, contour_data_list)
  }

  # --- 5. Plot using ggplot2 ---
  p <- ggplot(grid_data, aes(x = x, y = y)) +
    geom_raster(aes(fill = dominant_comp), interpolate = TRUE)

  # <--- MODIFIED SECTION: Conditionally apply color scale
  if (!is.null(component_colors)) {
    p <- p + scale_fill_manual(values = component_colors, name = "Dominant Component")
  } else {
    p <- p + scale_fill_viridis_d(name = "Dominant Component")
  }
  # --->

  p <- p +
    labs(title = title, x = x_axis_label, y = y_axis_label) +
    coord_fixed(ratio = 1) +
    theme_minimal() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.position = "none"
    )

  # Add contours if requested
  if (!is.null(contour_levels) && !is.null(contour_data)) {
    p <- p + geom_contour(data = contour_data, aes(z = z, color = component),
                          breaks = contour_levels, linewidth = 0.5, alpha = 0.7) +
      scale_color_brewer(palette = "Set1", name = "Component Contours", guide = "none")
  }

  # Add original data points if requested
  if (show_points && !is.null(data_points)) {
    p <- p + geom_point(data = data_points, aes(x = x, y = y),
                        color = point_color, size = 1, alpha = 0.6)
  }

  return(p)
}


plot_multiple_gaussian_density_zones <- function(list_of_means, list_of_covariances,
                                                 list_of_titles,
                                                 x_range = c(-5, 10), y_range = c(-5, 10),
                                                 grid_res = 100, contour_levels = NULL,
                                                 list_of_show_points = FALSE,
                                                 list_of_data_points = NULL,
                                                 component_colors = NULL,
                                                 list_of_point_colors = "black",
                                                 list_of_x_axis_labels = "X-axis", # Added new parameter
                                                 list_of_y_axis_labels = "Y-axis") { # Added new parameter
  


  # --- Input Validation for multiple plots ---
  num_plots <- length(list_of_means)
  if (num_plots == 0) {
    stop("Error: 'list_of_means' cannot be empty.")
  }
  if (length(list_of_covariances) != num_plots) {
    stop("Error: 'list_of_covariances' must have the same length as 'list_of_means'.")
  }
  if (length(list_of_titles) != num_plots) {
    stop("Error: 'list_of_titles' must have the same length as 'list_of_means'.")
  }

  # Handle show_points, data_points, point_colors potentially being single values
  if (length(list_of_show_points) == 1) {
    list_of_show_points <- rep(list_of_show_points, num_plots)
  } else if (length(list_of_show_points) != num_plots) {
    stop("Error: 'list_of_show_points' must be a single logical value or a vector of length equal to the number of plots.")
  }

  if (list_of_show_points[1] && is.null(list_of_data_points)) {
    # If the first show_points is TRUE and data_points is NULL, it means no data points provided
    # for any plot. This is an error if show_points is TRUE for any plot.
    if (any(list_of_show_points)) {
      stop("Error: If 'list_of_show_points' is TRUE for any plot, 'list_of_data_points' cannot be NULL.")
    }
  } else if (!is.null(list_of_data_points) && length(list_of_data_points) != num_plots) {
    stop("Error: 'list_of_data_points' must be NULL or a list of length equal to the number of plots.")
  }

  if (length(list_of_point_colors) == 1) {
    list_of_point_colors <- rep(list_of_point_colors, num_plots)
  } else if (length(list_of_point_colors) != num_plots) {
    stop("Error: 'list_of_point_colors' must be a single character string or a vector of length equal to the number of plots.")
  }

  # Handle axis labels potentially being single values
  if (length(list_of_x_axis_labels) == 1) {
    list_of_x_axis_labels <- rep(list_of_x_axis_labels, num_plots)
  } else if (length(list_of_x_axis_labels) != num_plots) {
    stop("Error: 'list_of_x_axis_labels' must be a single character string or a vector of length equal to the number of plots.")
  }

  if (length(list_of_y_axis_labels) == 1) {
    list_of_y_axis_labels <- rep(list_of_y_axis_labels, num_plots)
  } else if (length(list_of_y_axis_labels) != num_plots) {
    stop("Error: 'list_of_y_axis_labels' must be a single character string or a vector of length equal to the number of plots.")
  }


  all_plots <- list()

  for (i in 1:num_plots) {
    current_means <- list_of_means[[i]]
    current_covariances <- list_of_covariances[[i]]
    current_title <- list_of_titles[i]
    current_show_points <- list_of_show_points[i]
    current_data_points <- if (current_show_points) list_of_data_points[[i]] else NULL
    current_point_color <- list_of_point_colors[i]
    current_x_axis_label <- list_of_x_axis_labels[i] # Get current x-axis label
    current_y_axis_label <- list_of_y_axis_labels[i] # Get current y-axis label

    # Validate inner lists for means and covariances
    if (!is.list(current_means) || !all(sapply(current_means, function(m) is.numeric(m) && length(m) == 2))) {
      stop(paste0("Error: Element ", i, " of 'list_of_means' is not a list of 2-element numeric vectors."))
    }
    if (!is.list(current_covariances) || !all(sapply(current_covariances, function(s) is.matrix(s) && all(dim(s) == 2)))) {
      stop(paste0("Error: Element ", i, " of 'list_of_covariances' is not a list of 2x2 numeric matrices."))
    }
    if (length(current_means) != length(current_covariances)) {
      stop(paste0("Error: Inner lists of means and covariances for plot ", i, " must have the same length."))
    }

    # Validate covariance matrices for the current plot
    for (j in 1:length(current_means)) {
      if (!isSymmetric(current_covariances[[j]])) {
        stop(paste0("Error: Covariance matrix for component ", j, " in plot ", i, " is not symmetric."))
      }
      eigen_values <- eigen(current_covariances[[j]])$values
      if (any(eigen_values <= 0)) {
        stop(paste0("Error: Covariance matrix for component ", j, " in plot ", i, " is not positive definite."))
      }
    }

    # Create individual plot
    p <- plot_single_gaussian_density_zones(
      means = current_means,
      covariances = current_covariances,
      x_range = x_range,
      y_range = y_range,
      grid_res = grid_res,
      contour_levels = contour_levels,
      title = current_title,
      show_points = current_show_points,
      data_points = current_data_points,
      point_color = current_point_color,
      component_colors = component_colors,
      x_axis_label = current_x_axis_label, # Pass x-axis label
      y_axis_label = current_y_axis_label  # Pass y-axis label
    )
    all_plots[[i]] <- p
  }

  # Combine plots using patchwork
  # You can specify layout (e.g., `ncol = 2`) for more control
  combined_plot <- wrap_plots(all_plots)

  return(combined_plot)
}
