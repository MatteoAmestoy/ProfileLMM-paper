# Plot results from real data

setwd('/Lifeline_server/')
storeMain = readRDS('profile_DiasStab.RData')
storeStab = readRDS('profile_DiasMain.RData')
# out = readRDS('DiasMain.RData')

source("paper_utility.R")
library(openxlsx)
library(ggplot2)



varNames = read.xlsx('VarNames.xlsx')
expLabel = c('main','stab')



idxClusMain = which(table(storeMain$optClus)>100)
idxClusStab = which(table(storeStab$optClus)>100)
cenMain = storeMain$cen[,idxClusMain]
cenStab = storeStab$cen[,idxClusStab]
simMat = t(cenMain)%*%cenStab

storeMain$order = c(5,4,2,1,3)
storeStab$order = unname(idxClusStab[c(2,4,1,3)])
table(storeMain$optClus)[storeMain$order]

table(storeStab$optClus)[storeStab$order]

for(store in list(storeMain,storeStab)){
  store$gamma = store$gamma[,store$order,]
  store$centroids = store$centroids[,store$order,]
  store$coVar = store$coVar[,,store$order]
  store$cen = store$cen[,store$order]
  store$gam = store$gam[,store$order]

  nC = length(store$order)
  cen = store$cen
  coVar =  store$coVar
  gam =  store$gam

  pct = 5
  apply(store$BetaFE, 1, median, na.rm=TRUE)[
    apply(store$BetaFE, 1, quantile, (pct/2)/100)*apply(store$BetaFE, 1, quantile, (100-pct/2)/100)>0]


  # Fixed effect estimates ------------------------------------------


  coef_data <- data.frame(
    term = varNames$X2,
    estimate = unname(apply(store$BetaFE, 1, median, na.rm=TRUE)),
    conf.low = unname(apply(store$BetaFE, 1, quantile, (pct/2)/100)),
    conf.high = unname(apply(store$BetaFE, 1, quantile, 1-(pct/2)/100))
  )

  coef_data$term <- factor(coef_data$term, levels = coef_data$term)

  windows()
  ggplot(coef_data[2:26,], aes(x = estimate, y = term,
                               xmin = conf.low, xmax = conf.high)) +
    geom_pointrange(
      size = 1,       # Thickness of the line
      linewidth = 2,
      color = "steelblue" # Color of the points and lines
    ) +
    geom_vline(
      xintercept = 0,   # Add a vertical line at 0 (for no effect)
      linetype = "dashed",
      color = "gray50",
      size = 0.9
    ) +
    labs(
      title = "Fixed effect Coefficient Estimates with 95% Credible Intervals",
      x = "Estimate",
      y = "Coefficient"
    ) +
    theme_minimal() + # A clean theme
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"), # Center and bold title
      axis.text.y = element_text(size = 10), # Adjust y-axis label size
      axis.title.x = element_text(margin = margin(t = 10)) # Add margin to x-axis title
    )




  # CLuster representations -------------------------------------------
  library(RColorBrewer)
  colors_for_labels <- brewer.pal(n = nC, name = "Set3") # if more than 9 clusters
  colors_for_labels <- brewer.pal(n = nC, name = "Set1")




  lim =1.4*c(min(store$centroids),max(store$centroids))
  lay = c()
  cum = 0
  xxx={}
  for(j in 1:4){
    lay = c(lay,rep(0,j-1),(1:(4-j+1))+cum)
    cum =cum+(4-j+1)
    xxx[[j]]=(-100:100)/100*lim[2]*2
  }
  windows()
  layout.matrix <- matrix(lay, nrow = 4, ncol = 4)
  layout.matrix = layout.matrix+(layout.matrix>0)
  layout.matrix[1,4] = 1
  layout(mat = layout.matrix,
         heights = rep(1,4), # Heights of the two rows
         widths = rep(1,4))
  par(mar = c(0, 0, 0, 0))
  plot(
    1,
    1,
    type = "n", # 'n' prevents any plotting, creating a blank canvas
    ann = FALSE,
    axes = FALSE
  )

  # --- 3. Add the legend to the empty plot ---
  # The 'legend()' function is used to place the legend.
  legend(
    "center",                  # Position the legend in the center of the plot
    legend = 1:nC,    # The text for each label
    col = colors_for_labels,       # The colors associated with each label
    pch = 15,                  # The point character to use (15 is a filled square)
    pt.cex = 3,                # Size of the legend points
    bty = "n",                 # 'n' means no border box around the legend
    title = "Cluster number",    # Title for the legend
    cex = 2                 # Size of the legend text
  )


  for (i in 1:4){
    for (j in i:4){
      par(mar = rep(2, 4))
      first = T
      if(i==j){

        for (c in 1:nC){
          if(first){
            plot(xxx[[i]],dnorm(xxx[[i]],cen[i,c],sqrt(coVar[i,i,c])),
                 col=colors_for_labels[c],
                 main =rownames(cen)[i],
                 xlim=2*lim,
                 type = 'l')
            first = F
          }else{lines(xxx[[i]],dnorm(xxx[[i]],cen[i,c],sqrt(coVar[i,i,c])),col=colors_for_labels[c])}
        }


      }else{
        plot(cen[i,],t(cen[j,]),col=colors_for_labels,xlim=lim,ylim=lim)
        for (c in 1:nC){
          aaa = ellipse(
            cen[ c(i,j), c],
            coVar[ c(i,j),c(i,j) ,c]
            ,0.3,100)
          lines(aaa$xx,aaa$yy,
                col = colors_for_labels[c],
                lwd = 3)
        }}

    }}




  # Cluster effect estimation-----------------------------------------------

  # windows()
  # lim = {}
  # for (i in 1:4){
  #   lim[[i]] = c(min(store$gamma[i,,]),max(store$gamma[i,,]))
  # }
  #
  #
  # par(mfrow = c(1,4), col.sub="blue", cex.sub=2)
  # for (i in 1:4){
  #   first = T
  #   for (c in 1:nC){
  #     density_obj <- density(store$gamma[i,c,])
  #     if(first){
  #       plot(density_obj,
  #            xlim=lim[[i]],
  #            col=colors_for_labels[c],
  #            main =rownames(store$gam)[i] ,
  #            type = 'l')
  #       first = F
  #     }else{lines(density_obj,col=colors_for_labels[c])}
  #   }
  # }
  #

  bbb = array(0,dim = dim(store$gamma)-c(0,1,0))
  for(c in 2:nC){
    bbb[,c-1,] =store$gamma[,c,]-store$gamma[,1,]
  }
  custom_colors <- c(
    "1" = colors_for_labels[1],
    "2" = colors_for_labels[2],
    "3" = colors_for_labels[3],
    "4" = colors_for_labels[4],
    "5" = colors_for_labels[5]
  )
  pct = 10
  plo={}
  for (i in 1:4){
    coef_data <- data.frame(
      term = as.character( 2:nC),
      estimate = unname(apply(bbb[i,,], 1, median, na.rm=TRUE)),
      conf.low = unname(apply(bbb[i,,], 1, quantile, (pct/2)/100)),
      conf.high = unname(apply(bbb[i,,], 1, quantile, 1-(pct/2)/100))
    )

    coef_data$term <- factor(coef_data$term, levels = coef_data$term)


    plo[[i]]=ggplot(coef_data, aes(x = estimate, y = term,
                                   xmin = conf.low, xmax = conf.high,
                                   color = term)) +
      geom_pointrange(
        size = 1.5,       # Thickness of the line
        linewidth = 1.5      # Size of the point
      ) +
      geom_vline(
        xintercept = 0,   # Add a vertical line at 0 (for no effect)
        linetype = "dashed",
        color = "gray50",
        size = 0.5
      ) +
      labs(
        title = paste0("Coefficient Estimates with ",100-pct,"% Credible Intervals for covariate ", rownames(store$gam)[i]),
        x = "Estimate",
        y = "Cluster"
      ) +
      theme_minimal() + # A clean theme
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold"), # Center and bold title
        axis.text.y = element_text(size = 10), # Adjust y-axis label size
        axis.title.x = element_text(margin = margin(t = 10)), # Add margin to x-axis title
        legend.position = "none" # Hide the legend as terms are on y
      )+
      scale_color_manual(values = custom_colors)
  }
  windows()
  plo[[1]]+plo[[2]]+plo[[3]]+plo[[4]]+
    plot_annotation(title = "(Cluster 1 as reference)")&
    theme(plot.title = element_text(hjust = 0.5))
}
