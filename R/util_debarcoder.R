########################
#debarcode_1 function, for debarcoding first level using mixture modeling.
########################
#fcb_df: barcodeded, compensated, data.frame

#bc_single_level: data.frame, sample singly stained with one level of channel to be debarcoded, usually compensation control, compensated, gated
#if no bc_single_level is provided, fcb_df will be substitued in (within regression script)

#channel: string, the channel to be debarcoded as a (eg: "Pacific Orange-A")

#levels: numeric, number of levels to find (eg; 6)

#uccutoff: numeric between 0 and 1, threshold for discarding uncertain cells, lower is more selective.

#opt: different mixture modeling algorithms
#opt1: equal variance, gaussian
#opt2: unequal variance, fits models with from levels:levels+2 and then combines clusters to result in the correct number of levels, good for contaminated normal
#opt3: depreciated
#opt4: fitting skew.normal distrubtion, seems to work best at the moment
#opt5: fitting skewt distrubtion, may be better than skew.normal but more parameters to fit so it takes longer
###############################################################################

#bc_single_level can't be null, change default?
debarcode_1 <- function(fcb_df, bc_single_level = NULL, channel, levels,
                        uccutoff = 0.05, opt = 4, subsample = 10e3, trans = 'arcsinh', 
                        updateProgress = NULL, cofactor_bc1 = NULL) {

  if (is.function(updateProgress)) {
    updateProgress(detail = "Mapping cellular Density")
  }

  area_density <- selectDenseScatterArea(bc_single_level, subsample = subsample)

  if (is.function(updateProgress)) {
    updateProgress(detail = "Performing Morphology Correction")
  }
  
  regression.output <- doRegressContrained(bc_single_level, fcb_df, Loc = area_density$loc, weight = area_density$c,
                                  trans = trans, columns = c(channel), monodir = c(1,1), cofactor = cofactor_bc1)

  cor.data <- regression.output[[1]]

  fcb_df2 <- fcb_df
  fcb_df2[, channel]<- cor.data[, channel]

  #prob better way to handle passing plot
  #mix_hist_plot <- hist(1)
  print(trans)
  if (trans == "log10") {
    print(trans)
    print(cofactor_bc1)
    print("opt was 4")
    vec <- log10(fcb_df2[,channel])
  } else if (trans == 'arcsinh') {
    print("opt was asinh")
    vec <- asinh(fcb_df2[,channel]/cofactor_bc1)
  }
  
  {
    vecss <- sample(vec, subsample)
    if (is.function(updateProgress)) {
      updateProgress(detail = "Initializing Model")
    }
    summary(vecss)
    mod.int <- classInt::classIntervals(vecss, levels, style = "fisher")
    classif <- sapply(vecss, function(x) pracma::findintervals(x, mod.int$brks))
    classif <- levels + 1 - classif
    mu.i <- as.numeric(unlist(lapply(split(vecss, classif), median))[-1])
    if (is.function(updateProgress)) {
      updateProgress(detail = "Optimizing Mixture Model")
    }


    Snorm.analysis <- mixsmsn::smsn.mix(vecss, nu = 3, g = levels, criteria = TRUE,
                               get.init = TRUE, group = TRUE, family = "Skew.normal", calc.im = FALSE, obs.prob = TRUE,
                               kmeans.param = list(iter.max = 20, n.start = 10, algorithm = "Hartigan-Wong"))


    loc <- Snorm.analysis$mu
    scale <- sqrt(Snorm.analysis$sigma2)
    shape <- Snorm.analysis$shape

    Snorm.df <- data.frame(loc, scale, shape)
    probs <- data.frame(x = vec)
    for (i in (1:nrow(Snorm.df))){
      probs[,as.character(i)] <- sn::dsn(vec, dp = as.numeric(Snorm.df[i,]))
    }
    probs.scaled <- apply(probs[,-1], 1, function(vec) { vec * Snorm.analysis$pii})

    probs.scaled.df <- as.data.frame(cbind(x = vec, t(probs.scaled)))


    probs.scaled.norm <- t(apply(t(probs.scaled), 1, function(vec) {vec/sum(vec)}))
    probs.scaled.norm.df <- as.data.frame(cbind(x = vec, t(probs.scaled)))

    #####
    #For plotting
    ####
    if (is.function(updateProgress)) {
      updateProgress(detail = "Generating Plots")
    }

    mysec <- seq(from = 0, to = 5, by = 0.01)
    model <- data.frame(x = mysec)
    for (i in (1:nrow(Snorm.df))){
      model[,as.character(i)] <- sn::dsn(mysec, dp = as.numeric(Snorm.df[i,]))

    }
    
    model.scaled <- apply(model[,-1], 1, function(vec) { vec * Snorm.analysis$pii})
    model.scaled<- cbind(x = model[,1], as.data.frame(t(as.matrix(model.scaled))))
    num.colnames <- as.numeric(colnames(model.scaled[-1]))
    colnames.order <- rev(order(Snorm.df$loc))

    colnames(model.scaled)[-1] <- match(num.colnames, colnames.order)

    melt.model.scaled <- reshape2::melt(model.scaled, id.vars = "x")
  

    vec.df <- data.frame(value = vec)
    if(trans == "log10"){
      vec.bt.df <- data.frame(value = 10^vec) #bt = back transform
      melt.model.scaled$x10 <- 10^melt.model.scaled$x
      mytrans <- "log10"
    } else if (trans =="arcsinh"){
      vec.bt.df <- data.frame(value = sinh(vec)*cofactor_bc1)
      melt.model.scaled$x10 <- sinh(melt.model.scaled$x)*cofactor_bc1
      mytrans <- asinh_trans(cofactor_bc1)
    }
    melt.model.scaled$variable<- factor(melt.model.scaled$variable,
                                        as.character(1:length(melt.model.scaled$variable)))
    
    print(colnames(melt.model.scaled))
    
    print(colnames(vec.bt.df))
    
    print(summary(melt.model.scaled))
    print(summary(vec.bt.df))
    
    #(melt.model.scaled$variable)
    
  
    mix.model.plot <- ggplot2::ggplot(melt.model.scaled,
                                      ggplot2::aes_string(x = "x10",
                                                          y = "value",
                                                          fill = "variable")) +
      ggplot2::geom_histogram(data = vec.bt.df, ggplot2::aes(x = value, y = ..density..),
                     inherit.aes = F, bins = 100, col = "black", fill = NA) +
      ggplot2::geom_area(alpha = 0.5, col = "grey50") +
      ggplot2::scale_x_continuous(trans = mytrans,
                                  breaks = major.ticks,
                    labels = major.ticks,
                    minor_breaks = minor.ticks,
                    name = channel) +
      ggplot2::coord_cartesian(xlim= quantile(vec.bt.df$value, c(0.0005,0.9995))) +
      ggplot2::ylab("Density") +
      ggplot2::scale_fill_discrete(name = "Population") +
      ggplot2::theme_classic()
    
    classif<- apply((probs.scaled.norm), 1, which.max)
    classif.uc <- apply(probs.scaled.norm, 1, max)
    classif <- match(classif, rev(order(Snorm.analysis$mu)))
    classif[classif.uc < (1- uccutoff)] <- 0 #

    fcb_df$bc1 <-  classif

    #don't need to return classif
    debarcoded_data <- list('df' = fcb_df,
                            "df2" = fcb_df2,
                            "channel" = channel,
                            'snorm' = Snorm.analysis,
                            "plot" = mix.model.plot,
                            'regressionmodel' = regression.output[["coefs"]])

    return(debarcoded_data)


  }

  print("Mixture modeling complete")

  fcb_df$bc1 <-  classif

  return(fcb_df)
}


debarcode.2 <- function(fcb_df, prevchannel, channel, levels, uccutoff = 0.05,
                        subsample = 10e3, trans = "arcsinh",
                        cofactor_bc1 = NULL,
                        cofactor_bc2 = NULL) {
  print('running db2')
  #previous channel should always be bc1
  prevlevel <- "bc1"

  fcb_df$bc2 <- 0
  bc1table <- table(fcb_df$bc1)

  bc1table.clean <- as.numeric(names(bc1table)[!(names(bc1table) == 0)])

  #if level = 1, exit function assign all bc2 to 1 (for single level barcodes)
  if (levels == 1) {
    fcb_df["bc2"] <- 1
    return(fcb_df)
  }

  for ( i in bc1table.clean) {

    #print(paste("Debarcoding Level:", i))
    ind <- fcb_df[,prevlevel] == i

    fcb_df.i <- fcb_df[which(ind), ]
    if (trans == "log10"){
      Y <- log10(fcb_df.i[,channel])
      X <- log10(fcb_df.i[,prevchannel])
    } else if (trans == "arcsinh"){
      Y <- asinh(fcb_df.i[,channel]/cofactor_bc2)
      X <- asinh(fcb_df.i[,prevchannel]/cofactor_bc1)
    }
    
    resids <- Y - X + median(Y, na.rm = TRUE)


    if(length(resids) > subsample) {
      resids.ss <- sample(resids, subsample)
    } else {
      resids.ss <- resids
    }

    mod.int <- classInt::classIntervals(resids.ss, levels, style = "fisher")

    classif <- lapply(resids, FUN = function(x) {findInterval(x, mod.int$brks)})

    classif <- unlist(classif)

    classif <- levels + 1 - classif
    classif[classif > levels] <- 0          
    fcb_df[which(ind), "bc2"] <- classif


  }
  return(fcb_df)
}
