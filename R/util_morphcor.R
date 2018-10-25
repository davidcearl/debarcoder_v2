

morphology.corr <- function(fcb_df, bc_single_level = NULL, channel,
                       opt = "regression", subsample = 10e3, trans = 'arcsinh', 
                       updateProgress = NULL,
                       cofactor_bc1 = NULL, cofactor_uptake = NULL,
                       uptake_channel = NULL) {
  # fcb_df
  # bc_single_level <- fcb_df
  # channel <- "Pacific-Orange-A"
  # levels <- 6
  # cofactor_bc1 <- 150
  # opt <- "regression"
  # updateProgress <- NULL
  # subsample <- 10e3
  # trans  <- "arcsinh"
  if(opt == "regression") {
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
    ####################################
  } else if(opt == "controldye") {
    if(is.null(uptake_channel)){
      warning("No uptake control channel specified")
      return()
    } 
    if(is.null(cofactor_uptake) & trans == "arcsinh"){
      warning("No uptake control cofactor specified")
      return()
    }
    if (trans == "log10"){
      Y <- log10(fcb_df[,channel])
      X <- log10(fcb_df[,uptake_channel])
    } else if (trans == "arcsinh"){
      Y <- asinh(fcb_df[,channel]/cofactor_bc1)
      X <- asinh(fcb_df[,uptake_channel]/cofactor_uptake)
    }
    
    vec <- Y - X + median(Y, na.rm = TRUE)
    
    if (trans == "log10"){
      vec <- vec^10
    } else if (trans == "arcsinh"){
      vec <- sinh(vec)*cofactor_bc1
    }
    fcb_df2 <- fcb_df
    fcb_df2[, channel]<- vec
    

  }
  return(fcb_df2[, channel])
}

#neeed to make sure levels are in order!
fit.models <- function(vec, #vector of barcoding intensities, output of morphology.corr
                       levels, #number of levels
                       opt = "mixture", #mixture (guassian mixture models) or fisher (univariate k-means)
                       dist = NULL, #for gaussian mixture models, Skew.normal, normal, T.dist
                       trans = "arcsinh", #asinh or log10
                       cofactor = NULL, 
                       subsample = 10e3, 
                       updateProgress = NULL){#cofactor for asinh transofrmation

  if (trans == "log10"){
    vec <- log10(vec)
  } else if (trans == "arcsinh"){
    vec <- asinh(vec/cofactor)
  }
  
  vecss <- sample(vec, subsample, replace = TRUE)
  
  if(opt == "mixture") {
    #cofactor <- 150

    
    if (levels > 1) {
      mod.int <- classInt::classIntervals(vecss, levels, style = "fisher")
      classif <- sapply(vecss, function(x) pracma::findintervals(x, mod.int$brks))
      classif <- levels + 1 - classif
      mu.i <- as.numeric(unlist(lapply(split(vecss, classif), median))[-1])
    }
    if (is.function(updateProgress)) {
      updateProgress(detail = "Optimizing Mixture Model")
    }
    # dist <- "Skew.normal"
    # levels <- 4
    # ?mixsmsn::smsn.mix
    # ptm <- proc.time()
     Snorm.analysis <- mixsmsn::smsn.mix(vecss, nu = 3, g = levels, criteria = TRUE,
                                        get.init = TRUE, group = TRUE, family = dist, calc.im = FALSE, obs.prob = TRUE,
                                        kmeans.param = list(iter.max = 20, n.start = 10, algorithm = "Hartigan-Wong"))
    # proc.time() - ptm
    # mix.hist(vecss, Snorm.analysis, breaks = 50)
    # (Snorm.analysis)
    loc <- Snorm.analysis$mu
    scale <- sqrt(Snorm.analysis$sigma2)
    shape <- Snorm.analysis$shape
    
    Snorm.df <- data.frame(loc, scale, shape)
    probs <- data.frame(x = vec)
    
    for (i in (1:nrow(Snorm.df))){
      probs[,as.character(i)] <- sn::dsn(vec, dp = as.numeric(Snorm.df[i,]))
    }
    probs.scaled <- apply(probs[,-1], 1, function(vec) { vec * Snorm.analysis$pii})
    
    probs.scaled.df <- as.data.frame(t(probs.scaled))

    
  } else if(opt == "fisher") {
    mod.int <- classInt::classIntervals(vecss, levels, style = "fisher")
    
    classif <- lapply(vec, FUN = function(x) {findInterval(x, mod.int$brks)})
    
    classif <- unlist(classif)
    
    classif <- levels + 1 - classif
    classif[classif > levels] <- 0          
    
    
    # ggplot(max(vec.split[[2]]))
    # preplot <-ggplot(data.frame(x = vec.split[[3]]), aes(x = x)) + 
    #   geom_histogram(bins = 100, col = "black", fill = "grey99") + 
    #   geom_vline(data = data.frame(x = mod.int$brks), aes(xintercept = x), 
    #              linetype = 2, size = 1)
    # preplot <- ggplot_build(preplot)
    # preplot$data[[1]][["count"]]
    # ggplot(data.frame(x = vec.split[[3]]), aes(x = x)) + 
    #   geom_histogram(bins = 100, col = "black", fill = "grey99") + 
    #   geom_vline(data = data.frame(x = mod.int$brks), aes(xintercept = x), 
    #              linetype = 2, size = 1) + 
    #   geom_hline(yintercept = max(preplot$data[[1]][["count"]])/8,
    #              linetype = 4, size = 1, col = "blue") + 
    #   theme_classic()

    #calculate the empircal probability for each cell belonging to each 
    #population based on the histogram
    
    vec.split <- split(vec, classif)
    
    hist.probs <- list()
    for (i in as.character(1:max(names(vec.split)))){
      myhist <- hist(vec.split[[i]],100, plot = FALSE)
      binprobs <- myhist$counts/sum(myhist$counts)
      hist.probs.i<- rep(0, times = length(vec))
      bin.assingments <- findInterval(vec, myhist$breaks)
      hist.probs.i[which(bin.assingments != 0)] <- binprobs[bin.assingments]
      hist.probs.i[which(is.na(hist.probs.i))] <- 0
      hist.probs[[i]] <- hist.probs.i
    }
    hist.probs.m <- do.call(cbind, hist.probs)
    probs.scaled.df <- as.data.frame(hist.probs.m)  
    
  }
  colMax <- apply(probs.scaled.df, 2, max)
  probs.rescale.col <- sweep(probs.scaled.df, 2, colMax, FUN="/")
  return(probs.rescale.col)
} 


assign.cells <- function(fcb_df, probs, likelihoodcut, ambiguitycut, output = "classif"){
  likelihoodcut <- 8
  ambiguitycut <- 0.05
  classif <- rep(0, nrow(fcb_df))
  probs.norm.row <- t(apply(probs, 1, function(vec) {vec/sum(vec)}))
  classif <- unlist(apply(probs.norm.row, 1, which.max))
  likely <- probs > 1/likelihoodcut
  likely.sum <- (apply(likely, 1, sum))
  classif[likely.sum != 1] <- 0
  non.ambigious <- apply(probs.norm.row, 1, max) > (1 - ambiguitycut)
  classif[!non.ambigious] <- 0
  if(output == "classif") {
    return(classif)
  } else if (output == "plot") {
    preplot<- ggplot(mydf, aes(x = vec, fill = classif)) + 
      geom_histogram(bins = 400) + 
      scale_x_continuous(trans = asinh_trans(150)) + 
      scale_y_continuous(expand = c(0,0)) +
      theme_classic()
    preplot.build<- ggplot_build(preplot)
    preplot.data <- preplot.build$data[[1]]
    preplot.data.list <- split(preplot.data, preplot.data$group)
    
    bounds <- list()
    for (i in 1:    ncol(probs)) {
      bounds[[i]]<- c(classif = i, min = min(vec[likely[,i]]), max = max(vec[likely[,i]]))
    }
    
    thresholds<- data.frame(do.call(rbind, bounds))
    thresholds$classif <- as.factor(thresholds$classif)
    thresholds$y<- do.call(rbind, lapply(preplot.data.list, function(mydf) max(mydf[,"count"])/likelihoodcut))[-1]
    
    ggplot(mydf, aes(x = vec, fill = classif)) + 
      geom_histogram(bins = 400) + 
      geom_segment(data = thresholds, aes(x = min, xend = max,
                                         y = y, yend = y),
                   col = "black", inherit.aes = FALSE) + 
      scale_x_continuous(trans = asinh_trans(150)) + 
      scale_y_continuous(expand = c(0,0)) +
      theme_classic()
    
    
    
    str(preplot.build$data)
  }
  
  # hist(apply(probs.norm.row,1, max), n = 21)
}

fcb_df <- read.csv("fcb_df.csv", stringsAsFactors = FALSE)
colnames(fcb_df) <- gsub("\\.","-", colnames(fcb_df))
colnames(fcb_df)
head(fcb_df[2:11])
nrow(fcb_df)
hist(log10(vec), n = 200)

vec1.d <- morphology.corr(fcb_df,
                       fcb_df,
                       channel = "Pacific-Orange-A",
                       cofactor_bc1 = 150, 
                       opt = "controldye",
                       uptake_channel = "APC-H7-A", 
                       cofactor_uptake = 150)
vec2.d <- morphology.corr(fcb_df,
                        fcb_df,
                        channel = "Pacific-Blue-A",
                        cofactor_bc1 = 150, 
                        opt = "controldye",
                        uptake_channel = "APC-H7-A", 
                        cofactor_uptake = 150)

vec1 <- morphology.corr(fcb_df,
                        fcb_df,
                        channel = "Pacific-Orange-A",
                        cofactor_bc1 = 150, 
                        opt = "regression")
vec2 <- morphology.corr(fcb_df,
                        fcb_df,
                        channel = "Pacific-Blue-A",
                        cofactor_bc1 = 150, 
                        opt = "regression")


mynewdf <- data.frame("PacBlue" = vec1, "PacOrange" = vec2)
mynewdf <- data.frame("PacBlue" = vec1,
                      "PacOrange" = vec2,
                      "Ax750" = fcb_df$`APC-H7-A`,
                      "SSC" = fcb_df$`SSC-A`, 
                      "PacOrangeOrig" = fcb_df$`Pacific-Orange-A`,
                      "PacBlueOrig" = fcb_df$`Pacific-Blue-A`)


ggplot(mynewdf, aes(x= PacBlue, y = PacOrange)) + 
  geom_bin2d(bins = 256) + 
  scale_x_continuous(trans = asinh_trans(150)) + 
  scale_y_continuous(trans = asinh_trans(150)) + 
  viridis::scale_fill_viridis(option = "A")
ggplot(mynewdf, aes(x= Ax750, y = PacOrange)) + 
  geom_bin2d(bins = 256) + 
  scale_x_continuous(trans = asinh_trans(150)) + 
  scale_y_continuous(trans = asinh_trans(150)) + 
  viridis::scale_fill_viridis(option = "A")

vec <- morphology.corr(fcb_df,
                       fcb_df,
                       channel = "Pacific-Orange-A",
                       cofactor_bc1 = 150, 
                       opt = "regression")


probs <- fit.models(vec,
                    levels = 6,
                    opt = "mixture",
                    dist = "Skew.normal",
                    cofactor = 150)

assignments <- assign.cells(fcb_df,
                            probs,
                            likelihoodcut = 8,
                            ambiguitycut = 0.02)
table(assignments)
hist(log10(vec), n = 200, xlim = c(1.8, 4.2))
hist(log10(vec[assignments == 6]), n = 20, xlim = c(1.8, 4.2))
