list.files()
fcb_df <- read.csv("fcb_df.csv")
fcb_df
#colnames(fcb.df)
fcb_df$bc1 <- 0

uptake_channel <- "APC.H7.A"
channel <- "Pacific.Orange.A"
cofactor_bc1 <- 150
cofactor_uptake <- 150
fcb_df.i <- fcb_df
trans <- "arcsinh"
  if (trans == "log10"){
    Y <- log10(fcb_df.i[,channel])
    X <- log10(fcb_df.i[,uptake_chanel])
  } else if (trans == "arcsinh"){
    Y <- asinh(fcb_df.i[,channel]/cofactor_bc1)
    X <- asinh(fcb_df.i[,uptake_channel]/cofactor_uptake)
  }
  
  resids <- Y - X + median(Y, na.rm = TRUE)
  
  hist(resids, n = 300)
  vecss <- resids[sample(1:length(resids), 10e3, replace = TRUE)]
  levels <-6
  Snorm.analysis <- mixsmsn::smsn.mix(vecss, nu = 3, g = levels, criteria = TRUE,
                                      get.init = TRUE, group = TRUE, family = "Skew.normal", calc.im = FALSE, obs.prob = TRUE,
                                      kmeans.param = list(iter.max = 20, n.start = 10, algorithm = "Hartigan-Wong"))
  
  #library(mixsmsn)
  
  loc <- Snorm.analysis$mu
  scale <- sqrt(Snorm.analysis$sigma2)
  shape <- Snorm.analysis$shape
  Snorm.df <- data.frame(loc, scale, shape)
  probs <- data.frame(x = resids)
  for (i in (1:nrow(Snorm.df))){
    probs[,as.character(i)] <- sn::dsn(resids, dp = as.numeric(Snorm.df[i,]))
  }
  probs.scaled <- apply(probs[,-1], 1, function(vec) { vec * Snorm.analysis$pii})
  summary(probs.scaled)
  probs.scaled.df <- as.data.frame(cbind(x = resids, t(probs.scaled)))
  
  
  probs.scaled.norm <- t(apply(t(probs.scaled), 1, function(vec) {vec/sum(vec)}))
  probs.scaled.norm.df <- as.data.frame(cbind(x = resids, t(probs.scaled)))
  summary(probs.scaled.norm.df)
  
  
  
  
  vec<-resids
  
  mysec <- seq(from = min(round(vec,2)), to = max(round(vec,2)), by = 0.01)
  model <- data.frame(x = mysec)
  for (i in (1:nrow(Snorm.df))){
    model[,as.character(i)] <- sn::dsn(mysec, dp = as.numeric(Snorm.df[i,]))
    
  }
  
  # write.csv(model, "model.csv")
  # write.csv(Snorm.analysis$pii, "snormpii.csv")
  # write.csv(vec, "vec.csv")
  model.scaled <- apply(model[,-1], 1, function(vec) { vec * Snorm.analysis$pii})
  #model.scaled <- model.scaled/max(model.scaled)
  model.scaled<- cbind(x = model[,1], as.data.frame(t(as.matrix(model.scaled))))
    str(resids)
  model.scaled.melted <- reshape2::melt(model.scaled, id.vars = "x")
  ggplot(model.scaled.melted, aes( x= x, y = value, fill = variable)) +
    geom_area() + 
    geom_histogram(data = data.frame(x = resids), aes(x = x, y = ..density..),
                   inherit.aes = F, fill = "NA", col = "black", bins = 100) + 
    theme_classic()
  
  
  
  
  ###
  debarcode_1 <- c(fcb_df, bc_single_level = NULL, channel, levels,
                          uccutoff = 0.05, opt = 4, subsample = 10e3, trans = 'arcsinh', 
                          updateProgress = NULL,
                          cofactor_bc1 = NULL, cofactor_uptake = NULL,
                          uptakecx = NULL,
                          likelihoodcut = 1)
  dbc1 <- debarcode_1(fcb_df, channel = "Pacific.Orange.A", levels = 6,
              cofactor_bc1 = 150, cofactor_uptake = 150,
              uptakecx = "APC.H7.A")
  colnames(fcb_df)
  dbc1 <- debarcode_1(fcb_df, fcb_df, channel = "Pacific-Orange-A", levels = 6,
                      cofactor_bc1 = 150, subsample = 10e3)
  dbc1$plot
  
  
  
##########3
  
 
      subsample <- 10e3
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
      resids.split <- split(resids, classif)
      
      hist.probs <- list()
      for (i in as.character(1:max(names(resids.split)))){
            myhist <- hist(resids.split[[i]],100, plot = FALSE)
            binprobs <- myhist$counts/sum(myhist$counts)
            hist.probs.i<- rep(0, times = length(resids))
            bin.assingments <- findInterval(resids, myhist$breaks)
            hist.probs.i[which(bin.assingments != 0)] <- binprobs[bin.assingments]
            hist.probs.i[which(is.na(hist.probs.i))] <- 0
            hist.probs[[i]] <- hist.probs.i
      }
      hist.probs.m <- do.call(cbind, hist.probs)
      hist.probs.m
      
      likelihoodcut <- 8
      ind <- myhist$counts > max(myhist$counts)/likelihoodcut
      
      
      resids
      valid.ints <- cbind(myhist$breaks[which(ind)], myhist$breaks[which(ind) + 1])
      ind
      ints <- findInterval(resids, myhist$breaks)
      resids[ints %in% (which(ind))]
      myhist2 <- hist(resids[ints %in% (which(ind))], n = 100, xlim = c(5.5, 5.95))
      myhist2
      max(myhist2$counts) / min(myhist2$counts)
      ?hist  
      
      #valid.ints <- numeric()
      valid.ints
      myhist$counts
      ?findInterval
      
      hist(resids.ss, n = 300)
      
      classif
      
      fcb_df[which(ind), "bc2"] <- classif

  