source("./debarcoder/regression.R")
#source("gates_transforms.R")
require("mclust")
require("MASS")
#install.packages("mixsmsn")
require("classInt")
require("mixsmsn")
require("ggplot2")
require("viridis")
require("reshape2")
require(sn)

#opt <- 1


if (1 == 2) { #for debugging
  mydf.pure <- po.df
  channel <- "Pacific Orange-A"
  levels <- 6
  uccutoff <- 0.05
}

########################
#debarcode.1 function, for debarcoding first level using mixture modeling. 
########################
#mydf: barcodeded, compensated, data.frame

#mydf.pure: data.frame, sample singly stained with one level of channel to be debarcoded, usually compensation control, compensated, gated
#if no mydf.pure is provided, mydf will be substitued in (within regression script)

#channel: string, the channel to be debarcoded as a (eg: "Pacific Orange-A")

#levels: numeric, number of levels to find (eg; 6)

#uccutoff: numeric between 0 and 1, threshold for discarding uncertain cells, lower is more selective. 

#opt: different mixture modeling algorithms
#opt1: equal variance, gaussian
#opt2: unequal variance, fits models with from levels:levels+2 and then combines clusters to result in the correct number of levels, good for contaminated normal
#opt3: depreciated
#opt4: fitting skew.normal distrubtion, seems to work best at the moment
#opt5: fitting skewt distrubtion, may be better than skew.normal but more parameters to fit so it takes longer

debarcode.1 <- function(mydf, mydf.pure = NULL, channel, levels, uccutoff = 0.05, opt = 4, subsample = 10e3,
                        updateProgress = NULL) {
  mydata <-mydf.pure
  mydata2 <- mydf

  if (is.function(updateProgress)) {  
    updateProgress(detail = "Mapping cellular Density")
  }
  
  k.da <- selectDenseScatterArea(mydata, subsample = subsample)
  
  if (is.function(updateProgress)) {  
    updateProgress(detail = "Performing Morphology Correction")
  }
  regression.output <- doRegressContrained(mydata, mydata2, k.da$loc, k.da$c,
                                  opt='logF', plot_flag = 0,
                                  columns = c(channel), MonoDir = c(1,1))
  
  cor.data <- regression.output[[1]]
  #channel
  #hist(log10(mydata2$`Pacific Orange-A`), n = 200)
  #mydata$`Pacific Orange-A`
  #hist(log10(cor.data[[channel]]), n = 200)
  # str(cor.data)
  # ?sample()
  # #hist(, n = 100)
  # mysample <- sample(1:nrow(cor.data), 5000, replace = TRUE)
  # 
  # skewmm <- smsn.mix(log10(cor.data$`Pacific Orange-A`[mysample]), nu = 1,
  #                    g = 6, family = "Normal", obs.prob = TRUE, error = 0.001)
  # 
  # ?smsn.mix
  # skewt <- fmmst(g = 6, dat = cbind(log10(cor.data$`Pacific Orange-A`[mysample],rnorm(5000))))
  # cbind(cor.data$`Pacific Orange-A`[mysample])
  # 
  # mu1 <- 5; mu2 <- 20; mu3 <- 35
  # sigma2.1 <- 9; sigma2.2 <- 16; sigma2.3 <- 9
  # lambda1 <- 5; lambda2 <- -3; lambda3 <- -6
  # nu = 5
  # 
  # mu <- c(mu1,mu2,mu3)
  # sigma2 <- c(sigma2.1,sigma2.2,sigma2.3)
  # shape <- c(lambda1,lambda2,lambda3)
  # pii <- c(0.5,0.2,0.3)
  # 
  # arg1 = c(mu1, sigma2.1, lambda1, nu)
  # arg2 = c(mu2, sigma2.2, lambda2, nu)
  # arg3 = c(mu3, sigma2.3, lambda3, nu)
  # y <- rmix(n=1000, p=pii, family="Skew.t", arg=list(arg1,arg2,arg3))
  # 
  # hist(y, n = 30)
  # Norm.analysis <- smsn.mix(log10(cor.data$`Pacific Orange-A`), nu = 3, g = 6,
  #                           get.init = TRUE, criteria = TRUE, 
  #                           group = TRUE, family = "Normal", calc.im=FALSE)
  # mix.hist(log10(cor.data$`Pacific Orange-A`),Norm.analysis)
  # mix.print(Norm.analysis)
  # mix.dens(log10(cor.data$`Pacific Orange-A`),Norm.analysis)
  # 
  # ptm <- proc.time()
  # Snorm.analysis <- smsn.mix(log10(cor.data$`Pacific Orange-A`), nu = 3, g = 6,
  #                           get.init = TRUE, criteria = TRUE, 
  #                           group = TRUE, family = "Skew.normal", calc.im = FALSE, obs.prob = TRUE)
  # proc.time() - ptm
  # 
  # ptm <- proc.time()
  # St.analysis <- smsn.mix(log10(cor.data$`Pacific Orange-A`), nu = 3, g = 6,
  #                            get.init = TRUE, criteria = TRUE, 
  #                            group = TRUE, family = "Skew.t", calc.im = FALSE, obs.prob = TRUE)
  # proc.time() - ptm
  # 
  # mix.hist(log10(cor.data$`Pacific Orange-A`),Snorm.analysis)
  # mix.hist(log10(cor.data$`Pacific Orange-A`),St.analysis)
  # 
  # mix.print(Norm.analysis)
  # mix.print(Snorm.analysis)
  # mix.print(St.analysis)
  # mix.dens(log10(cor.data$`Pacific Orange-A`),Norm.analysis)
  # library(ggplot2)
  # qplot(log10(cor.data$`Pacific Orange-A`),   1- apply(Snorm.analysis$obs.prob, 1, max))
  # str(log10(cor.data$`Pacific Orange-A`))
  # skewt
  # ?fmmst
  # skewmm
  # ?smsn.mix
  #str(mydf)
  #str(cor.data)
  mydf2 <- mydf
  mydf2[, channel]<- cor.data[, channel]
  if (opt == 1) { #option 1: Mixture modeling on log transformed data, 6 clusters,  equal variances
    mix.model <- Mclust(log(mydf2[,channel]), G = 6 ,model = "E")
    #plot(mix.model)
    classif <- mix.model$classification
    mydf.split2 <- split(mydf, classif)
    classif <- match(classif, rev(order(tail(unlist(lapply(mydf.split2, function(df) mean(df[,channel])))))))
    classif[mix.model$uncertainty > uccutoff] <- 0 #
  }

  if (opt ==2) { #option2: mixture modeling on log transformed data, 6:8 clusters, combinded afterwards
    #good option for peaks which contain multiple populations with same mean, different variances
    mix.model <- Mclust(log(mydf2[,channel]), G =   levels:(levels+2))
    combi <- clustCombi(log(mydf2[,channel]), mix.model)

    #plot(combi, log(mydf2[,channel]))
    
    classif <- combi$classification[[6]]
    classif.uc <- apply(combi$combiz[[6]],1, max)
    mydf.split2 <- split(mydf, classif)
    classif <- match(classif, rev(order(tail(unlist(lapply(mydf.split2, function(df) mean(df[,channel])))))))
    plot(classif.uc ~ log(mydf2[,channel]))
    classif[classif.uc < (1- uccutoff)] <- 0 #
  
  }
  
  if (opt ==3) { #option3: work in progress
    mix.model <- Mclust(log(mydf2[,channel]), G =   levels)
    plot(mix.model, what = "density")
    plot(mix.model, what = "BIC")
    plot(mix.model, what = "classification")
    
    combi <- clustCombi(log(mydf2[,channel]), mix.model)
    classif <- combi$classification[[6]]
    classif.uc <- apply(combi$combiz[[6]],1, max)
    mydf.split2 <- split(mydf, classif)
    classif <- match(classif, rev(order(tail(unlist(lapply(mydf.split2, function(df) mean(df[,channel])))))))
    classif[classif.uc < (1- uccutoff)] <- 0 #
    
  }
  
  #prob better way to handle passing plot
  #mix_hist_plot <- hist(1)
  if (opt == 4) {
    vec <- log10(mydf2[,channel])

    vecss <- sample(vec, subsample)
    if (is.function(updateProgress)) {
      updateProgress(detail = "Initializing Model")
    }
    
    mod.int <- classIntervals(vecss, levels, style = "fisher")
    classif <- sapply(vecss, function(x) findintervals(x, mod.int$brks))
    classif <- levels + 1 - classif
    mu.i <- as.numeric(unlist(lapply(split(vecss, classif), median))[-1])
    if (is.function(updateProgress)) {
      updateProgress(detail = "Optimizing Mixture Model")
    }
    
    
    Snorm.analysis <- smsn.mix(vecss, nu = 3, g = levels, criteria = TRUE,
                               get.init = TRUE, group = TRUE, family = "Skew.normal", calc.im = FALSE, obs.prob = TRUE, 
                               kmeans.param = list(iter.max = 20, n.start = 10, algorithm = "Hartigan-Wong"))
    
    #mix.hist(log10(mydf2[,channel]),Snorm.analysis, 120)
    
    loc <- Snorm.analysis$mu 
    scale <- sqrt(Snorm.analysis$sigma2)
    shape <- Snorm.analysis$shape

    Snorm.df <- data.frame(loc, scale, shape)
    probs <- data.frame(x = vec)
    for (i in (1:nrow(Snorm.df))){
      probs[,as.character(i)] <- dsn(vec, dp = as.numeric(Snorm.df[i,]))
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
      model[,as.character(i)] <- dsn(mysec, dp = as.numeric(Snorm.df[i,]))
      
    }
    model.scaled <- apply(model[,-1], 1, function(vec) { vec * Snorm.analysis$pii})
    model.scaled<- cbind(x = model[,1], as.data.frame(t(as.matrix(model.scaled))))
    num.colnames <- as.numeric(colnames(model.scaled[-1]))
    colnames.order <- rev(order(Snorm.df$loc))
    
    colnames(model.scaled)[-1] <- match(num.colnames, colnames.order)
    #model.scaled[2:7] <- model.scaled[as.character(1:6)]
    #str(model.scaled)
    # colnames.order
    # model.classif <- apply(model.scaled[2:7], 1, function(vec) {which.max(vec)})
    # model.classif <- match(model.classif, rev(order(Snorm.df$loc)))
    # 
    # model.sum.scaled <- apply(model.scaled[2:7], 1, sum)
    # model.scaled.df <- data.frame(x = model[,1], value = model.sum.scaled, variable = model.classif)
    # 
    
    melt.model.scaled <-melt(model.scaled, id.vars = "x")
    melt.model.scaled$x10 <- 10^melt.model.scaled$x
    
    vec.df <- data.frame(value = vec)
    vec10.df <- data.frame(value = 10^vec)
    
    melt.model.scaled$variable<- factor(melt.model.scaled$variable, as.character(1:length(melt.model.scaled$variable)))
    #(melt.model.scaled$variable)
    mix.model.plot <- ggplot(melt.model.scaled, aes(x = x10, y = value, fill = variable)) + 
      geom_histogram(data = vec10.df, aes (x = value, y = ..density..),
                     inherit.aes = F, bins = 100, col = "black", fill = NA) + 
      geom_area(alpha = 0.5, col = "grey50") + 
      
      scale_x_log10(breaks = major.ticks,
                    labels = major.ticks,
                    minor_breaks = minor.ticks,
                    name = channel) + 
      coord_cartesian(xlim= quantile(vec10.df$value, c(0.0005,0.9995))) +  
      ylab("Density") + 
      scale_fill_discrete(name = "Population") + 
      theme_classic()
    
    
    
    #ggsave(file = "plot.svg", mix.model.plot, width = 6, height = 4, units = "in")
    #mix.model.plot
    #####
    
    
    classif<- apply((probs.scaled.norm), 1, which.max)
    classif.uc <- apply(probs.scaled.norm, 1, max)

    #classif.uc
    #plot(classif.uc ~ vec)
    classif <- match(classif, rev(order(Snorm.analysis$mu)))
    #classif
    classif[classif.uc < (1- uccutoff)] <- 0 #
    #nrow(mydf2)
    #table(is.na(classif))
    #length(log10(mydf2[classif,channel]))
    #split(log10(mydf2[classif,channel]), classif)
    
    mydf$bc1 <-  classif
    
    #don't need to return classif
    debarcoded_data <- list('df' = mydf, "df2" = mydf2, "channel" = channel,
                            'snorm' = Snorm.analysis, "classif" = classif, "plot" = mix.model.plot,
                            'rglobj' = regression.output[["rglobj"]],
                            'regressionmodel' = regression.output[["coefs"]])
    
    return(debarcoded_data)
    
    
  }
  
  if (opt == 5) {
    vec <- log10(mydf2[,channel])
    Norm.analysis <- smsn.mix(vec, nu = 3, g = levels,
                               get.init = TRUE, criteria = TRUE, 
                               group = TRUE, family = "Skewt", calc.im = FALSE, obs.prob = TRUE)
    
    mix.hist(log10(mydf2[,channel]),Norm.analysis)
    classif <- Snorm.analysis$group
    classif.uc <- apply(Snorm.analysis$obs.prob, 1, max)
    classif <- match(classif, rev(order(Snorm.analysis$mu)))
    classif[classif.uc < (1- uccutoff)] <- 0 #
  }
  #qplot(classif.uc)

  #plot(mix.model, what = "classification")
  
  #plot(mix.model)
  #mix.model <- Mclust(mydf2[,channel], G = levels + 1)
  #plot(mix.model)
  #classif <- mix.model$classification
  #levels + 1
  #as.numeric(classi)
  #classif <- levels + 1 - classif #reverses order of levels
  #classif[classif.uc < (1- uccutoff)] <- 0 #
  print("Mixture modeling complete")
  

  ggplot(mydf, aes_string(y = paste0("`", channel, "`"), x =   paste0("`","FSC-A","`"))) + 
    scale_y_continuous(trans = logicle_trans()) + 
    geom_hex(bins = 100) + ggtitle("PO Assignments") + scale_fill_viridis()


  ggplot(mydf[classif != 0,], aes(y = `Pacific Orange-A`, x = `SSC-A`)) + 
    scale_y_continuous(trans = logicle_trans()) + 
    geom_point(data = mydf[classif == 0,], col = "black", alpha = 0.2) + 
    geom_point(aes(col = as.factor(classif[classif != 0])), alpha = 0.2) +
    ggtitle("PO Assignments")
  
  ggsave(paste0(channel,".jpg"),
         ggplot(mydf, aes(y = `Pacific Orange-A`, x = `SSC-A`)) + 
           scale_y_continuous(trans = logicle_trans()) + 
           geom_point(aes(col = as.factor(classif)), alpha = 0.5) + ggtitle("PO Assignments")
  )
  
  ggsave(paste0(channel,"2.jpg"),
         ggplot(mydf, aes(y = `Pacific Orange-A`, x = `FSC-A`)) + 
           scale_y_continuous(trans = logicle_trans()) + 
           geom_point(aes(col = as.factor(classif)), alpha = 0.5) + ggtitle("PO Assignments")
  )
  #table(classif)
  mydf$bc1 <-  classif
  
  
  
  #modified return value
  return(mydf)
}


debarcode.2 <- function(mydf, prevchannel, channel, levels, uccutoff = 0.05, subsample = 10e3) {
  print('running db2')
  #previous channel should always be bc1
  prevlevel <- "bc1"
  
  mydf$bc2 <- 0
  bc1table <- table(mydf$bc1) 
  
  bc1table.clean <- as.numeric(names(bc1table)[!(names(bc1table) == 0)])

  #if level = 1, exit function assign all bc2 to 1 (for single level barcodes)
  if (levels == 1) {
    mydf["bc2"] <- 1
    return(mydf)
  }
  
  for ( i in bc1table.clean) {

    print(paste("Debarcoding Level:", i))
    ind <- mydf[,prevlevel] == i

    mydf.i <- mydf[which(ind), ]
    
    Y <- log10(mydf.i[,channel])
    #mydf.i
    X <- log10(mydf.i[,prevchannel])

    resids <- Y - X + median(Y, na.rm = TRUE)


    if(length(resids) > subsample) {
      resids.ss <- sample(resids, subsample)
    } else {
      resids.ss <- resids
    }
    
    mod.int <- classIntervals(resids.ss, levels, style = "fisher")
    
    classif <- lapply(resids, FUN = function(x) {findInterval(x, mod.int$brks)})
    
    classif <- unlist(classif)

    classif <- levels + 1 - classif

    mydf[which(ind), "bc2"] <- classif
    #ggsave(paste0("PB_PO_",i,".jpg"),
    #)
    
    #should this be passed????
    #ggplot(data.frame(Y= Y, X = X, classif = as.factor(classif)), aes(y= Y, x = X)) + 
    #  geom_point(aes(col = classif)) + ggtitle(paste("PB Assignments", i))
    
  }
  return(mydf)
}
