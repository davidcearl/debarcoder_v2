library(pracma)
#david earl's regression correction
deCorr <- function(mydata, opt = 'log') {
  
  lo <- 1
  if (opt == 'lin'){
    mydata <- mydata
  } else if (opt == 'log') {
    mydata<- log(mydata + lo)
  } else if (opt == 'logF') {
    mydata$`Pacific Orange-A` <- log(mydata$`Pacific Orange-A` + lo)
    mydata$`Pacific Blue-A` <- log(mydata$`Pacific Blue-A` + lo)
  }
  
  fscA_mean <- mean(mydata$`FSC-A`)
  sscA_mean <- mean(mydata$`SSC-A`)
  event_mean <- fscA_mean + sscA_mean
  
  mydata$`Pacific Orange-A` <- (event_mean/(mydata$`FSC-A` + mydata$`SSC-A`))*mydata$`Pacific Orange-A`
  mydata$`Pacific Blue-A` <- (event_mean/(mydata$`FSC-A` + mydata$`SSC-A`))*mydata$`Pacific Blue-A`
  
  if (opt == 'lin'){
    return(mydata)
  } else if (opt == 'log') {
    mydata <- exp(mydata) + lo
    return (mydata)
  } else if (opt == 'logF') {
    mydata$`Pacific Orange-A` <- exp(mydata$`Pacific Orange-A`) - lo
    mydata$`Pacific Blue-A` <- exp(mydata$`Pacific Blue-A`) - lo
    return(mydata)
  }
}

## Implementation of morphology correction from 


## PreProccessing
fcPreprocess <- function (Data.ff, t.trim = 0.2) {
  if(!class(Data.ff) == "flowFrame"){
    stop("Argument must be Data.frame")
  }
  
  Data.df <- as.data.frame(exprs(Data.ff))
  as.numeric(description(Data.ff)$`$TIMESTEP`)
  p<- t.trim / as.numeric(description(Data.ff)$`$TIMESTEP`)
  Data.df <- Data.df[Data.df$Time > min(Data.df$Time) + p, ]#remove first .2 seconds
  Data.df <- Data.df[Data.df$Time < max(Data.df$Time) - p, ]#remove last .2 seconds
  
  #exlcude data with minimal or maximal values
  maxs <- (sapply(Data.df, max))
  mins <- (sapply(Data.df, min))
  Data.df2 <- Data.df[(apply(mapply('<', Data.df, maxs), 1, all)),]
  Data.df2 <- Data.df2[(apply(mapply('>', Data.df2, mins), 1, all)),]
  nrow(Data.df2) / nrow(exprs(Data.ff))
  
  #excude top and bottom 5% of FSC/SSC data
  Data.df3 <- Data.df2[Data.df2$`FSC-A` > quantile(Data.df2$`FSC-A`, 0.05),]
  Data.df3 <- Data.df3[Data.df3$`FSC-A` < quantile(Data.df2$`FSC-A`, 0.95),]
  Data.df3 <- Data.df3[Data.df3$`SSC-A` > quantile(Data.df2$`SSC-A`, 0.05),]
  Data.df3 <- Data.df3[Data.df3$`SSC-A` < quantile(Data.df2$`SSC-A`, 0.95),]
  
  #exlucde negative values for Barcoding channels
  Data.df4 <- Data.df3[Data.df3$`Pacific Orange-A` > 0,]
  Data.df4 <- Data.df4[Data.df4$`Pacific Blue-A` > 0,]
  
  return(Data.df4)
}


## Selecting Densest Area
#----------------------------------------
selectDenseScatterArea <- function (Data, fsc = 'FSC-A', ssc = 'SSC-A', subsample = 10e3) {
  if(!is.null(subsample) & (nrow(Data) > subsample)){
    Data<- Data[sample(1:nrow(Data), subsample),]
  }
  val1 <- c(min(Data[fsc]), max(Data[fsc]))
  val2 <- c(min(Data[ssc]), max(Data[ssc]))
  S <- 8
  area <- 0.95
  FGA <- matrix(0, 2^S, 2^S)
  AP <- matrix(1, 2^S, 2^S)   
  N <- nrow(data)
  
  P <-kde2d(Data[,fsc], Data[,ssc], n=2^S, lims= c(val1, val2))
  
  x.d <- matrix(rep(P[[1]],2^S), ncol = 2^S, nrow = 2^S, byrow = TRUE)
  y.d <- matrix(rep(P[[2]],2^S), ncol = 2^S, nrow = 2^S, byrow = FALSE)
  z.d <- P[[3]]
  z.d <- z.d/sum(z.d)
  image(x.d[1,], y.d[,1], z.d, xlab= fsc, ylab = ssc, main = "Cell Density", xlim = val1, ylim = val2)
  
  ft <- function(x) abs(sum(z.d[z.d>x])- area)
  
  
  #pracma::fminsearch(ft, 1/(2^S)^2)
  #optim((1/(2^S)^2), ft)
  # guess = 1/(2^S)^2
  # optimize(ft, c(0,guess), tol = 1e-9)
  # optimize(ft, c(0,10*guess), tol = 1e-9)
  tmin <- optimize(ft, c(0,1e-4), tol = 1e-9)$objective
  #image(x.d[1,], y.d[,1], z.d>tmin, xlab= 'FSC-A', ylab = 'SSC-A', main = "95% of density")
  
  AP <- AP*(z.d<tmin)
  FGA <- FGA + z.d
  
  Y1 <- x.d[!AP]
  Y2 <- y.d[!AP]
  Loc <- cbind(Y1, Y2)
  c <- FGA[!AP]
  
  
  return(list("loc" = Loc , "c" = c))
}

#Perfomring the Constrained regression
#----------------------------------------



doRegressContrained <- function(mydata, mydata2 = NULL, Loc, Weight, opt, val3 = NULL, constrained_flag = 1, rm = NULL, plot_flag = 1, columns = NULL, MonoDir = NULL) {
  #channel <- "Pacific Orange-A"
  # columns <-c(channel)
  #  opt <- 'logF'
  #  MonoDir <- c(1,1)
  #  constrained_flag <- 1
  #  val3 <- NULL
  #  plot_flag <- 1
  #   Weight <- k.da$c
  #   Loc <- k.da$loc
  lo<- 1 #log ofset
  val1 <- c(min(mydata['FSC-A']), max(mydata['FSC-A']))
  val2 <- c(min(mydata['SSC-A']), max(mydata['SSC-A']))
  data.corr.ls <- vector("list", length = length(columns))
  #names(data.corr.ls) <- columns
  #n <- 1
  for (n in seq_along(columns)){
    ind <- append(c("FSC-A", "SSC-A"), columns[n])
    D <- mydata[,ind]
    if (is.null(mydata2)) {
      D2 <- D
    } else {
      D2 <- mydata2[,ind]
    }
    nrow(D2)
    nrow(D)
    #summary(D)
    rm <- matrix(nrow=2, ncol =8)
    rm[1,] <- c(1,1,1,2,2,2,3,3) #regression model (1= FSC, 2= SSC, 3= FSC*SSC)
    rm[2,] <- c(1, 0.5, 2, 1, 0.5, 2, 1, 0.5) #regression model (indicates power to which each predictor is raised)
    switch (opt,
            lin = {
              D <- D
              D2 <- D2
            },
            log = {
              Loc <- log(Loc + lo)
              val1 <- log(val1 + lo)
              val2 <- log(val2 + lo)
              
              D <- log(D+lo)
              D2 <- log(D2+lo)
            },
            logF = {
              D <- cbind(D[1:2], log(D[3]+lo))
              D2 <- cbind(D2[1:2], log(D2[3]+lo))
            }
    )
    #generate regressors
    Y <- D[,3]
    Y2 <- D2[,3]
    #summary(Y2)
    D[,3] <- D[,1]*D[,2]
    D2[,3] <- D2[,1]*D2[,2]
    
    X <- matrix(0, nrow= length(Y), ncol = ncol(rm))  
    X2 <-  matrix(0, nrow= length(Y2), ncol = ncol(rm)) 
    for (m in 1:ncol(rm)){
      X[,m] = D[,rm[1,m]]^rm[2,m]
    }
    
    for (m in 1:ncol(rm)){
      X2[,m] = D2[,rm[1,m]]^rm[2,m]
    }
    #unconstrained regression
    B <- lm(Y ~ X)
    OFFSET <- B$coefficients[1]
    Bx <- B$coefficients[2:9]
    mod.resid <- as.numeric(B$residuals) 
    
    if (constrained_flag == 1) {
      S = 5
      area = 0.5
      xlimr <- c(floor(min(val1)), ceiling(max(val1)))
      ylimr <- c(floor(min(val2)), ceiling(max(val2)))
      P <-kde2d(D[,1], D[,2], n=2^S, lims= c(xlimr, ylimr))
      XX <- matrix(rep(P[[1]],2^S), ncol = 2^S, nrow = 2^S, byrow = TRUE)
      YY <- matrix(rep(P[[2]],2^S), ncol = 2^S, nrow = 2^S, byrow = FALSE)
      PP <- P[[3]]
      PP <- PP/sum(PP)
      
      ft <- function(t) {abs(sum(PP[PP>t])- area)}
      tmin <- optimize(ft, c(0,max(PP)))$minimum
      
      ZZ= matrix(nrow=sum(PP>tmin), ncol=4)
      ZZ[,1] <- XX[PP>tmin]
      ZZ[,2] <- YY[PP>tmin]
      ZZ[,3] <- ZZ[,1]*ZZ[,2]
      ZZ[,4] <- OFFSET
      for (i in 1:ncol(rm)){
        ZZ[,4] <- ZZ[,4] + Bx[i]*ZZ[,rm[1,i]]^rm[2,i]
      }
      
      dir1 <- lm(ZZ[,4] ~ ZZ[,1])
      dir2 <- lm(ZZ[,4] ~ ZZ[,2])
      
      if(is.null(MonoDir)){
        MonoDir <- +c(dir1$coefficients[2]>=0, dir2$coefficients[2]>=0)
      } else {
        MonoDir <- MonoDir
      }
      #2. Select points to evaluate on monotonicity
      S= 3
      P <-kde2d(D[,1], D[,2], n=2^S, lims= c(xlimr, ylimr))
      XX <- matrix(rep(P[[1]],2^S), ncol = 2^S, nrow = 2^S, byrow = TRUE)
      YY <- matrix(rep(P[[2]],2^S), ncol = 2^S, nrow = 2^S, byrow = FALSE)
      PP <- P[[3]]
      
      noc <- 2*(2^S-1)*2^S
      Q <- matrix(nrow=noc, ncol = 4)
      p <- 0
      for (i in 1:2^S) {
        for (j in 1:(2^S-1)) {
          p <- p + 1
          if (MonoDir[1] == 1) {
            Q[p,] <- c(XX[1,j], YY[i,1], XX[1,j+1], YY[i,1])
          } else {
            Q[p,] <- c(XX[1, j+1], YY[i, 1], XX[1, j], YY[i,1])
          }
        }
      }
      
      for (i in 1:2^S) {
        for (j in 1:(2^S-1)) {
          p <- p + 1
          if (MonoDir[2] == 1) {
            Q[p,] <- c(XX[1,i], YY[j,1], XX[1,i], YY[j+1,1])
          } else {
            Q[p,] <- c(XX[1, i], YY[j+1, 1], XX[1,i], YY[j,1])
          }
        }
      }
      
      XX1 <- matrix(nrow = noc, ncol = 3)
      XX1[,1] <- Q[,1]
      XX1[,2] <- Q[,2]
      XX1[,3] <- XX1[,1]*XX1[,2]
      A1 <- matrix(nrow=noc, ncol = ncol(rm))
      for (m in 1:ncol(rm)) {
        A1[,m] <- XX1[,rm[1,m]]^rm[2,m]
      }
      XX2 <- matrix(nrow = noc, ncol = 3)
      XX2[,1] <- Q[,3]
      XX2[,2] <- Q[,4]
      XX2[,3] <- XX2[,1]*XX2[,2]
      A2 <- matrix(nrow=noc, ncol = ncol(rm))
      for (m in 1:ncol(rm)) {
        A2[,m] <- XX2[,rm[1,m]]^rm[2,m]
      }
      !is.null(val3)
      if (!is.null(val3)){
        A <- matrix(nrow = noc + 2, ncol = ncol(rm))
        A[1:noc,] <- A1-A2
        
        XBP <- matrix(nrow = 2, ncol=3)
        if (MonoDir[1]==1 & MonoDir[2]==1){
          XBP[,1]<- XX[1, c(1, 2^S)]
          XBP[,2]<- YY[c(1, 2^S), 1]
          "FSC up, SSC up"
        } else if (MonoDir[1]==0 & MonoDir[2]==1){
          XBP[,1]<- XX[1, c(2^S,1)]
          XBP[,2]<- YY[c(1, 2^S), 1]
          "FSC down, SSC up"
        } else if (MonoDir[1]==1 & MonoDir[2]==0){
          XBP[,1]<- XX[1, c(1, 2^S)]
          XBP[,2]<- YY[c(2^S, 1), 1]
          "FSC up, SSC down"
        } else if (MonoDir[1]==0 & MonoDir[2]==0){
          XBP[,1]<- XX[1, c(2^S, 1)]
          XBP[,2]<- YY[c(2^S, 1), 1]
          "FSC down, SSC down"
        }
        
        XBP[,3] <- XBP[,1]*XBP[,2]
        AX <- matrix(nrow = 2, ncol= ncol(rm))
        for (m in 1:ncol(rm)){
          AX[,m] <- XBP[,rm[1,m]]^rm[2,m]
        }
        
        A[noc+1, ] <- -AX[1,]
        A[noc+2, ] <- AX[2,]
        
        b <- matrix(0, nrow=noc+2, 1)
        b(noc+1) <- val3[1]
        b(noc+2) <- val3[2]
        
      } else {
        A <- A1-A2
        b <- matrix(0, nrow=noc, ncol=1)
      }
      
      #constrained regression analysis
      NB <- numeric(length = length(Bx))
      Bc <- lsqlincon(cbind(X, 1),Y,cbind(A,0),b)
      NOFFSET <- Bc[9]
      NB <- Bc[1:8]
      Nresidual <- Y - (X%*%Bc[1:8] + NOFFSET)
      Nresidual2 <- Y2 - (X2%*%Bc[1:8] + NOFFSET)
    } else {
      NB <- Bx
      NOFFSET <- OFFSET
      Nresidual <- mod.resid
    }
    
    
    if (plot_flag == 1) {
      my_surface <- function(f, n=10, ...) { 
        x <- seq(min(x1), max(x1), length=n)
        y <- seq(min(x2), max(x2), length=n)
        z <- outer(x,y,f)
        surface3d(x, y, z, ...)
      }
      library(rgl)
      fc <- function(x1, x2) {
        NOFFSET + NB[1]*x1 + NB[2]*x1^0.5 + NB[3]*x1^2 + NB[4]*x2 + NB[5]*x2^0.5 + NB[6]*x2^2 + NB[7]*(x1*x2) + NB[8]*(x1*x2)^0.5
      }
      
      f <- function(x1, x2) {
        OFFSET + Bx[1]*x1 + Bx[2]*x1^0.5 + Bx[3]*x1^2 + Bx[4]*x2 + Bx[5]*x2^0.5 + Bx[6]*x2^2 + Bx[7]*(x1*x2) + Bx[8]*(x1*x2)^0.5
      }
      
      sampler <- sample(nrow(X),1e4, replace= TRUE)
      X.downsampled <- X2[sampler,]
      x1 <- X.downsampled[,1]
      x2 <- X.downsampled[,4]
      y <- Y2[sampler]
      plot3d(x1,x2,y, type="p", col="red", xlab="FSC", ylab="SSC", zlab="Dye", site=5, lwd=1)
      my_surface(fc, alpha=.2, col = "blue")
      my_surface(f, alpha=.2, col = "green")
    }
    
    XO <- matrix(nrow= length(Weight), ncol = 4)
    XO[,1] <- Loc[,1]
    XO[,2] <- Loc[,2]
    XO[,3] <- Loc[,1]*Loc[,2]
    XO[,4] <- NOFFSET
    
    for (i in 1:ncol(rm)){
      XO[,4] <- XO[,4] + NB[i]*XO[,rm[1,i]]^rm[2,i]
    }

    R1 <- as.numeric(Nresidual) + as.vector(t(Weight)%*%XO[,4]/sum(Weight))
    R2 <- as.numeric(Nresidual2) + as.vector(t(Weight)%*%XO[,4]/sum(Weight))


    if(opt == 'log') {
      data.corr <- exp(R2) - lo
    } else if (opt == 'logF') {
      data.corr <- exp(R2) - lo
    } else if (opt == 'lin') {
      data.corr <- R2
    }
    print(paste(columns[n], "complete"))
    data.corr.ls[[n]] <- data.corr
  }
  data.corr.df <- as.data.frame(data.corr.ls)
  colnames(data.corr.df) <- columns
  
  if (plot_flag == 0){
    my_surface <- function(f, n=10, ...) { 
      x <- seq(min(x1), max(x1), length=n)
      y <- seq(min(x2), max(x2), length=n)
      z <- outer(x,y,f)
      surface3d(x, y, z, ...)
    }
    #library(rgl)
    fc <- function(x1, x2) {
      NOFFSET + NB[1]*x1 + NB[2]*x1^0.5 + NB[3]*x1^2 + NB[4]*x2 + NB[5]*x2^0.5 + NB[6]*x2^2 + NB[7]*(x1*x2) + NB[8]*(x1*x2)^0.5
    }
    
    f <- function(x1, x2) {
      OFFSET + Bx[1]*x1 + Bx[2]*x1^0.5 + Bx[3]*x1^2 + Bx[4]*x2 + Bx[5]*x2^0.5 + Bx[6]*x2^2 + Bx[7]*(x1*x2) + Bx[8]*(x1*x2)^0.5
    }
    
    sampler <- sample(nrow(X),1e4, replace= TRUE)
    X.downsampled <- X2[sampler,]
    x1 <- X.downsampled[,1]
    x2 <- X.downsampled[,4]
    y <- Y2[sampler]
    #plot3d(x1,x2,y, type="p", col="red", xlab="FSC", ylab="SSC", zlab="Dye", site=5, lwd=1, main = "regression w/ monotonicity constraint")
    #my_surface(fc, alpha=.2, col = "blue")
    #my_surface(f, alpha=.2, col = "green")
  }
  return(list(df = data.corr.df, rglobj = list(x1 = x1, x2 = x2, y = y, fc = fc, f = f), coefs = Bc))
}
