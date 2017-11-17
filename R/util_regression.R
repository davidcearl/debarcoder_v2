## Selecting Densest Area
#----------------------------------------
selectDenseScatterArea <- function (data, fsc = 'FSC-A', ssc = 'SSC-A',
                                    subsample = 10e3) {

  if(!is.null(subsample) & (nrow(data) > subsample)){
    data <- data[sample(1:nrow(data), subsample),]
  }

  fsc_limits <- c(min(data[fsc]), max(data[fsc]))
  ssc_limits <- c(min(data[ssc]), max(data[ssc]))
  S <- 8
  area <- 0.95
  FGA <- matrix(0, 2^S, 2^S)
  AP <- matrix(1, 2^S, 2^S)
  N <- nrow(data)
  P <- MASS::kde2d(data[,fsc], data[,ssc], n = 2^S, lims = c(fsc_limits, ssc_limits))
  x_d <- matrix(rep(P[[1]],2^S), ncol = 2^S, nrow = 2^S, byrow = TRUE)
  y_d <- matrix(rep(P[[2]],2^S), ncol = 2^S, nrow = 2^S, byrow = FALSE)
  z_d <- P[[3]]
  z_d <- z_d/sum(z_d)
  ft <- function(x) abs(sum(z_d[z_d>x])- area)
  tmin <- optimize(ft, c(0,1e-4), tol = 1e-9)$objective
  AP <- AP*(z_d<tmin)
  FGA <- FGA + z_d
  Y1 <- x_d[!AP]
  Y2 <- y_d[!AP]
  Loc <- cbind(Y1, Y2)
  c <- FGA[!AP]

  return(list("loc" = Loc , "c" = c))

}

#Performing the Constrained regression
#----------------------------------------

#vars (1= FSC, 2= SSC, 3= FSC*SSC)
#powers (indicates power to which each predictor is raised)
regression_model_matrix <- function(nrow = 2, ncol = 8,
                                    vars = c(1, 1, 1, 2, 2, 2, 3, 3),
                                    powers = c(1, 0.5, 2, 1, 0.5, 2, 1, 0.5)){

    rm <- matrix(nrow = nrow, ncol = ncol)
    rm[1,] <- vars
    rm[2,] <- powers

    return(rm)
}

generate_regressors <- function(D, rm){

    Y <- D[,3]
    D[,3] <- D[,1]*D[,2]
    X <- matrix(0, nrow = length(Y), ncol = ncol(rm))

    for (m in 1:ncol(rm)){
        X[,m] = D[,rm[1,m]]^rm[2,m]
    }

    return(list('Y' = Y, 'D' = D, 'X' = X))
}

constrained_regression <- function(X, Y, fsc_limits, ssc_limits, val3, D,
                                   OFFSET, Bx, rm, monodir){
    S = 5
    area = 0.5
    xlimr <- c(floor(min(fsc_limits)), ceiling(max(fsc_limits)))
    ylimr <- c(floor(min(ssc_limits)), ceiling(max(ssc_limits)))
    P <- MASS::kde2d(D[,1], D[,2], n=2^S, lims= c(xlimr, ylimr))
    XX <- matrix(rep(P[[1]],2^S), ncol = 2^S, nrow = 2^S, byrow = TRUE)
    YY <- matrix(rep(P[[2]],2^S), ncol = 2^S, nrow = 2^S, byrow = FALSE)
    PP <- P[[3]]
    PP <- PP/sum(PP)

    ft <- function(t) {abs(sum(PP[PP>t])- area)}
    tmin <- optimize(ft, c(0,max(PP)))$minimum

    ZZ= matrix(nrow=sum(PP>tmin), ncol=4)
    ZZ[,1] <- as.numeric(XX[PP>tmin])
    ZZ[,2] <- as.numeric(YY[PP>tmin])
    ZZ[,3] <- ZZ[,1]*ZZ[,2]
    ZZ[,4] <- OFFSET

    for (i in 1:ncol(rm)){
        ZZ[,4] <- ZZ[,4] + Bx[i]*ZZ[,rm[1,i]]^rm[2,i]
    }

    dir1 <- lm(ZZ[,4] ~ ZZ[,1])
    dir2 <- lm(ZZ[,4] ~ ZZ[,2])

    if(is.null(monodir)){
        monodir <- +c(dir1$coefficients[2]>=0, dir2$coefficients[2]>=0)
    } else {
        monodir <- monodir
    }

    #2. Select points to evaluate on monotonicity
    S= 3
    P <- MASS::kde2d(D[,1], D[,2], n=2^S, lims= c(xlimr, ylimr))
    XX <- matrix(rep(P[[1]],2^S), ncol = 2^S, nrow = 2^S, byrow = TRUE)
    YY <- matrix(rep(P[[2]],2^S), ncol = 2^S, nrow = 2^S, byrow = FALSE)
    PP <- P[[3]]

    noc <- 2*(2^S-1)*2^S
    Q <- matrix(nrow=noc, ncol = 4)
    p <- 0
    for (i in 1:2^S) {
        for (j in 1:(2^S-1)) {
            p <- p + 1
            if (monodir[1] == 1) {
                Q[p,] <- c(XX[1,j], YY[i,1], XX[1,j+1], YY[i,1])
            } else {
                Q[p,] <- c(XX[1, j+1], YY[i, 1], XX[1, j], YY[i,1])
            }
        }
    }

    for (i in 1:2^S) {
        for (j in 1:(2^S-1)) {
            p <- p + 1
            if (monodir[2] == 1) {
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


    if (!is.null(val3)){
        A <- matrix(nrow = noc + 2, ncol = ncol(rm))
        A[1:noc,] <- A1-A2

        XBP <- matrix(nrow = 2, ncol=3)
        if (monodir[1]==1 & monodir[2]==1){
            XBP[,1]<- XX[1, c(1, 2^S)]
            XBP[,2]<- YY[c(1, 2^S), 1]
            "FSC up, SSC up"
        } else if (monodir[1]==0 & monodir[2]==1){
            XBP[,1]<- XX[1, c(2^S,1)]
            XBP[,2]<- YY[c(1, 2^S), 1]
            "FSC down, SSC up"
        } else if (monodir[1]==1 & monodir[2]==0){
            XBP[,1]<- XX[1, c(1, 2^S)]
            XBP[,2]<- YY[c(2^S, 1), 1]
            "FSC up, SSC down"
        } else if (monodir[1]==0 & monodir[2]==0){
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
    Bc <- pracma::lsqlincon(cbind(X, 1),Y,cbind(A,0),b)
    NOFFSET <- Bc[9]
    NB <- Bc[1:8]


    return(list('NB' = NB, 'Bc' = Bc, 'NOFFSET' = NOFFSET))
}

###############################################################################

doRegressContrained <- function(single_level_bc, fcb_df = NULL, Loc, weight, trans,
                                val3 = NULL, constrained_flag = 1,
                                columns = NULL, monodir = NULL, cofactor = NULL) {
  #print(cofactor)
  lo <- 1 #log offset
  fsc_limits <- c(min(single_level_bc['FSC-A']), max(single_level_bc['FSC-A']))
  ssc_limits <- c(min(single_level_bc['SSC-A']), max(single_level_bc['SSC-A']))
  data.corr.ls <- vector("list", length = length(columns))

  for (n in seq_along(columns)){
    ind <- append(c("FSC-A", "SSC-A"), columns[n])
    D <- single_level_bc[,ind]
    if (is.null(fcb_df)) {
      D2 <- D
    } else {
      D2 <- fcb_df[,ind]
    }

    rm <- regression_model_matrix()

    switch (trans,
            none = {
              D <- D
              D2 <- D2
            },
            logF = {
              D <- cbind(D[1:2], log(D[3] + lo))
              D2 <- cbind(D2[1:2], log(D2[3] + lo))
            },
            arcsinh = {
              D <- cbind(D[1:2], asinh(D[3]/cofactor))
              D2 <- cbind(D2[1:2], asinh(D2[3]/cofactor))
            }
    )

    #generate regressors
    regressors_1 <- generate_regressors(D, rm)
    regressors_2 <- generate_regressors(D2, rm)

    D <- regressors_1[['D']]
    D2 <- regressors_2[['D']]

    Y <- regressors_1[['Y']]
    Y2 <- regressors_2[['Y']]

    X <- regressors_1[['X']]
    X2 <- regressors_2[['X']]

    #intial unconstrained regression
    B <- lm(Y ~ X)
    OFFSET <- B$coefficients[1]
    Bx <- B$coefficients[2:9]
    mod.resid <- as.numeric(B$residuals)

    if (constrained_flag == 1) {
        con_reg <- constrained_regression(X, Y, fsc_limits, ssc_limits, val3, D,
                                          OFFSET, Bx, rm, monodir)
        Bc <- con_reg[['Bc']]
        NOFFSET <- con_reg[['NOFFSET']]
        NB <- con_reg[['NB']]
        Nresidual <- Y - (X%*%Bc[1:8] + NOFFSET)
        Nresidual2 <- Y2 - (X2%*%Bc[1:8] + NOFFSET)
    }
    else {
      NB <- Bx
      NOFFSET <- OFFSET
      Nresidual <- mod.resid
    }

    XO <- matrix(nrow= length(weight), ncol = 4)
    XO[,1] <- Loc[,1]
    XO[,2] <- Loc[,2]
    XO[,3] <- Loc[,1]*Loc[,2]
    XO[,4] <- NOFFSET

    for (i in 1:ncol(rm)){
      XO[,4] <- XO[,4] + NB[i]*XO[,rm[1,i]]^rm[2,i]
    }

    R2 <- as.numeric(Nresidual2) + as.vector(t(weight)%*%XO[,4]/sum(weight))

    if (trans == 'logF') {
        data.corr <- exp(R2) - lo
    } else if (trans == 'none') {
        data.corr <- R2
    } else if (trans == 'arcsinh') {
        data.corr <- sinh(R2)*cofactor
    }
  

    print(paste(columns[n], "complete"))
    data.corr.ls[[n]] <- data.corr
  }
  data.corr.df <- as.data.frame(data.corr.ls)
  colnames(data.corr.df) <- columns

  return(list(df = data.corr.df, coefs = Bc))
}
