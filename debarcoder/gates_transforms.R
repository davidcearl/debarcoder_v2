
require(ggplot2)

minor.ticks <- (rep(2:9, 6))*(rep(10^(0:5), each = 8))
major.ticks <- 10^(1:10)


gatein <- function(data, gate) { #this function subsets the data within the specified gate
  Subset(data, filter(data, gate))
}
?Subset
plotgate <- function(what, cells, gate, xlim = c(0,250e3), ylim = c(0, 250e3), ..., transforms = NULL) {
  xyplot(what, data=cells,
         ylim=ylim, xlim=xlim, displayFilter=TRUE,
         filter=gate, smooth=F, stat=T, xbin = 256,
         pos=0.5, abs=T)
}

asinh_trans = function(cofactor) {
  transform =   function(x,...) {
    asinh(x/cofactor)/log(10) + log(cofactor/2)/log(10)
  }
  inverse = function(x,...) {
    sinh(x*log(10) - log(cofactor/2))*cofactor
  }
  breaks = function(x) {
    minor.ticks <- (rep(2:9, 6))*(rep(10^(0:5), each = 8))
    major.ticks <- 10^(2:6)
    breaks <-(c(0, major.ticks,- major.ticks))
    breaks[order(breaks)]
  }
  scales::trans_new("asinh", transform, inverse, breaks)
} 


logicle_trans = function(w=0.5, t=262144, m=4.5) {
  logicleTransform = flowCore::logicleTransform(w=w, t=t, m=m, a=0)
  transform = function(x) {
    logicleTransform(x)
  }
  inverse = function(x) {
    flowCore::inverseLogicleTransform("inverseLogicle", logicleTransform)(x)
  }
  breaks = function(x) {
    lim = 10**(1+w)
    linear = scales::pretty_breaks(n=3, min.n=3)(c(x[1], lim))
    logs = scales::log_breaks(10)(c(lim, x[2]))
    unique(c(linear[linear <= lim], logs[logs > lim]))
  }
  scales::trans_new("logicle", transform, inverse, breaks)
}

logicle_trans_rev = function(w=0.5, t=262144, m=4.5) {
  logicleTransform = flowCore::logicleTransform(w=w, t=t, m=m, a=0)
  transform = function(x) {
    -logicleTransform(x)
  }
  inverse = function(x) {
    -flowCore::inverseLogicleTransform("inverseLogicle", logicleTransform)(x)
  }
  breaks = function(x) {
    lim = 10**(1+w)
    linear = scales::pretty_breaks(n=3, min.n=3)(c(x[1], lim))
    logs = scales::log_breaks(10)(c(lim, x[2]))
    -unique(c(linear[linear <= lim], logs[logs > lim]))
  }
  scales::trans_new("logicle", transform, inverse, breaks)
}
remove_negatives <- function (mydf, channels) {
  #mydf is a dataframe
  #channels is a vector of channel names
  for (i in 1:length(channels)) {
    mydf <- mydf[mydf[channels[i]] >0,]
  }
  return(mydf)
}

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}


fancy_scientific <- function(l) {
  # turn in to character string in scientific notation
  l <- format(l, scientific = TRUE)
  # quote the part before the exponent to keep all the digits
  # turn the 'e+' into plotmath format
  l <- gsub("e", "%*%10^", l)
  l <- gsub("1%\\*%","", l)
  l <- gsub("\\+0", "", l)
  # return this as an expression
  parse(text=l)
}

stairstepn <- function( data, direction="hv", yvars="y" ) {
  direction <- match.arg( direction, c( "hv", "vh" ) )
  data <- as.data.frame( data )[ order( data$x ), ]
  n <- nrow( data )
  
  if ( direction == "vh" ) {
    xs <- rep( 1:n, each = 2 )[ -2 * n ]
    ys <- c( 1, rep( 2:n, each = 2 ) )
  } else {
    ys <- rep( 1:n, each = 2 )[ -2 * n ]
    xs <- c( 1, rep( 2:n, each = 2))
  }
  
  data.frame(
    x = data$x[ xs ]
    , data[ ys, yvars, drop=FALSE ]
    , data[ xs, setdiff( names( data ), c( "x", yvars ) ), drop=FALSE ]
  ) 
}

stat_stepribbon <- 
  function(mapping = NULL, data = NULL, geom = "ribbon", position = "identity", inherit.aes = TRUE) {
    ggplot2::layer(
      stat = Stepribbon, mapping = mapping, data = data, geom = geom, 
      position = position, inherit.aes = inherit.aes
    )
  }

StatStepribbon <- 
  ggproto("stepribbon", Stat,
          compute_group = function(., data, scales, direction = "hv", yvars = c( "ymin", "ymax" ), ...) {
            stairstepn( data = data, direction = direction, yvars = yvars )
          },                        
          required_aes = c( "x", "ymin", "ymax" )
  )
