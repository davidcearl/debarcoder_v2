minor.ticks <- (rep(2:9, 6))*(rep(10^(0:5), each = 8))
major.ticks <- 10^(1:10)


gatein <- function(data, gate) { #this function subsets the data within the specified gate
  flowCore::Subset(data, flowCore::filter(data, gate))
}

plotgate <- function(what, cells, gate, xlim = c(0,250e3), ylim = c(0, 250e3), ..., transforms = NULL) {
  lattice::xyplot(what, data=cells,
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
