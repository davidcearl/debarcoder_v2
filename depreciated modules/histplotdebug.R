list.files()
model <- read.csv("model.csv")[-1]
vec.bt.df <-read.csv("vecbtdf.csv")[-1]
str(vec.bt.df)
snorm.pii <-unlist(read.csv("snormpii.csv")[-1])
vec <- read.csv("vec.csv")[-1]

ggplot2::qplot(vec, bins = 200)

model.scaled <- apply(model[,-1], 1, function(vec) { vec * snorm.pii})
#model.scaled <- model.scaled/max(model.scaled)
model.scaled<- cbind(x = model[,1], as.data.frame(t(as.matrix(model.scaled))))
num.colnames <- as.numeric(colnames(model.scaled[-1]))
colnames.order <- rev(order(Snorm.df$loc))

colnames(model.scaled)[-1] <- match(num.colnames, colnames.order)

melt.model.scaled <- reshape2::melt(model.scaled, id.vars = "x")

melt.model.scaled <- reshape2::melt(model.scaled, id.vars = "x")

library(ggplot2)
str(vec)
ggplot(melt.model.scaled, aes(x = x, y = value, fill = variable)) + 
  geom_histogram(data = vec, aes (x = x, y = ..density..), bins = 100, inherit.aes = F) + 
  geom_area()

mytrans <- "identity"
channel <- "Pacific O"
colnames(melt.model.scaled)
quantile(vec.bt.df$x)

ggplot2::ggplot(melt.model.scaled,
                ggplot2::aes_string(x = "x",
                                    y = "value",
                                    fill = "variable")) +
   ggplot2::geom_histogram(data = data.frame(x = vec), ggplot2::aes(x = x, y = ..density..),
                           inherit.aes = F, bins = 100, col = "black", fill = NA) +
   ggplot2::geom_area(alpha = 0.5, col = "grey50") +
   ggplot2::scale_x_continuous(trans = mytrans,
                               breaks = major.ticks,
                               labels = major.ticks,
                               minor_breaks = minor.ticks,
                               name = channel) +
   ggplot2::coord_cartesian(xlim= quantile(vec, c(0.0005,0.9995))) +
   ggplot2::ylab("Density") +
   ggplot2::scale_fill_discrete(name = "Population") +
   ggplot2::theme_classic()
#

trans <- "arcsinh"
vec.df <- data.frame(value = vec)
cofactor_bc1 <- 150

if(trans == "log10"){
  vec <- unname(unlist(vec))
  vec.bt.df <- data.frame(value = 10^vec) #bt = back transform
  melt.model.scaled$x10 <- 10^melt.model.scaled$x
  print(colnames(vec.bt.df))
  mytrans <- "log10"
} else if (trans =="arcsinh"){
  vec <- unname(unlist(vec))
  vec.bt.df <- data.frame(value = sinh(vec)*cofactor_bc1)
  print(colnames(vec.bt.df))
  melt.model.scaled$x10 <- sinh(melt.model.scaled$x)*cofactor_bc1
  mytrans <- asinh_trans(cofactor_bc1)
}
#vec.bt.df <- NULL

# melt.model.scaled$variable<- factor(melt.model.scaled$variable,
#                                     as.character(1:length(melt.model.scaled$variable)))

melt.mo

length(vec)

melt.model.scaled$pseudocount <-melt.model.scaled$value
binobject <- ggplot2::stat_bin(data = vec.bt.df, ggplot2::aes(x = value, y = ..density..),
                  inherit.aes = F, bins = 100, col = "black", fill = NA)

preplot <- ggplot(vec.bt.df, aes(x = value, y = ..density..)) + 
  stat_bin(bins = 100) + 
  ggplot2::scale_x_continuous(trans = mytrans,
                              breaks = major.ticks,
                              labels = major.ticks,
                              minor_breaks = minor.ticks) + 
  ggplot2::coord_cartesian(xlim= quantile(vec.bt.df$value, c(0.0005,0.9995)))
  
preplot.build <- ggplot_build(preplot)

melt.model.scaled$value_rescaled <- melt.model.scaled$value/
  max(melt.model.scaled$value)*max(preplot.build$data[[1]][,"y"]) 

ggplot2::ggplot(melt.model.scaled,
                ggplot2::aes_string(x = "x10",
                                    y = "pseudocount",
                                    fill = "variable")) +
  ggplot2::stat_bin(data = vec.bt.df, ggplot2::aes(x = value, y = ..density..),
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

ggplot2::ggplot(melt.model.scaled,
                ggplot2::aes_string(x = "x10",
                                    y = "value_rescaled",
                                    fill = "variable")) +
  ggplot2::stat_bin(data = vec.bt.df, ggplot2::aes(x = value, y = ..density..),
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


#(melt.model.scaled$variable)

