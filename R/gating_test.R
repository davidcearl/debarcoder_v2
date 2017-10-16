library(flowCore)
library(CytobankAPI)

token <- "eyJ0eXAiOiJKV1QiLCJhbGciOiJIUzI1NiJ9.eyJqdGkiOiJjZGMyZWZhNTlhNzM1ZjY4MDRiMDRhOWI3ZmJlMjU4NSIsImV4cCI6MTUwODEyNjg4NywiYXVkIjoiY3l0b2JhbmtfYXBpX3YxX3VzZXJzIiwidXNlcl9pZCI6MTQ3LCJpYXQiOjE1MDgwOTgwODcsImlzcyI6Imh0dHBzOi8vdmFuZGVyYmlsdC5jeXRvYmFuay5vcmcvIiwibmJmIjoxNTA4MDk4MDg3LCJzdWIiOiJjeXRvYmFua19hcGlfdjEifQ.grig6yecT4QqkIeKATeZrsFPx6NHSFFxg4qQ2UV2mrU"
cyto_session <- authenticate(site="irishlab", username="benjamin.reisman@vanderbilt.edu", auth_token  = token)

experiment.id <- 24389
sample <- 349727
experiment.name <- unlist(experiments.show(cyto_session, experiment.id)$experimentName)
experiment.name

fcs_files.download(cyto_session, experiment.id, sample)
my.ff <- read.FCS("H3K18Ac_FCB Only_001.fcs")
my.ff <- compensate(my.ff, my.ff@description$SPILL)



#get exp population list from cytobank api
get_populations <- function(cyto_session, exp_id){
  return(CytobankAPI::populations.list(cyto_session, exp_id,  output = "default"))
}


#get exp gate list from cytobank api
get_gates <- function(cyto_session, exp_id){
  return(CytobankAPI::gates.list(cyto_session, exp_id, output = "default"))
}

get_compensations <- function(cyto_session, exp_id) {
  return(CytobankAPI::compensations.list(cyto_session, exp_id, output = 'default'))
}

#name lookup table for consistent naming?
# "Panel 1" is hardcoded
#needs exp_id from shiny ui
get_lut <- function(cyto_session, exp_id ) {
  scales <- CytobankAPI::scales.list(cyto_session,
                                     exp_id,
                                     output = "default")
  mypanel <- CytobankAPI::panels.list(cyto_session,
                                      exp_id,
                                      output = "default")[["Panel 1"]][["channels"]]
  scales_df <- as.data.frame(lapply(scales, function(X) unname(unlist(X))))
  mypanel_df <- as.data.frame(lapply(mypanel, function(X) unname(unlist(X))))
  lut <- merge.data.frame(mypanel_df,
                          scales_df,
                          by.x= "normalizedShortNameId",
                          by.y ="normalizedShortNameId")
  extrarow <- data.frame(-1, "Null", "Null", -1, 0, 0, 0, 0, 0, 1)
  colnames(extrarow) <- colnames(lut)
  lut <- rbind(lut, extrarow)

  return(lut)
}

#data.frame(-1, "Null", "Null", -1, 0, 0, 0, 0, 0, 0)
exp_pops <- get_populations(cyto_session, experiment.id)
exp_gates <- get_gates(cyto_session, experiment.id)
exp_comps <- get_compensations(cyto_session, experiment.id)
exp_lut <- get_lut(cyto_session, experiment.id)


#exp_gates$definition[[7]]
  
exp_info <- list('exp_id' = experiment.id,
                 'exp_pops' = exp_pops,
                 'exp_comps' = exp_comps,
                 'exp_gates' = exp_gates,
                 'exp_lut' = exp_lut)




gate_defs <- define_gates(exp_info$exp_gates, exp_info$exp_lut)

#gates.gatingML_download(cyto_session, experiment.id)

define_gates <- function(gates, lut) {
  #exp_info$exp_gates
  #gates <- exp_info$exp_gates
  #lut <- exp_info$exp_lut
  #gates$name
  mygates <- vector("list", length =  max(unlist(gates$gateId)))
  names(mygates)[unlist(gates$gateId)] <- unlist(gates$name)
  # gates$name[[4]]
  # i <- 4

  #mygates[[j]]
  for(i in 1:length(unlist(gates$gateId))) {
    j <- gates$gateId[[i]]
    #j
    #i <- 3
  #  gates$xNormalizedShortNameId[[i]]
    #gates$yNormalizedShortNameId[[i]]
  #  mygates[[j]]
    #as.character(lut[match(mygates[[i]][["channels"]],
     #                      lut$normalizedShortNameId),"shortName"])
    #mygates[[j]]
    #j <- 6
    mygates[[j]][["channels"]] <- c(gates$xNormalizedShortNameId[[i]], gates$yNormalizedShortNameId[[i]])
    mygates[[j]][["channels"]] <- as.character(lut[match(mygates[[j]][["channels"]],
                                                         lut$normalizedShortNameId),"shortName"])
    
    mygates[[j]][["type"]] <- gates$type[[i]]
    mygates[[j]]
    #mygates[[i]][["type"]]
    #print(str(mygates))
    if(mygates[[j]][["type"]] == "PolygonGate") {
      print("waspolygon")
      mygates[[j]][["coords"]] <- do.call(rbind,
                                          lapply(
                                            gates$definition[[i]][[1]][["polygon"]][["vertices"]],
                                            as.numeric))
      colnames(mygates[[j]][["coords"]]) <- mygates[[j]][["channels"]]
    } else if(mygates[[j]][["type"]] == "RectangleGate") {
      print("wasrecta")
      mygates[[j]][["coords"]] <- matrix(
        unlist(gates[["definition"]][[i]][[1]][["rectangle"]]), ncol = 2, byrow = TRUE)
      colnames(mygates[[j]][["coords"]]) <- mygates[[j]][["channels"]]
      
    } else if(mygates[[j]][["type"]] == "EllipseGate") {
      print("acquiring elipsegate")
      #procedure for computing ellipsegate parameters from : https://support.bioconductor.org/p/35360/
      
      mydef <- (gates[["definition"]][[i]])[[1]]
        angle <- mydef$ellipse$angle
        major <- mydef$ellipse$major/2
        minor <- mydef$ellipse$minor/2
        
        m1 <- cos(angle)^2/(major^2) + (sin(angle)^2)/(minor^2)
        m2 <- sin(angle)*cos(angle)*(((1/major^2)) - (1/(minor^2)))
        m3 <- m2
        m4 <- sin(angle)^2/(major^2) + (cos(angle)^2)/(minor^2)
      mygates[[j]][["cov_matrix"]] <- solve(matrix(c(m1, m2, m3, m4), nrow = 2))
      mygates[[j]][["coords"]] <- matrix(unlist(mydef$ellipse$center), ncol = 2)
      colnames(mygates[[j]][["coords"]]) <- mygates[[j]][["channels"]]
      
    } else if(mygates[[j]][["type"]]  == "SplitGate") {
      mydef <- gates$definition[[i]][[1]]
      
      #left
      ind <- mydef[["split"]][["L"]]
      names(mygates)[ind] <- paste0(names(mygates)[j],"_","low") #define gate name
      mygates[[ind]][["channels"]] <- mygates[[j]][["channels"]]
      mygates[[ind]][["coords"]]<- data.frame(-Inf, mydef[["split"]][["x"]])
      colnames(mygates[[ind]][["coords"]]) <- mygates[[ind]][["channels"]]
      mygates[[ind]][["type"]] <- "RangeGate"
      
      
      #right
      ind <- mydef[["split"]][["R"]]
      names(mygates)[ind] <- paste0(names(mygates)[j],"_","high") #define gate name
      mygates[[ind]][["channels"]] <- mygates[[j]][["channels"]]
      mygates[[ind]][["coords"]]<- data.frame(mydef[["split"]][["x"]], Inf)
      colnames(mygates[[ind]][["coords"]]) <- mygates[[ind]][["channels"]]
      mygates[[ind]][["type"]] <- "RangeGate"

    } else if(mygates[[j]][["type"]]  == "RangeGate") {
      mydef <- gates$definition[[i]][[1]]
      mygates[[j]][["coords"]]<- data.frame(mydef$range$x1, mydef$range$x2)
      colnames(mygates[[j]][["coords"]]) <- mygates[[j]][["channels"]]
      
    } else if(mygates[[j]][["type"]]  == "QuadrantGate") {
      mydef <- gates$definition[[i]][[1]]
      
      #UR
      ind <- mydef$quadrant$UR
      mygates[[ind]] <- list()
      mygates[[ind]][["channels"]] <- mygates[[j]][["channels"]]
      mygates[[ind]][["coords"]] <- data.frame(matrix(c(mydef$quadrant$x, Inf, mydef$quadrant$y, Inf), nrow = 2))
      colnames(mygates[[ind]][["coords"]]) <- mygates[[ind]][["channels"]]
      mygates[[ind]][["type"]] <- "QuadrantGate"
      
      #UL
      ind <- mydef$quadrant$UL
      mygates[[ind]] <- list()
      mygates[[ind]][["channels"]] <- mygates[[j]][["channels"]]
      mygates[[ind]][["coords"]] <- data.frame(matrix(c(-Inf, mydef$quadrant$x, mydef$quadrant$y, Inf), nrow = 2))
      colnames(mygates[[ind]][["coords"]]) <- mygates[[ind]][["channels"]]
      mygates[[ind]][["type"]] <- "QuadrantGate"
      
      #LL
      ind <- mydef$quadrant$LL
      mygates[[ind]] <- list()
      mygates[[ind]][["channels"]] <- mygates[[j]][["channels"]]
      mygates[[ind]][["coords"]] <- data.frame(matrix(c(-Inf, mydef$quadrant$x, -Inf, mydef$quadrant$y), nrow = 2))
      colnames(mygates[[ind]][["coords"]]) <- mygates[[ind]][["channels"]]
      mygates[[ind]][["type"]] <- "QuadrantGate"
      
      #LR
      ind <- mydef$quadrant$LR
      mygates[[ind]] <- list()
      mygates[[ind]][["channels"]] <- mygates[[j]][["channels"]]
      mygates[[ind]][["coords"]] <- data.frame(matrix(c(mydef$quadrant$x, Inf, -Inf, mydef$quadrant$y), nrow = 2))
      colnames(mygates[[ind]][["coords"]]) <- mygates[[ind]][["channels"]]
      mygates[[ind]][["type"]] <- "QuadrantGate"
      
    }else  {
      warning(paste(mygates[[j]][["type"]], "are not currently supported!"))
      break
      
    }
    #print(str(mygates[[i]]))
    #mygates[[j]][["coords"]]
    #print(mygates[[j]][["coords"]])
    #colnames(mygates[[j]][["coords"]]) <- mygates[[j]][["channels"]]
  }
  return(mygates)
}


exp_pops
gate_defs <- define_gates(exp_info$exp_gates, exp_info$exp_lut)
my.gated.ff <- gate_population(my.ff,
                exp_pops$name[8],
                exp_info$exp_pops,
                exp_info$exp_gates,
                gate_defs,
                exp_info$exp_lut)
my.gated.ff
exp_pops$name[8]
gate_defs[8]
#flowViz::xyplot(`FSC-A` ~ `SSC-A`, my.ff)
#library(ggplot2)
mydf <- as.data.frame(my.gated.ff@exprs)
#mydf
ggplot(mydf, aes(x = `Pacific Blue-A`, y = `Pacific Orange-A`)) + 
  geom_hex(bins = 128) + 
  scale_x_log10(limits = c(1,250000)) + 
  scale_y_log10(limits = c(1,250000))

gate_defs[7]
?scale_x_log10
?geom_hex()
flowViz::xyplot(`Pacific Orange-A` ~ `Pacific Blue-A`, my.gated.ff)
flowViz:x
summary(my.gated.ff)

colnames(gate_defs$elipsoid$cov_matrix) <- gate_defs$elipsoid$channels
gate.i<- flowCore::ellipsoidGate(gate_defs$elipsoid$cov_matrix, mean = as.numeric(gate_defs$elipsoid$coords))
my.gated.ff
my.gated.ff)
exp_pops$name[5]


exp_po
exp_pops$name[1]
gate_defs

mygates[[4]]




gate_population <- function(flow_frame, population, poplist, gatelist, gate_defs, lut) { #need to pass flow_frame (fcb flowCore file), gatelist (get_gates), population name (string), poplist (get_populations), lut
    # flow_frame <- my.ff
    # population <- exp_pops$name[8]
    # population
    # poplist <- exp_info$exp_pops
    # gatelist <- exp_info$exp_gates
    # lut <- exp_info$exp_lut
  #  gate_defs[[4]]
  # exp_pops
   #gate_defs changes, here is refering to gate_defs object
  
  pop.gates <- unlist(poplist[[match(population, poplist$name), "definition"]][["gates"]]) # gets the sequence of gates (numeric ids) that defines a population
  pop.gates
  # pop.gates
  # exp_gates$gateId
  for(j in pop.gates) { #for gate in pop_gates
    #gate_defs[[j]]
    i<- j
    
    gatelist$gateId
    #print(gatelist$name[[i]]) #can delete
    gate_defs[[i]]
    axis <- as.formula(paste0("`",gate_defs[[i]][["channels"]][2],"`", "~" ,"`",gate_defs[[i]][["channels"]][1],"`"))
    
    channel.ind<- match(gate_defs[[i]][["channels"]],   as.character(lut$shortName))
    channel.char <- gate_defs[[i]][["channels"]]
    
    if (any(lut[channel.ind, "scaleType"] == 4)) { #transforms arcsinh channels approraitely
      ind <- which(lut[channel.ind, "scaleType"] == 4)
      print(ind)
      flowCore::exprs(flow_frame)[,channel.char[ind]] <- asinh(flowCore::exprs(flow_frame)[,channel.char[ind]]/lut[channel.ind[ind], "cofactor"])
    }
    
    print(gate_defs[[i]][["type"]])
    if(gate_defs[[i]][["type"]] == "RectangleGate") {
      print("rectangle gate gated")
      gate.i <- flowCore::rectangleGate(gate_defs[[i]][["coords"]])
      flow_frame <- gatein(flow_frame, gate.i)
    } else if (gate_defs[[i]][["type"]] == "PolygonGate") {
      print("polygon gate gated")
      gate.i <- flowCore::polygonGate(gate_defs[[i]][["coords"]])
      flowViz::xyplot(axis, flow_frame, filter=gate.i)
      flow_frame <- gatein(flow_frame, gate.i)
      
    } else if (gate_defs[[i]][["type"]] == "EllipseGate") {
      print("elipse gate gated")
      colnames(gate_defs[[i]][["cov_matrix"]]) <- gate_defs[[i]][["channels"]]
      gate.i<- flowCore::ellipsoidGate(gate_defs[[i]][["cov_matrix"]], mean = as.numeric(gate_defs[[i]][["coords"]]))
      flow_frame<- gatein(flow_frame, gate.i)
    }  else if (gate_defs[[i]][["type"]] == "RangeGate") {
      print("range gate gated")
      flow_frame <- flow_frame[
        flow_frame@exprs[,gate_defs[[i]][["channels"]][1]] > gate_defs[[i]][["coords"]][[1]] & 
          flow_frame@exprs[,gate_defs[[i]][["channels"]][1]] < gate_defs[[i]][["coords"]][[2]],
        ]
    } else if (gate_defs[[i]]["type"] == "QuadrantGate") {
      #colnames(gate_defs[[i]][["coords"]])[[1]]
      flow_frame<- flow_frame[
        flow_frame@exprs[,gate_defs[[i]][["channels"]][1]] > gate_defs[[i]][["coords"]][[1,1]] & 
          flow_frame@exprs[,gate_defs[[i]][["channels"]][1]] < gate_defs[[i]][["coords"]][[2,1]] &
          flow_frame@exprs[,gate_defs[[i]][["channels"]][2]] > gate_defs[[i]][["coords"]][[1,2]] &
          flow_frame@exprs[,gate_defs[[i]][["channels"]][2]] < gate_defs[[i]][["coords"]][[2,2]]
          ,]
    }
    

    #plotgate(what = axis, flow_frame, gate.i, ylim = c(-2, 12))
    #plotgate(what = axis, flow_frame, gate.i)
    
    if (any(lut[channel.ind, "scaleType"] == 4)) { #backtransforms arcsinh channels back to linear scale
      ind <- which(lut[channel.ind, "scaleType"] == 4)
      flowCore::exprs(flow_frame)[,channel.char[ind]] <- sinh(flowCore::exprs(flow_frame)[,channel.char[ind]])*lut[channel.ind[ind], "cofactor"]
    }
    #print(nrow(exprs(flow_frame)))
  }
  return(flow_frame)
}

