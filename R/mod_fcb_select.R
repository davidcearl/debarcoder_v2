## helper functions
compensate_files <- function(comp_matrix, flowframes) {
    flowCore_uncomps <- c()
    flowCore_comps <- c()
    f_names <- c()
    for (frame in names(flowframes)){
        flowCore_uncomps <- c(flowCore_uncomps, flowframes[[frame]])
        flowCore_comps <- c(flowCore_comps, flowCore::compensate(flowframes[[frame]], as.matrix(comp_matrix)))
        flowframes[[frame]] <- flowCore::compensate(flowframes[[frame]], as.matrix(comp_matrix))
        f_names <- c(f_names, simplify_file_name(flowframes[[frame]]@description$FILENAME))
    }
    names(flowCore_comps) <- f_names
    names(flowCore_uncomps) <- f_names
    return(list(flowCore_comps, orig = flowCore_uncomps))
}

# define_gates <- function(gates, lut) {
#     mygates <- vector("list", length = length(unlist(gates$id)))
#     names(mygates) <- unlist(gates$name)
# 
#     for(i in 1:length(mygates)) {
#         mygates[[i]][["channels"]] <- c(gates$xNormalizedShortNameId[[i]], gates$yNormalizedShortNameId[[i]])
#         mygates[[i]][["channels"]] <- as.character(lut[match(mygates[[i]][["channels"]],
#                                                              lut$normalizedShortNameId),"shortName"])
# 
#         mygates[[i]][["type"]] <- gates$type[[i]]
# 
#         #print(str(mygates))
#         if(mygates[[i]][["type"]] == "PolygonGate") {
#             mygates[[i]][["coords"]] <- do.call(rbind,
#                                                 lapply(
#                                                     gates$definition[[i]][[1]][["polygon"]][["vertices"]],
#                                                     as.numeric)
#             )
#         }
#         if(mygates[[i]][["type"]] == "RectangleGate") {
#             mygates[[i]][["coords"]] <- matrix(
#                 unlist(gates[["definition"]][[i]][[1]][["rectangle"]]), ncol = 2, byrow = TRUE)
#         }
# 
#         mygates[[i]][["coords"]]
#         colnames(mygates[[i]][["coords"]]) <- mygates[[i]][["channels"]]
#     }
#     return(mygates)
# }

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
    
    mygates[[j]][["channels"]] <- c(gates$xNormalizedShortNameId[[i]], gates$yNormalizedShortNameId[[i]])
    mygates[[j]][["channels"]] <- as.character(lut[match(mygates[[j]][["channels"]],
                                                         lut$normalizedShortNameId),"shortName"])
    
    mygates[[j]][["type"]] <- gates$type[[i]]
    
    if(mygates[[j]][["type"]] == "PolygonGate") {
      mygates[[j]][["coords"]] <- do.call(rbind,
                                          lapply(
                                            gates$definition[[i]][[1]][["polygon"]][["vertices"]],
                                            as.numeric))
      colnames(mygates[[j]][["coords"]]) <- mygates[[j]][["channels"]]
    } else if(mygates[[j]][["type"]] == "RectangleGate") {
      mygates[[j]][["coords"]] <- matrix(
        unlist(gates[["definition"]][[i]][[1]][["rectangle"]]), ncol = 2, byrow = TRUE)
      colnames(mygates[[j]][["coords"]]) <- mygates[[j]][["channels"]]
      
    } else if(mygates[[j]][["type"]] == "EllipseGate") {
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
      break #all gate types are supported at this point, future implementations will also include tailored gates. 
      
    }
    #print(str(mygates[[i]]))
    #mygates[[j]][["coords"]]
    #print(mygates[[j]][["coords"]])
    #colnames(mygates[[j]][["coords"]]) <- mygates[[j]][["channels"]]
  }
  return(mygates)
}

######

#requires gates_transforms.R
#takes ungated flowrame, population, population list, gate list, gate definitions, and lookup table
#returns gate_defs flowframe for poupulation
#requires gates_transforms.R
#takes ungated flowrame, population, population list, gate list, gate definitions, and lookup table
#returns gate_defs flowframe for poupulation
gate_population <- function(flow_frame, population, poplist, gatelist, gate_defs, lut) { #need to pass flow_frame (fcb flowCore file), gatelist (get_gates), population name (string), poplist (get_populations), lut
  
  pop.gates <- unlist(poplist[[match(population, poplist$name), "definition"]][["gates"]]) # gets the sequence of gates (numeric ids) that defines a population
  for(j in pop.gates) { #for gate in pop_gates
    i<- j
    
    axis <- as.formula(paste0("`",gate_defs[[i]][["channels"]][2],"`", "~" ,"`",gate_defs[[i]][["channels"]][1],"`"))
    
    channel.ind<- match(gate_defs[[i]][["channels"]],   as.character(lut$shortName))
    channel.char <- gate_defs[[i]][["channels"]]
    
    if (any(lut[channel.ind, "scaleType"] == 4)) { #transforms arcsinh channels approraitely
      ind <- which(lut[channel.ind, "scaleType"] == 4)
      flowCore::exprs(flow_frame)[,channel.char[ind]] <- asinh(flowCore::exprs(flow_frame)[,channel.char[ind]]/lut[channel.ind[ind], "cofactor"])
    }
    
    if(gate_defs[[i]][["type"]] == "RectangleGate") {
      gate.i <- flowCore::rectangleGate(gate_defs[[i]][["coords"]])
      flow_frame <- gatein(flow_frame, gate.i)
    } else if (gate_defs[[i]][["type"]] == "PolygonGate") {
      gate.i <- flowCore::polygonGate(gate_defs[[i]][["coords"]])
      flow_frame <- gatein(flow_frame, gate.i)
      
    } else if (gate_defs[[i]][["type"]] == "EllipseGate") {
      colnames(gate_defs[[i]][["cov_matrix"]]) <- gate_defs[[i]][["channels"]]
      gate.i<- flowCore::ellipsoidGate(gate_defs[[i]][["cov_matrix"]], mean = as.numeric(gate_defs[[i]][["coords"]]))
      flow_frame<- gatein(flow_frame, gate.i)
    }  else if (gate_defs[[i]][["type"]] == "RangeGate") {
      flow_frame <- flow_frame[
        flow_frame@exprs[,gate_defs[[i]][["channels"]][1]] > gate_defs[[i]][["coords"]][[1]] & 
          flow_frame@exprs[,gate_defs[[i]][["channels"]][1]] < gate_defs[[i]][["coords"]][[2]],
        ]
    } else if (gate_defs[[i]]["type"] == "QuadrantGate") {
      flow_frame<- flow_frame[
        flow_frame@exprs[,gate_defs[[i]][["channels"]][1]] > gate_defs[[i]][["coords"]][[1,1]] & 
          flow_frame@exprs[,gate_defs[[i]][["channels"]][1]] < gate_defs[[i]][["coords"]][[2,1]] &
          flow_frame@exprs[,gate_defs[[i]][["channels"]][2]] > gate_defs[[i]][["coords"]][[1,2]] &
          flow_frame@exprs[,gate_defs[[i]][["channels"]][2]] < gate_defs[[i]][["coords"]][[2,2]]
        ,]
    }
    
    if (any(lut[channel.ind, "scaleType"] == 4)) { #backtransforms arcsinh channels back to linear scale
      ind <- which(lut[channel.ind, "scaleType"] == 4)
      flowCore::exprs(flow_frame)[,channel.char[ind]] <- sinh(flowCore::exprs(flow_frame)[,channel.char[ind]])*lut[channel.ind[ind], "cofactor"]
    }
  }
  return(flow_frame)
}
#get dataframe from flowFrame
#removes events that fall below the bc_cutoff
#returns cleaned dataframe
get_cleaned_df <- function(pop_gated_ff, bc_channels, bc1_cutoff = 1e2, bc2_cutoff = 1e2){

    #print(bc_channels)

    myff <- pop_gated_ff
    mydf <- as.data.frame(flowCore::exprs(myff))

    file_name <- simplify_file_name(myff@description$FILENAME)

    mydf <- mydf[which(mydf[,c(bc_channels[1])] > bc1_cutoff),]
    mydf <- mydf[which(mydf[,c(bc_channels[2])] > bc2_cutoff),]


    return(list("df" = mydf, "file_name" = file_name))
}


## ui

fcb_select_ui <- function(id){
    ns <- NS(id)
    tagList(
        fluidRow(
            h3("Select FCB population and BC comps"),
            column(width = 4,
                   selectInput(ns("db_fcb_fcs"),
                               label = 'FCB File',
                               c('none'))
            ),
            column(width = 3,
                   selectInput(ns("db_pop"),
                               label = 'FCB Population',
                               c('none'))
            ),
            column(width = 4,
                   selectInput(ns("db_fcb_comp"),
                               label = 'Comp Matrix',
                               c('none'))
            )
        ),
        fluidRow(
            column(width = 4,
                   selectInput(ns("db_bc1_channel"),
                               label = 'BC1 Channel',
                               c('none')),
                   selectInput(ns("db_bc1_fcs"),
                               label = 'BC1 Comp (optional)',
                               c('none'))
            ),
            column(width = 4,
                   selectInput(ns("db_bc2_channel"),
                               label = 'BC2 Channel',
                               c('none')),
                   selectInput(ns("db_bc2_fcs"),
                               label = 'BC2 Comp (optional)',
                               c('none'))
            ),
            column(width = 3,
                   selectInput(ns("db_bc_pop"),
                               label = 'Comp Population',
                               c('none'))
            )
        ),
        sliderInput(ns('bc1_sig_cutoff'), label = 'BC1 Signal Cutoff', min = -500, max = 2000, value = 100, step = 50),
        sliderInput(ns('bc2_sig_cutoff'), label = 'BC2 Signal Cutoff', min = -500, max = 2000, value = 100, step = 50),

        actionButton(ns('df_button'), label = 'Done'),
        plotOutput(ns('cutoff_plot'))
    )
}


fcb_select <- function(input, output, session, setup){

    fcs_flowframes <- reactive(
        setup()[['fcs_flowframes']]
    )

    exp_info <- reactive(
        setup()[['exp_info']]
    )

    observe({
        fcs_files <- names(fcs_flowframes())
        updateSelectInput(session, 'db_fcb_fcs', choices = fcs_files, selected = head(fcs_files, 1))
        updateSelectInput(session, 'db_bc1_fcs', choices = c(fcs_files, "none"), selected = "none")
        updateSelectInput(session, 'db_bc2_fcs', choices = c(fcs_files, "none"), selected = "none")

        selected_fcb <- fcs_flowframes()[[input$db_fcb_fcs]]

        if(!is.null(selected_fcb)) {
            channels <- unname(colnames(flowCore::exprs(selected_fcb)))
            updateSelectInput(session, 'db_bc1_channel',
                              choices = channels,
                              selected = head(channels, 2))
            updateSelectInput(session, 'db_bc2_channel',
                              choices = channels,
                              selected = head(channels, 1))
        }

        pop_ch <- exp_info()[['exp_pops']]['name']
        updateSelectInput(session, 'db_pop', choices = pop_ch)
        updateSelectInput(session, 'db_bc_pop', choices = pop_ch)

        comps <- names(exp_info()[['exp_comps']])
        updateSelectInput(session, 'db_fcb_comp', choices = c("internal compensation", comps))
    })

    #returns the matrix of the selected comp
    #maybe should move, might break if input$df_button is clicked before fcs files are downloaded
    #called by fcb_dfs
    comp_to_use <- eventReactive(input$db_fcb_comp, {
        comp_name <- input$db_fcb_comp
        if(comp_name == 'none'){
            comp_name <- names(exp_info()[['exp_comps']])[1]
        } else if (comp_name == "internal compensation") { #pulls internal comp
          print(names(flowCore::description(fcs_flowframes()[[1]])))
          comp.matrix <- flowCore::description(fcs_flowframes()[[1]])[["SPILL"]]
        } else {
          comp.matrix <- exp_info()[['exp_comps']][[comp_name]][['compensation_matrix']]
        }
        return(comp.matrix)
    })

    fcb_dfs <- eventReactive(input$df_button, {

        comp <- comp_to_use()
        comp.output <- compensate_files(comp, fcs_flowframes())
        comped_flowCore_files <- comp.output[[1]]
        uncomp <- comp.output[[2]]
        gate_defs <- define_gates(exp_info()[['exp_gates']],
                                  exp_info()[['exp_lut']])

        print(comped_flowCore_files)

        pop_gated_flowFrame <- gate_population(comped_flowCore_files[[input$db_fcb_fcs]],
                                               input$db_pop,
                                               exp_info()[['exp_pops']],
                                               exp_info()[['exp_gates']],
                                               gate_defs,
                                               exp_info()[['exp_lut']])
        pop_df <- get_cleaned_df(pop_gated_flowFrame,
                                 bc_channels = c(input$db_bc1_channel,
                                                 input$db_bc2_channel),
                                 bc1_cutoff = input$bc1_sig_cutoff,
                                 bc2_cutoff  = input$bc2_sig_cutoff)
        uncompsubset <- uncomp[[input$db_fcb_fcs]]
        if(input$db_bc1_fcs != 'none') {
            bc1_gated_flowFrame <- gate_population(comped_flowCore_files[[input$db_bc1_fcs]],
                                                   input$db_bc_pop,
                                                   exp_info()[['exp_pops']],
                                                   exp_info()[['exp_gates']],
                                                   gate_defs,
                                                   exp_info()[['exp_lut']])
            bc1_df <- get_cleaned_df(bc1_gated_flowFrame,
                                     bc_channels = c(input$db_bc1_channel,
                                                     input$db_bc2_channel),
                                     bc1_cutoff = input$bc1_sig_cutoff,
                                     bc2_cutoff  = input$bc2_sig_cutoff)
        }
        else {
            bc1_df <- list()
            bc1_df[['file_name']] <- 'none'
            bc1_df[['df']] <- NULL
        }

        if(input$db_bc2_fcs != "none") {
            bc2_gated_flowFrame <- gate_population(comped_flowCore_files[[input$db_bc2_fcs]],
                                                   input$db_bc_pop,
                                                   exp_info()[['exp_pops']],
                                                   exp_info()[['exp_gates']],
                                                   gate_defs,
                                                   exp_info()[['exp_lut']])
        bc2_df <- get_cleaned_df(bc2_gated_flowFrame,
                                     bc_channels = c(input$db_bc1_channel,
                                                     input$db_bc2_channel),
                                     bc1_cutoff = input$bc1_sig_cutoff,
                                     bc2_cutoff  = input$bc2_sig_cutoff)
        }
        else  {
            bc2_df <- list()
            bc2_df[['file_name']] <- 'none'
            bc2_df[['df']] <- NULL
        }

        dfs_for_debarcoding <- list(pop_df[['df']],
                                    bc1_df[['df']],
                                    bc2_df[['df']],
                                    pop_gated_flowFrame,
                                    comp, uncompsubset,
                                    input$db_bc1_channel,
                                    input$db_bc2_channel)
        names(dfs_for_debarcoding) <- c(pop_df[['file_name']],
                                        bc1_df[['file_name']],
                                        bc2_df[['file_name']],
                                        "pop_ff",
                                        "comp",
                                        "uncomp",
                                        "bc1_channel",
                                        "bc2_channel")

        fcslog <- (
            list("popfcs" = flowCore::description(pop_gated_flowFrame),
                 "bc1fcs" = tryCatch({flowCore::description(pop_gated_flowFrame)}, error = NULL),
                 "bc2fcs" = tryCatch({flowCore::description(pop_gated_flowFrame)}, error = NULL))
        )
        modulelog <- (
            list("comp" = comp,
                 "gate_defs" = gate_defs,
                 "exp_info" =  exp_info(),
                 "db_pop" = input$db_pop,
                 "db_bc_pop" = input$db_bc_pop,
                 "db_fcb_fcs" = input$db_fcb_fcs,
                 "db_bc1_fcs" = input$db_bc1_fcs,
                 "db_bc2_fcs" = input$db_bc2_fcs,
                 "db_bc1_channel" = input$db_bc1_channel,
                 "db_bc2_channel" = input$db_bc2_channel,
                 "bc1_cutoff" =  input$bc1_sig_cutoff,
                 "bc2_cuttoff" =input$bc2_sig_cutoff
            ))

        return(list("dfs" = dfs_for_debarcoding,
                    "modulelog"= modulelog,
                    "fcslog" = fcslog))
    })

    #tie to button to avoid Error in strRep: Argument 's' must be a character vector.
    output$cutoff_plot <- renderPlot({
        df <- fcb_dfs()[['dfs']][[input$db_fcb_fcs]]

        names(df) <- pracma::strRep(pracma::strRep(names(df), "-", '_'), " ", "_") #sanitize for ggplot
        yaxis <- pracma::strRep(pracma::strRep(input$db_bc2_channel, "-", '_'), " ", "_")
        xaxis <- pracma::strRep(pracma::strRep(input$db_bc1_channel, "-", '_'), " ", "_")

        bc1cut <- input$bc1_sig_cutoff
        bc2cut <- input$bc2_sig_cutoff


        ggplot2::ggplot(df, ggplot2::aes_string(x=xaxis , y = yaxis)) +
          ggplot2::geom_hex(bins = 200) +
          ggplot2::scale_y_continuous(trans= logicle_trans(), name = yaxis,
                               breaks = major.ticks,
                               labels = major.ticks,
                               minor_breaks = minor.ticks) +
          ggplot2::scale_x_continuous(trans= logicle_trans(), name = xaxis,
                               breaks = major.ticks,
                               labels = major.ticks,
                               minor_breaks = minor.ticks) +
          viridis::scale_fill_viridis(option = "magma", begin = 0, end = 1) +
          ggplot2::coord_cartesian(xlim = c(-1e2, 2e5), ylim = c(-1e2, 2e5)) +
          ggplot2::geom_vline(xintercept = bc1cut, col = "orange", size = 1) +
          ggplot2::geom_hline(yintercept = bc2cut, col = "blue", size = 1) +
          ggplot2::ggtitle("Removing non-barcoded cells")
    }, width = 600, height = 600)


    return(fcb_dfs)
}


