## ui

run_debarcoder_ui <- function(id) {
    ns <- NS(id)

    tagList(
        fluidRow(
            column(6,
            wellPanel(
              radioButtons(ns('morpho_method'), label = "1. Morphology Correction Method", 
                           choices = c("Regression on Scatter" = "regression", "Normalize to Uptake Control" = "controldye"), 
                           selected = "regression"),
              sliderInput(ns("subsample"), label = "Cells to subsample", min = 100, max = 50000, value = 10000, step = 100, ticks = TRUE),
              conditionalPanel(
                  condition = paste0("input['", ns("morpho_method"), "'] == 'controldye'"),
                  selectInput(ns('uptake_channel'), label = "Uptake Control",
                              choices = c("Gathering Channels")),
                  radioButtons(ns('split_on_bc'), label = "Split by Previous Level?", 
                               choices = c("Yes", "No"),
                               selected = "Yes")
              )),
            wellPanel(
              radioButtons(ns('modeling_method'), label = "2. Population Modeling Method", 
                           choices = c("Guassian Mixture Modeling" = "mixture", "Jenks Natural Breaks" = "fisher"), 
                           selected = "mixture"),
              sliderInput(ns('bc1_levels'), label = "Number of BC Levels", min = 1, max = 12, value = 6, step = 1, ticks = FALSE),
              
              conditionalPanel(
                condition = paste0("input['", ns("modeling_method"), "'] == 'mixture'"),
                selectInput(ns('bc1_model'), label = "Model to use",
                            c("Normal" = "Normal", 
                              "Student's T" = "t",
                              "Skew-Normal" = "Skew.normal"), 
                            selected = "Skew.normal") 
              )
            ),
            wellPanel(
              h5(strong("3. Cell Assignment Cutoffs")),
              sliderInput(ns('bc1_likelihood'), label = "Likelihood Cutoff", min = 1, max = 32, value = 8, step = 1),
              conditionalPanel(
                condition = paste0("input['", ns("modeling_method"), "'] == 'mixture'"),
                sliderInput(ns('bc1_unccutoff'), label = "Abiguity Cutoff", min = 0, max = 0.5, value = 0.05, step = 0.01)
              )
            ),
            actionButton(ns('submit_dbc1'), label = 'Run Debarcoder on BC1')
            ),
            column(6,
              h4('BC1 assignment plot'),
              plotOutput(ns("assignment_plot.d"))
            )
        ),
        fluidRow()
    )
}

###
# modulue reactive fuction
###

run_debarcoder <- function(input, output, session, fcb_dfs, x, channel, bc_prev.df = NULL) {
  # morpho_uptake_ui <- eventReactive(input$morpho_method == "controldye"){ 
  #   renderUI(
  #       selectInput(ns('uptake_channel'), label = "Uptake Control",
  #                   choices = c("Gathering Channels")),
  #       radioButtons(ns('split_on_bc'), label = "Split by Previous Level?", 
  #                    choices = c("Yes", "No"),
  #                    selected = "Yes"),
  #       sliderInput(ns("subsample"), label = "Cells to subsample", min = 100, max = 50000, value = 10000, step = 100, ticks = TRUE)
  #     )
  #   )
  # }
  
  observeEvent(input$morpho_method, {
    if(input$morpho_method == "controldye") {
      channels <- colnames(fcb_dfs()[[1]][[1]])
  
      if(is.null(bc_prev.df)) {
        selected <- head(channels, 1)
      } else {
        selected <- fcb_dfs()[[1]][['bc1_channel']]
      }
      updateSelectInput(session, 'uptake_channel', choices = channels, selected = selected) 
    }
  })

   split_on_bc <- eventReactive(input$submit_dbc1, {
     ch <- names(fcb_dfs()[[1]])
     #fcb_file <- ch[1]
     bc1_channel <- channel 
     #print(str(bc_prev.df))
     if(is.null(bc_prev.df)) {
       fcb_df <- fcb_dfs()[[1]][[1]]
     } else {
       fcb_df <- bc_prev.df()[['db']]
     }
     #write.csv(fcb_df, "fcb_df.csv")
     bc1_df <- fcb_dfs()[[1]][[2]]
     if(is.null(bc1_df)){
       bc1_df <- fcb_df
     }
     #print(nrow(fcb_df))
     fcb_df$ind <- 1:nrow(fcb_df)
     ###

     
     
     #split up the new channels by the older channel
     if(input$split_on_bc == "Yes" & input$morpho_method == "controldye") {
       #print("split_on_bc was TRUE")
       fcb_df.list <- split(fcb_df, fcb_df[, "bc1"])
     } else {
       fcb_df.list <- list("1" = fcb_df)
     }
     return(fcb_df.list)
   })
   
   debarcoded_output <- eventReactive(split_on_bc(), {
     
     lut <- fcb_dfs()[["modulelog"]][["exp_info"]][["exp_lut"]]
     cofactor_bc1<- lut[(which(lut$shortName == channel)), "cofactor"]
     
     if(input$morpho_method == "controldye") {
       cofactor_uptake <- lut[(which(lut$shortName == input$uptake_channel)), "cofactor"]
       
     } else {
       cofactor_uptake <- NULL
     } 
     fcb_df.list <- split_on_bc()
     progress <- shiny::Progress$new()
     progress$set(message = "Debarcoding BC1: ", value = 0)
     on.exit(progress$close())
     n <- 5
     updateProgress <- function(detail = NULL) {
       progress$inc(amount = 1/n, detail = detail)
     }
     
      assignments.list <- vector("list", length = length(fcb_df.list))
      #print(str(fcb_df.list))
      for(i in 1:length(fcb_df.list)) {
        #print(length(fcb_df.list))
        #print("Names:")
        #print(names(fcb_df.list))
        #print(names(fcb_df.list)[[i]])
        fcb_df.i <- fcb_df.list[[i]]
        
        if(names(fcb_df.list)[i] == "0") {
          assignments.list[[i]][[1]] <- data.frame("pop" = rep(0, nrow(fcb_df.list[[i]])),
                                                   "ind" = fcb_df.i[,"ind"])
          next()
        }
        
        if(TRUE){ #will need to make this rective when the singly stained controlls are added
          bc1_df.i <- fcb_df.i
        }
        vec <- morphology_corr(fcb_df.i,
                               bc_single_level = bc1_df.i ,
                               channel = channel, 
                               opt = input$morpho_method,
                               subsample = input$subsample,
                               trans = 'arcsinh',
                               updateProgress = updateProgress,
                               cofactor_bc1 = cofactor_bc1,
                               cofactor_uptake = cofactor_uptake,
                               uptake_channel = input$uptake_channel
                               )
        
        #print("line161")
        probs <- fit_models(vec,
                            levels = input$bc1_levels,
                            opt = input$modeling_method, 
                            dist = input$bc1_model,
                            cofactor = cofactor_bc1)
        
        updateProgress <- function(detail = NULL) {
          progress$inc(amount = 1/n, detail = "Assigning Cells")
        }
        
  
        fcb_df.i[,channel] <- vec
        print(str(assignments.list[[i]]))
        assignments.list[[i]]<- assign_cells(fcb_df.i,
                                             probs,
                                             likelihoodcut = input$bc1_likelihood,
                                             ambiguitycut = input$bc1_unccutoff, 
                                             output = "both", 
                                             channel = channel)
        #print(summary(fcb_df.i[,"ind"]))
        assignments.list[[i]][[1]] <- data.frame("pop" = assignments.list[[i]][[1]],
                                                 "ind" = fcb_df.i[,"ind"]) #for later matching
        print(nrow(assignments.list[[i]][[1]] ))
        #print(summary(assignments.list[[i]][[1]][,"ind"]))
        #print(str(assignments.list[[i]][[1]]))
        #ggplot2::ggsave(paste0("level", i, ".png"), assignments.list[[i]][[2]])

      }
      
      return(assignments.list)
    })
    
    output$assignment_plot.h <- renderPlot({
      debarcoded_output()[[1]][["plot"]]
    })
    
    output$assignment_plot.d <- renderPlot({
      lut <- fcb_dfs()[["modulelog"]][["exp_info"]][["exp_lut"]]
      
      plot.df <- assigned_cells()
      
      bc_n <- tail(colnames(plot.df),1)
      names(plot.df) <- pracma::strRep(pracma::strRep(colnames(plot.df), "-", '_'), " ", "_") #sanitize for ggplot
      yaxis <- pracma::strRep(pracma::strRep(channel, "-", '_'), " ", "_")
      ytrans <- asinh_trans(lut[(which(lut$shortName == channel)), "cofactor"])
      if(input$morpho_method == "regression") {
        xaxis <- 'SSC_A'
        xtrans <- "identity"
      } else {
        xaxis <- pracma::strRep(pracma::strRep(input$uptake_channel, "-", '_'), " ", "_")
        xtrans <- asinh_trans(lut[(which(lut$shortName == input$uptake_channel)), "cofactor"])
        
      }
      
      
      plot.df[,bc_n] <- factor(plot.df[,bc_n], levels = c(names(tail(table(plot.df[,bc_n]),-1)), "0"), ordered = TRUE)
      levels(plot.df[,bc_n])[length(levels(plot.df[,bc_n]))] <- "Uncertain"
      #print(table(plot.df$bc1))
      pop.yields <- as.data.frame(table(plot.df[,bc_n])) #obtaining counts
      #print(colnames(plot.df))
      df.split <- split(plot.df, plot.df[,bc_n])
      mydf <- lapply(df.split, #obtain medians for each population
                     function(plot.df) {
                       c(median(plot.df[,xaxis]),
                         median(plot.df[,yaxis]))
                     }
      ) 
      
      
      medians <- as.data.frame(do.call(rbind, mydf))
      pop.yields.df <- cbind(pop.yields, medians)
      colnames(pop.yields.df) <- c(bc_n, "count", xaxis, yaxis)
      
      pop.yields.df$precent <- pop.yields.df$count/sum(pop.yields.df$count)*100 #percent of cells in each pop
      pop.yields.df$label <- paste0(prettyNum(pop.yields.df$count, big.mark=","), #pretty labels
                                    " ", "(", round(pop.yields.df$precent,1), "%)")
      
      #these lines clean-up the uncertain cells population label
      pop.yields.df[nrow(pop.yields.df), "label"] <- paste("Uncertain: \n", tail(pop.yields.df$label,1)) 
      #print(pop.yields.df)
      pop.yields.df[pop.yields.df[,bc_n] == "Uncertain",c(xaxis,yaxis)] <- c(quantile(plot.df[,xaxis], 0.99),quantile(plot.df[,yaxis], 0.01))
      
      #reset the colors such that uncertain is grey
      ggplotColours <- function(n = 6, h = c(0, 360) + 15){
        if ((diff(h) %% 360) < 1) h[2] <- h[2] - 360/n
        hcl(h = (seq(h[1], h[2], length = n)), c = 100, l = 65)
      }

        mycolors <- c(ggplotColours(input$bc1_levels), "#636363")

      # print(table(df$bc1))
      if (is.null(bc_prev.df)) {
        print("bc_prev.df was null")
        myplot <- ggplot2::ggplot(plot.df, ggplot2::aes_string(y = yaxis, x = xaxis)) +
          ggplot2::scale_y_continuous(trans = ytrans) +
          ggplot2::scale_x_continuous(trans = xtrans) + 
          ggplot2::geom_point(ggplot2::aes(col = as.factor(plot.df[,bc_n]))) +
          ggplot2::geom_label(data = pop.yields.df, ggplot2::aes_string(label = "label", colour = bc_n), show.legend = FALSE) + 
          ggplot2::scale_color_manual(name = "BC1 Levels", values=mycolors) + 
          ggplot2::ggtitle("BC 1 assignments")
        
      } else if(input$split_on_bc == "Yes" & input$morpho_method == "controldye") {
        fcb_df.unlist.split <- split(plot.df, plot.df$bc2)
        print(lapply(fcb_df.unlist.split, summary))
        
        
        print(bc_n)
        print(colnames(plot.df))
        myplot <- ggplot2::ggplot(plot.df, ggplot2::aes_string(y = yaxis, x = xaxis)) +
          ggplot2::scale_y_continuous(trans = ytrans) +
          ggplot2::scale_x_continuous(trans = xtrans) + 
          ggplot2::geom_point(ggplot2::aes_string(col = bc_n)) +
          # ggplot2::geom_label(data = pop.yields.df, ggplot2::aes(label = label, colour = bc1), show.legend = FALSE) + 
          ggplot2::scale_color_manual(name = "BC1 Levels", values=mycolors) + 
          ggplot2::ggtitle("BC 2 assignments") + 
          ggplot2::facet_grid(.~bc1, scales = "free_x")
      } else {
        print(bc_n)
        print(colnames(plot.df))
        myplot<- ggplot2::ggplot(plot.df, ggplot2::aes_string(y = yaxis, x = xaxis)) +
          ggplot2::scale_y_continuous(trans = ytrans) +
          ggplot2::scale_x_continuous(trans = xtrans) + 
          ggplot2::geom_point(ggplot2::aes_string(col = bc_n)) +
          ggplot2::scale_color_manual(name = "BC1 Levels", values=mycolors) + 
          ggplot2::ggtitle("BC 2 assignments")
      }
      return(myplot)
      
    })

    assigned_cells <- eventReactive(debarcoded_output(), {
      fcb_df.unlist <-do.call(rbind, split_on_bc())
      fcb_df.unlist <- fcb_df.unlist[order(fcb_df.unlist["ind"]),]
      fcb_df.unlist <- fcb_df.unlist[,-which(colnames(fcb_df.unlist) == "ind")]
      #print(str(fcb_df.unlist))

      assignments.unlist <-do.call(rbind, lapply(debarcoded_output(), '[[', 1))
      assignments.unlist <- assignments.unlist[order(assignments.unlist["ind"]),]
      #print(str(assignments.unlist))
      
      #print(head(assignments.unlist, 12))
      if (is.null(bc_prev.df)) {
        fcb_df.unlist[,"bc1"] <- as.numeric(assignments.unlist[,"pop"])
      } else {
        fcb_df.unlist[,"bc2"] <- as.numeric(assignments.unlist[,"pop"])
      }
      return(fcb_df.unlist)
    })
    
    assigned_cells_list <- eventReactive(assigned_cells(), {
      ch <- names(fcb_dfs()[[1]])
      
      
      channel1 <- fcb_dfs()[[1]][['bc1_channel']]
      channel2 <- fcb_dfs()[[1]][['bc2_channel']]
      lut <- fcb_dfs()[["modulelog"]][["exp_info"]][["exp_lut"]]
      cofactor_bc1<- lut[(which(lut$shortName == channel1)), "cofactor"]
      cofactor_bc2<- lut[(which(lut$shortName == channel2)), "cofactor"]
      print(table(assigned_cells()[,'bc1']))
      modulelog <- list("bc1" = list("bc1levels" = max(assigned_cells()[,'bc1'])),
                        "bc2" = tryCatch(list("bc2levels" = max(assigned_cells()[,'bc2'])),
                                         error = function(a) return(NULL))
      )
      print(colnames(assigned_cells()))
      return(list("db" = assigned_cells(), "modulelog" = modulelog))
    })
    
    
    observeEvent(input$proceed_button, {
      updateNavbarPage(x, "mainNavbarPage", "tab4")
    })
    return(assigned_cells_list)
}
