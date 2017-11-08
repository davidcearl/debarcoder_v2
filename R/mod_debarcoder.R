###############################################################################

## helper functions

##redundant?
debarcode_bc1 <- function(fcb_df, bc1_df, channel, levels, uccutoff, opt,
                          updateProgress = NULL, subsample = 10e3, cofactor =NULL){
    mydf <- debarcode_1(fcb_df = fcb_df, bc_single_level = bc1_df,
                        channel = channel, levels = levels,
                        uccutoff = uccutoff, opt = opt, trans = 'arcsinh',
                        subsample = subsample, updateProgress = updateProgress, cofactor_bc1 = cofactor)

    return(mydf)
}


#takes fcb_df_po_debarcoded, ??not pb_df??
#requires debarcoder.R
#which parameters should user be able to set in shiny ui?
debarcode_bc2 <- function(bc1_debarcoded, prevchannel, channel, levels, cofactor_bc1 = NULL, cofactor_bc2 = NULL){
    mydf2 <- debarcode.2(fcb_df = bc1_debarcoded,
                         prevchannel = prevchannel,
                         channel = channel,
                         levels = levels, 
                         cofactor_bc1 = cofactor_bc1,
                         cofactor_bc2 = cofactor_bc2)
    return(mydf2)
}

## ui


run_debarcoder_ui <- function(id) {
    ns <- NS(id)

    tagList(
        fluidRow(
            #selectInput(ns("bc1_file"), label = 'BC1', c('none')),
            #selectInput(ns("fcb_file"), label = 'FCB', c('none')),
            sliderInput(ns('bc1_levels'), label = "Number of BC Levels", min = 1, max = 12, value = 6, step = 1, ticks = FALSE),
            sliderInput(ns('bc1_unccutoff'), label = "Assignment Cutoff", min = 0, max = 0.5, value = 0.05, step = 0.01),
            selectInput(ns('bc1_model'), label = "Model to use", c("Fit skew.normal" = 4)),
            actionButton(ns('submit_dbc1'), label = 'Run Debarcoder on BC1'),
            sliderInput(ns("subsample"), label = "Subsample?", min = 100, max = 50000, value = 10000, step = 100, ticks = TRUE),
            h4('BC1 hist plot'),
            plotOutput(ns("hist_plot")),
            h4('BC1 assignment plot'),
            plotOutput(ns('assign_plot'))
        ),
        fluidRow(
            #selectInput(ns("bc2_file"), label = 'BC2', c('none')),
            sliderInput(ns('bc2_levels'), label = "Number of BC2 Levels", min = 1, max = 12, value = 8, step = 1, ticks = FALSE),
            helpText("Select `1` if only a single barcoding channel is present"),
            sliderInput(ns('bc2_unccutoff'), label = "Assignment Cutoff", min = 0, max = 0.5, value = 0.05, step = 0.01),
            helpText("Not currently implemented"),
            actionButton(ns('submit_dbc2'), label = 'Run Debarcoder on BC2'),
            actionButton(ns('proceed_button'), label = 'Proceed'),
            #actionButton(ns('skip'), label = 'Skip BC2'),
            h4('BC2 assignment plot'),
            plotOutput(ns('assign_plot2'))
        )
    )
}

###
# modulue reactive fuction
###

run_debarcoder <- function(input, output, session, fcb_dfs, x) {
    dbc1 <- eventReactive(input$submit_dbc1, {

        #####
        ch <- names(fcb_dfs()[[1]])
        #fcb_file <- ch[1]
        bc1_channel <- fcb_dfs()[[1]][['bc1_channel']]
        bc2_channel <- fcb_dfs()[[1]][['bc2_channel']]
        fcb_df <- fcb_dfs()[[1]][[1]]
        bc1_df <- fcb_dfs()[[1]][[2]]
        if(is.null(bc1_df)){
            bc1_df <- fcb_df
        }
        ###

        progress <- shiny::Progress$new()
        progress$set(message = "Debarcoding BC1: ", value = 0)
        # Close the progress when this reactive exits (even if there's an error)
        on.exit(progress$close())

        n <- 5
        updateProgress <- function(detail = NULL) {
            progress$inc(amount = 1/n, detail = detail)
        }

        #print(str(fcb_dfs()))
        ch <- names(fcb_dfs()[[1]])
        bc1_channel <- fcb_dfs()[[1]][['bc1_channel']]

        lut <- fcb_dfs()[["modulelog"]][["exp_info"]][["exp_lut"]]
        cofactor_bc1<- lut[(which(lut$shortName == bc1_channel)), "cofactor"]

        debarcoded_bc1 <- debarcode_bc1(fcb_df = fcb_df,
                                        bc1_df = bc1_df,
                                        channel = bc1_channel,
                                        levels = input$bc1_levels,
                                        uccutoff = input$bc1_unccutoff,
                                        opt = input$bc1_model,
                                        updateProgress= updateProgress,
                                        subsample = input$subsample, 
                                        cofactor = cofactor_bc1
        )
    })

    output$hist_plot <- renderPlot({
        dbc1()[["plot"]]
    })

    output$assign_plot <- renderPlot({
      df <- dbc1()[['df']]
      names(df) <- pracma::strRep(pracma::strRep(names(df), "-", '_'), " ", "_") #sanitize for ggplot
      yaxis <- pracma::strRep(pracma::strRep(dbc1()[['channel']], "-", '_'), " ", "_")
      xaxis <- 'SSC_A'
      
      #print(names(df))
      
      # print(str(df$bc1))
      # print(table(df$bc1))
      
      df$bc1 <- factor(df$bc1, levels = c(names(tail(table(df$bc1),-1)), "0"), ordered = TRUE)
      levels(df$bc1)[length(levels(df$bc1))] <- "Uncertain"
      
      pop.yields <- as.data.frame(table(df$bc1)) #obtaining counts
      
      df.split <- split(df, df$bc1)
      mydf <- lapply(df.split, #obtain medians for each population
                     function(df) {
                       c(median(df[,xaxis]),
                         median(df[,yaxis]))
                     }
      ) 
      
      
      medians <- as.data.frame(do.call(rbind, mydf))
      pop.yields.df <- cbind(pop.yields, medians)
      colnames(pop.yields.df) <- c("bc1", "count", xaxis, yaxis)
      
      pop.yields.df$precent <- pop.yields.df$count/sum(pop.yields.df$count)*100 #percent of cells in each pop
      pop.yields.df$label <- paste0(prettyNum(pop.yields.df$count, big.mark=","), #pretty labels
                                    " ", "(", round(pop.yields.df$precent,1), "%)")
      
      #these lines clean-up the uncertain cells population label
      pop.yields.df[nrow(pop.yields.df), "label"] <- paste("Uncertain: \n", tail(pop.yields.df$label,1)) 
      pop.yields.df["Uncertain",c(xaxis,yaxis)] <- c(quantile(df[,xaxis], 0.99), quantile(df[,yaxis], 0.01))
      
      #reset the colors such that uncertain is grey
      ggplotColours <- function(n = 6, h = c(0, 360) + 15){
        if ((diff(h) %% 360) < 1) h[2] <- h[2] - 360/n
        hcl(h = (seq(h[1], h[2], length = n)), c = 100, l = 65)
      }
      
      mycolors <- c(ggplotColours(length(table(df$bc1))-1), "#636363")
      
      # print(table(df$bc1))
      ggplot2::ggplot(df, ggplot2::aes_string(y = yaxis, x = xaxis)) +
        ggplot2::scale_y_continuous(trans = logicle_trans()) +
        ggplot2::geom_point(ggplot2::aes(col = as.factor(df$bc1))) +
        ggplot2::geom_label(data = pop.yields.df, ggplot2::aes(label = label, colour = bc1), show.legend = FALSE) + 
        ggplot2::scale_color_manual(name = "BC1 Levels", values=mycolors) + 
        ggplot2::ggtitle("BC 1 assignments")
    })

    dbc2 <- eventReactive(input$submit_dbc2, {
        ch <- names(fcb_dfs()[[1]])


        channel1 <- fcb_dfs()[[1]][['bc1_channel']]
        channel2 <- fcb_dfs()[[1]][['bc2_channel']]
        lut <- fcb_dfs()[["modulelog"]][["exp_info"]][["exp_lut"]]
        cofactor_bc1<- lut[(which(lut$shortName == channel1)), "cofactor"]
        cofactor_bc2<- lut[(which(lut$shortName == channel2)), "cofactor"]
        

        debarcoded_bc2 <- debarcode_bc2(bc1_debarcoded = dbc1()[['df']],
                                        prevchannel = channel1,
                                        channel = channel2,
                                        levels = input$bc2_levels, 
                                        cofactor_bc1 = cofactor_bc1,
                                        cofactor_bc2 = cofactor_bc2)

        modulelog <- list("bc1" = list("coefs"  = dbc1()[["regressionmodel"]],
                                       "snorm" = dbc1()[["snorm"]],
                                       "bc1levels" = input$bc1_levels,
                                       "unccutoff" = input$bc1_unccutoff),
                              "bc2" = list("bc2levels" = input$bc2_levels)
        )
        return(list("db" = debarcoded_bc2, "modulelog" = modulelog))
    })

    output$assign_plot2 <- renderPlot({
        df <- dbc2()[[1]]

        ch <- names(fcb_dfs()[[1]])

        channel1 <- fcb_dfs()[[1]][['bc1_channel']]
        channel2 <- fcb_dfs()[[1]][['bc2_channel']]

        names(df) <- pracma::strRep(pracma::strRep(names(df), "-", '_'), " ", "_") #sanitize for ggplot
        yaxis <- pracma::strRep(pracma::strRep(channel2, "-", '_'), " ", "_")
        xaxis <- pracma::strRep(pracma::strRep(channel1, "-", '_'), " ", "_")

        ggplot2::ggplot(df, ggplot2::aes(color = factor(bc2))) +
          ggplot2::scale_y_continuous(trans = logicle_trans()) +
          ggplot2::scale_x_continuous(trans = logicle_trans()) +
          ggplot2::geom_point(ggplot2::aes_string(y = yaxis, x = xaxis), alpha = 0.2, size = 0.5) +
          ggplot2::facet_grid(.~bc1, scales = "free") +
          ggplot2::scale_color_discrete(name = "BC2 Levels")
    })
    
    observeEvent(input$proceed_button, {
      updateNavbarPage(x, "mainNavbarPage", "tab4")
    })
    
    
    return(dbc2)
}
