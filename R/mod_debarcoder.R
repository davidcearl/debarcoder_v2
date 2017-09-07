###############################################################################

## helper functions

##redundant?
debarcode_bc1 <- function(fcb_df, bc1_df, channel, levels, uccutoff, opt,
                          updateProgress = NULL, subsample = 10e3){

    mydf <- debarcode_1(fcb_df = fcb_df, bc_single_level = bc1_df,
                        channel = channel, levels = levels,
                        uccutoff = uccutoff, opt = opt,
                        subsample = subsample, updateProgress = updateProgress)

    return(mydf)
}


#takes fcb_df_po_debarcoded, ??not pb_df??
#requires debarcoder.R
#which parameters should user be able to set in shiny ui?
debarcode_bc2 <- function(bc1_debarcoded, prevchannel, channel, levels){
    mydf2 <- debarcode.2(fcb_df = bc1_debarcoded,
                         prevchannel = prevchannel,
                         channel = channel,
                         levels = levels)
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
            #actionButton(ns('skip'), label = 'Skip BC2'),
            h4('BC2 assignment plot'),
            plotOutput(ns('assign_plot2'))
        )
    )
}

###
# modulue reactive fuction
###

run_debarcoder <- function(input, output, session, fcb_dfs) {

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

        debarcoded_bc1 <- debarcode_bc1(fcb_df = fcb_df,
                                        bc1_df = bc1_df,
                                        channel = bc1_channel,
                                        levels = input$bc1_levels,
                                        uccutoff = input$bc1_unccutoff,
                                        opt = input$bc1_model,
                                        updateProgress= updateProgress,
                                        subsample = input$subsample
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
        ggplot2::ggplot(df, ggplot2::aes_string(y = yaxis, x = xaxis)) +
          ggplot2::scale_y_continuous(trans = logicle_trans()) +
          ggplot2::geom_point(ggplot2::aes(col = as.factor(df$bc1))) +
          ggplot2::scale_color_discrete(name = "BC1 Levels") +
          ggplot2::ggtitle("BC 1 assignments")
    })

    dbc2 <- eventReactive(input$submit_dbc2, {
        ch <- names(fcb_dfs()[[1]])


        channel1 <- fcb_dfs()[[1]][['bc1_channel']]
        channel2 <- fcb_dfs()[[1]][['bc2_channel']]

        debarcoded_bc2 <- debarcode_bc2(bc1_debarcoded = dbc1()[['df']],
                                        prevchannel = channel1,
                                        channel = channel2,
                                        levels = input$bc2_levels)

        modulelog <- list("bc1" = list("coefs"  = dbc1()[["regressionmodel"]],
                                       "snorm" = dbc1()[["snorm"]],
                                       "bc1levels" = input$bc1_levels,
                                       "unccutoff" = input$bc1_unccutoff),
                              "bc2" = list("bc2levels" = "bc2_levelves")
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
    return(dbc2)
}