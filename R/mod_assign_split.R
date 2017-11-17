assign_wells <- function(bc1, bc2, debarcoded_df) {

    debarcoded_df <- debarcoded_df[!debarcoded_df$bc2 == 0,]

    colbc <- bc1
    rowbc <- bc2

    if(bc1[["layout"]] == 'row'){
        colbc <- bc2
        rowbc <- bc1
    }

    col_equ <- function(col_bc_level, highest_col_bc) {return((highest_col_bc + (col_bc_level-1)))}
    #col high > low
    if(colbc[["high"]] > colbc[["low"]]){
        col_equ <- function(col_bc_level, highest_col_bc) {return((highest_col_bc - (col_bc_level - 1)))}
        print(1)
    }

    row_equ <- function(row_bc_level, highest_row_bc) {return(12*((highest_row_bc + (row_bc_level-1))-1))}
    #row high > low
    if(rowbc[["high"]] > rowbc[["low"]]){
        row_equ <- function(row_bc_level, highest_row_bc) {return(12*((highest_row_bc - (row_bc_level-1))-1))}
        print(3)
    }

    colvals <- col_equ(debarcoded_df$bc1, colbc[[2]])
    rowvals <- row_equ(debarcoded_df$bc2, rowbc[[2]])

    if(bc1[["layout"]] == 'row'){
        colvals <- col_equ(debarcoded_df$bc2, colbc[[2]])
        rowvals <- row_equ(debarcoded_df$bc1, rowbc[[2]])
    }



    debarcoded_df$well <- paste0(rep(LETTERS[1:8], each = 12)[rowvals+colvals], stringr::str_pad(colvals, 2, pad  = "0"))

    return(debarcoded_df)
}

clean_ff <- function(flowframe) {
    #flowframe <- A1
    short.names <- character()
    long.desc <- character()
    fcsdesc <- flowframe@description
    if(TRUE) {
        short.names <- grep("\\$P\\d+N", names(fcsdesc), value = TRUE) #eg: "$P1N"
        long.desc <- (fcsdesc[short.names]) #list of short names eg: "FSC-A"
        name.numbers <- stringr::str_extract_all(short.names, "\\d")
        name.numbers.char <- unlist(lapply(name.numbers, function(vec) paste(vec, sep="", collapse="")))
        names(long.desc) <- paste0("$P", name.numbers.char, "S") #generates longname.names "$P1S"

        #long.desc<- lapply(long.desc, function(x) x <- "<NA>") #changes all short names to "<NA>"
        newdesc <- append(fcsdesc, long.desc) #adds longnames to the description parameter list
        flowCore::description(flowframe) <-  newdesc#updates the description

        flowframe@parameters@data$desc <- unlist(long.desc)[order(as.numeric(name.numbers.char))]
        flowframe
        return(flowframe) #return the flowframe
    } else {
        return(flowframe)
    }
}

split_and_write_fcs <- function(well_assigned_df, gated_flowFrame, comp_matrix, eID,
                                updateProgress = NULL, origin_flowFrame = NULL){
    dir.create(paste0(eID, "/debarcoded/"), recursive = TRUE)

    fcscolnames <- colnames(well_assigned_df)[1:(ncol(well_assigned_df)-3)]
    valid.wells <-unique(well_assigned_df$well)[!grepl("NA", unique(well_assigned_df$well))]
    valid.wells<- valid.wells[unlist(lapply(lapply(valid.wells, length), function(x) x > 0))] #removes empty wells. 
    #print(str(valid.wells))
    for (well_number in valid.wells){
        updateProgress(detail = paste0("well_", well_number, ".fcs"))
        pop <- well_assigned_df[which(well_assigned_df$well==well_number), colnames(flowCore::exprs(gated_flowFrame))]

        pop <- as.matrix(pop)
        #print(str(pop))
        #print(attributes(pop))
        colnames(pop) <- fcscolnames
        pop_ff <- flowCore::flowFrame(pop,
                                      description = flowCore::description(gated_flowFrame),
                                      parameters = flowCore::parameters(gated_flowFrame))
        pop_ff.uncomp <- flowCore::compensate(pop_ff, solve(comp_matrix))
        flowCore::description(pop_ff.uncomp) <- flowCore::description(gated_flowFrame)
        flowCore::parameters(pop_ff.uncomp) <- flowCore::parameters(gated_flowFrame)
        pop_ff.uncomp <- clean_ff(pop_ff.uncomp)

        suppressWarnings(flowCore::write.FCS(pop_ff.uncomp, paste0(eID,"/", "/debarcoded/well_",well_number,".fcs")))
    }
}



assign_split_ui <- function(id){
    ns <- NS(id)
    tagList(
      textInput(ns("well1"), label = 'Top Left Well', value = "A1", placeholder = "Format: A1"),
      selectInput(ns("nrow"), label = 'Number of barcoded rows',  c(1:8)),
      selectInput(ns("ncol"), label = 'Number of barcoded columns',  c(1:12)),
      selectInput(ns("bc_row"), label = 'Which barcode corresponds to rows?', c("BC1", "BC2")),
      actionButton(ns('check_layout'), label = "Verify Layout"),
      splitLayout(cellWidths = c("50%", "50%"), plotOutput(ns("BC1plot")), plotOutput(ns("BC2plot"))),
      #verbatimTextOutput(ns('bc2_details')),
      actionButton(ns('write_file'), label = 'save csv'),
      actionButton(ns('split_fcs'), label = 'write fcs files'),
      actionButton(ns('writelog'), label = 'save session'),
      br(), 
      br(), 
      actionButton(ns('proceed_button'), label = 'Proceed')
    )
}

assign_split <- function(input, output, session, setup, fcb, debarcoded, x){
  origin_well <- eventReactive(input$well1, {
    origin <- grep("^[[:upper:]]{1}[[:digit:]]+$", input$well1, value = TRUE) 
    validate(
      need(length(origin_well) != 0,
           "Please enter the starting well as a capital letter, followed by a number, ie: A1")
    )
    row1<- as.numeric(which(LETTERS == gsub("[[:digit:]]+", "", origin)))
    col1<- as.numeric(gsub("[[:upper:]]{1}", "", origin))
    return(c(row1 = row1, col1 = col1))
  })
  
  bc1_row <- eventReactive(input$bc_row, {
    if(input$bc_row == "BC1"){
      return(TRUE)
    } else {
      return(FALSE)
    }
  })
  
layouts <- eventReactive(input$check_layout, {
      row1 <- origin_well()[["row1"]]
      col1 <- origin_well()[["col1"]]
      nrow <- as.numeric(input$nrow)
      ncol <- as.numeric(input$ncol)
      bc1_row <- bc1_row()
      bc1_levels <- debarcoded()[["modulelog"]][["bc1"]][["bc1levels"]]
      bc2_levels <- debarcoded()[["modulelog"]][["bc2"]][["bc2levels"]]
      layout_bc_well.m <- matrix(paste0(
                              rep(LETTERS[row1:(row1+nrow-1)], each = ncol),
                              stringr::str_pad(rep(col1:(col1+ncol-1), times = nrow), 2, pad  = "0")
                              ),
                       ncol = bc2_levels,
                       byrow = bc1_row)
      
      
      layoutbc1 <- matrix(rep(1:bc1_levels, each = bc2_levels),
                          nrow = nrow, ncol = ncol, byrow = bc1_row)
      layoutbc2 <- matrix(rep(1:bc2_levels, times = bc1_levels),
                          nrow = nrow, ncol = ncol, byrow = bc1_row)
      layout_well_bc.df <- data.frame(row = rep(row1:(row1+nrow-1), times = ncol),
                              col = rep(col1:(col1+ncol-1), each = nrow),
                              colchar = rep(LETTERS[col1:(col1+ncol-1)], each = nrow),
                              bc1 = as.numeric(layoutbc1),
                              bc2 = as.numeric(layoutbc2)
      )
      return(list(layout_bc_well.m = layout_bc_well.m, #barcodes => wells
                  layout_well_bc.df = layout_well_bc.df #Wells => barcodes
                  )
      )
  })
  
  
  
  output$BC1plot <- renderPlot({
    layout_df <- layouts()[["layout_well_bc.df"]]
    row1 <- origin_well()["row1"]
    col1 <- origin_well()["col1"]
    nrow <- as.numeric(input$nrow)
    ncol <- as.numeric(input$ncol)
    bc1_row <- bc1_row()
    bc1_levels <- debarcoded()[["modulelog"]][["bc1"]][["bc1levels"]]
    bc2_levels <- debarcoded()[["modulelog"]][["bc2"]][["bc2levels"]]
    # print(paste("layout_df", str(layout_df)))
    # print(paste("row1", row1))
    # print(paste("col1", col1))
    # print(paste("nrow", nrow))
    # print(paste("ncol", ncol))
    # print(paste("bc1_row", bc1_row))
    # print(paste("bc1_levels", bc1_levels))
    # print(paste("bc2_levels", bc2_levels))
    
    
    myplot <- ggplot2::ggplot(layout_df, ggplot2::aes(x = col, y = row)) + 
      ggplot2::geom_tile(ggplot2::aes(fill = as.factor(-bc1)), col = "grey50") + 
      ggplot2::scale_x_continuous(breaks = col1:(col1+ncol-1), position = "top", expand = c(0,0)) + 
      ggplot2::scale_y_reverse(breaks = row1:(row1+nrow-1), labels = LETTERS[row1:(row1+nrow-1)], expand = c(0,0)) + 
      ggplot2::scale_fill_brewer(type = "seq", palette = "Oranges", guide = FALSE) + 
      ggplot2::geom_text(ggplot2::aes(label = bc1)) + 
      ggplot2::theme_classic()
    return(myplot)
  })
  
  output$BC2plot <- renderPlot({
    layout_df <- layouts()[["layout_well_bc.df"]]
    row1 <- origin_well()["row1"]
    col1 <- origin_well()["col1"]
    nrow <- as.numeric(input$nrow)
    ncol <- as.numeric(input$ncol)
    bc1_row <- bc1_row()
    bc1_levels <- debarcoded()[["modulelog"]][["bc1"]][["bc1levels"]]
    bc2_levels <- debarcoded()[["modulelog"]][["bc2"]][["bc2levels"]]
    print(paste("layout_df", str(layout_df)))
    print(paste("row1", row1))
    print(paste("col1", col1))
    print(paste("nrow", nrow))
    print(paste("ncol", ncol))
    print(paste("bc1_row", bc1_row))
    print(paste("bc1_levels", bc1_levels))
    print(paste("bc2_levels", bc2_levels))
    
    #grDevices::colorRampPalette(RColorBrewer::brewer.pal(9, "Blues"))(4)
    mypalette <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(9, "Blues"))(bc2_levels)
    myplot <- ggplot2::ggplot(layout_df, ggplot2::aes(x = col, y = row)) + 
      ggplot2::geom_tile(ggplot2::aes(fill = as.factor(-bc2)), col = "grey50") + 
      ggplot2::scale_x_continuous(breaks = col1:(col1+ncol-1), position = "top", expand = c(0,0)) + 
      ggplot2::scale_y_reverse(breaks = row1:(row1+nrow-1), labels = LETTERS[row1:(row1+nrow-1)], expand = c(0,0)) + 
      ggplot2::scale_fill_manual(values = mypalette, guide = FALSE) + 
      ggplot2::geom_text(ggplot2::aes(label = bc2)) + 
      ggplot2::theme_classic()
    ?ggplot2::scale_fill_manual()
    return(myplot)
  })
    

    #display rows to check well id assignments
    output$bc2_details <- renderPrint({
        head(assign_wells(bc1_assignment(), bc2_assignment(), debarcoded()[[1]])[,c('bc1', 'bc2', 'well')], 20)
    })

    #save well assigned FCB dataframe as csv
    observeEvent(input$write_file, {
        write.csv(x = assign_wells(bc1_assignment(), bc2_assignment(), debarcoded()[[1]]), file = 'well_assigned_exp.csv')
    })


    observeEvent(input$split_fcs, {
        layout.m <- layouts()[["layout_bc_well.m"]]
        well_assigned_df <- debarcoded()[[1]]
        print(layout.m)
        print(head(well_assigned_df$bc1))
        print(head(well_assigned_df$bc2))
        print(unique(well_assigned_df$bc1))
        print(unique(well_assigned_df$bc2))
        
        well_assigned_df$well <- mapply(function(x, y) layout.m[x, y], well_assigned_df$bc1, well_assigned_df$bc2)
        well_assigned_df[well_assigned_df$bc2 == 0, "well"] <- "NA"
        flow_frame <- fcb()[[1]][['pop_ff']]
        orig.flow.frame <- fcb()[[1]][['uncomp']]
        comp <- fcb()[[1]][['comp']]
        exp_id <- setup()[['exp_info']][['exp_id']]
        progress <- shiny::Progress$new()
        progress$set(message = "Writing FCS File: ", value = 0)
        # Close the progress when this reactive exits (even if there's an error)
        on.exit(progress$close())

        n <- length(unique(well_assigned_df$well))
        updateProgress <- function(detail = NULL) {
            progress$inc(amount = 1/n, detail = detail)
        }


        split_and_write_fcs(well_assigned_df, orig.flow.frame, comp, exp_id, updateProgress = updateProgress)

    })

    observeEvent(input$writelog, {
        sessionlog <- list("import" = fcb()[[3]],
                           "selectpop"= fcb()[[2]],
                           "debarcode" = debarcoded()[[2]],
                           "barcodeinfo" = list("bc1" = bc1_assignment(), "bc2" = bc2_assignment()))
        rlist::list.save(sessionlog, "sessionlog.rds")
    })
    
    observeEvent(input$proceed_button, {
      updateNavbarPage(x, "mainNavbarPage", "tab5")
    })
}



