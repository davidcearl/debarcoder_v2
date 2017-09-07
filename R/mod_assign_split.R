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
    for (well_number in valid.wells){
        updateProgress(detail = paste0("well_", well_number, ".fcs"))
        pop <- well_assigned_df[which(well_assigned_df$well==well_number), colnames(flowCore::exprs(gated_flowFrame))]

        pop <- as.matrix(pop)
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
        h3("BC1"),
        well_assignment_ui(ns('bc1')),
        h3("BC2"),
        well_assignment_ui(ns('bc2')),
         verbatimTextOutput(ns('bc2_details')),
         actionButton(ns('write_file'), label = 'save csv'),
         actionButton(ns('split_fcs'), label = 'write fcs files'),
         actionButton(ns('writelog'), label = 'save session')
    )
}

assign_split <- function(input, output, session, setup, fcb, debarcoded){
    bc1_assignment <- callModule(well_assignment, 'bc1')
    bc2_assignment <- callModule(well_assignment, 'bc2')

    #display rows to check well id assignments
    output$bc2_details <- renderPrint({
        head(assign_wells(bc1_assignment(), bc2_assignment(), debarcoded()[[1]])[,c('bc1', 'bc2', 'well')], 20)
    })

    #save well assigned FCB dataframe as csv
    observeEvent(input$write_file, {
        write.csv(x = assign_wells(bc1_assignment(), bc2_assignment(), debarcoded()[[1]]), file = 'well_assigned_exp.csv')
    })


    observeEvent(input$split_fcs, {
        well_assigned_df <- assign_wells(bc1_assignment(), bc2_assignment(), debarcoded()[[1]])
        #print(str(fcb()))
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
}



