###############################################################################
## api_mode module

## helper functions

download_fcs_files <- function(cyto_session, expID, fcsFileIDs,
                               lut = NULL, updateProgress = NULL){

    #removes existing fcs files before downloading newly selected files
    exp_dir <- unlist(list.dirs())[unlist(lapply(list.dirs(),
                                                 grepl,
                                                 pattern = '^./[0-9]+$'))]
    unlink(exp_dir, recursive = TRUE)
    unlist(list.dirs())[unlist(lapply(list.dirs(),
                                      grepl,
                                      pattern = '^./[0-9]+$'))]
    #save initial dir
    base_dir <- getwd()
    getwd()
    #create exp dir and debarcoded dir
    dir.create(file.path(base_dir, expID))
    dir.create(file.path(base_dir, expID, "debarcoded"))

    #create directories for each fcs file
    downloadlist <- list()
    for(fcsFileID in fcsFileIDs){
        dir.create(file.path(base_dir, expID, fcsFileID))
        setwd(file.path(base_dir, expID, fcsFileID))
        if(is.function(updateProgress)){
            updateProgress(detail = lut[lut$id== fcsFileID, 'filename'])
        }

        download_object <- CytobankAPI::fcs_files.download(cyto_session, expID, fcsFileID)
        downloadlist<- append(downloadlist, download_object)
        download_object
    }
    setwd(base_dir)
    return(downloadlist)
}

get_downloaded_exp_fcs_ids <- function(){
    exp_id <- unlist(list.dirs())[unlist(lapply(list.dirs(),
                                                grepl,
                                                pattern = '^./[0-9]+$'))]
    exp_id <- substr(exp_id, 3, nchar(exp_id))
    fcs_ids <- unlist(list.dirs())[unlist(lapply(list.dirs(),
                                                 grepl,
                                                 pattern = '^./[0-9]+/[0-9]+$'))]
    fcs_ids <- sapply(fcs_ids,
                      function(x){
                          substr(x, start = (nchar(exp_id)+4), stop = nchar(x))
                      },
                      USE.NAMES = FALSE)
    return(list(exp = exp_id, fcs = fcs_ids))
}

simplify_file_name <- function(file_name){

    if(grepl("Comp|comp", file_name)){
        file_name <- unlist(strsplit(file_name, "_"))[2]
        split_file_name <- unlist(strsplit(file_name, " "))
        file_name <- split_file_name[1:(length(split_file_name)-2)]
        file_name <- paste(file_name, collapse = " ")
    }
    else{
        file_name <- unlist(strsplit(file_name, '/'))
        file_name <- file_name[length(file_name)]
        file_name <- substr(file_name, 1, (nchar(file_name)-4)) #trim .fcs
    }
    return(file_name)
}

#loads downloaded files to named c() of flowFrames
downloaded_fcs_files_to_flowSet <- function(){
    #get_downloaded_exp_fcs_ids()[[1]]

    if(length(get_downloaded_exp_fcs_ids()[[1]] > 0)){

        flowCore_files <- c()
        f_names <- c()
        #is this matching pattern is too broad?
        for(fcs_file in list.files(pattern = '*.fcs$', recursive = T)){
            flowCore_file <- flowCore::read.FCS(file.path('.', fcs_file))
            flowCore_files <- c(flowCore_files, flowCore_file)
            f_names <- c(f_names,
                         simplify_file_name(flowCore_file@description$FILENAME))
        }

        names(flowCore_files) <- f_names
        return(flowCore_files)
    }
    else{
        return('no fcs files')
    }
}

## ui function

api_mode_ui <- function(id){
    ns <- NS(id)
    tagList(
        selectInput(ns("domain"),
                    label = 'Domain',
                    choices = c('vanderbilt',
                                'cellmass'),
                    selected = 'vanderbilt'),
        passwordInput(ns("token"),
                      label = 'API Token',
                      placeholder = "Copy + paste from Cytobank"),
        actionButton(ns("submit_connection"),
                     label = "Submit"),
        textOutput(ns("connected")),
        titlePanel("User's Experiments"),
        DT::dataTableOutput(ns("experiment_table")),
        verbatimTextOutput(ns("selected_row")),
        titlePanel("Selected Experiment FCS Files"),
        DT::dataTableOutput(ns("fcs_table")),
        verbatimTextOutput(ns("filepaths")),
        actionButton(ns("submit_download"), label = 'Download FCS files'),
        exp_info_ui(ns('exp_info'))
    )
}

## server function

api_mode <- function(input, output, session){
    login_info <- reactive({
        #will only run authenticate if user name and token are entered
        validate(c(
            need(input$token, message = "stil waiting for token")))
        token <- input$token
        return(token)
    })
    #funky behavior now, i think when renderUI is called
    #input$submit_connection resets
    #cyto_session must be triggered again
    cyto_session <- eventReactive(input$submit_connection, {
        userlogin <- login_info()
        attempt <- CytobankAPI::authenticate(site=input$domain,
                                auth_token = userlogin,
                                timeout = 360)
        return(attempt)
    })

    #not sure how this works
    check_if_connected <- function(cyto_session) {
        connection_status <- tryCatch({
            exp.list <- CytobankAPI::experiments.list(cyto_session, timeout = 1)
        },
        error = function(cond) {
            if(grepl("Timeout", cond)) {
                return(list('bool' = TRUE,
                            'status' = "connected"))
            }
            else {
                return(list('bool' = FALSE,
                            'status' = cond))
            }
        },
        warning = function(cond) {
            return(list('bool' = FALSE,
                        'status' = cond))
        })
        return(connection_status)
    }

    check_if_connected_reactive <- reactive({
        check_if_connected(cyto_session())
    })

    output$connected <- renderText({
        connection_status <- check_if_connected_reactive()
        if(connection_status[['bool']]) {
            return(paste("connected to", cyto_session()@site))
        }
        else {
            return(connection_status[['status']][[1]])
        }
    })

    exp_table <- reactive({
        connection_status <- check_if_connected_reactive()[['bool']]
        if(connection_status){
            withProgress(message = 'fetching experiments from cytobank',
                         detail = 'This may take a while...',
                         value = 0,
                         expr = {
                             exps <- tryCatch({
                               CytobankAPI::experiments.list(cyto_session())[, c('id', 'experimentName')]
                             })
                             setProgress(value = 1,
                                         message = 'Done!')
                             return(exps)
                         })
        }
        else{
            return(data.frame())
        }
    })


    output$experiment_table <- DT::renderDataTable({
        exp_table()
    }, selection = list(mode = 'single',
                        selected = c(1),
                        target = 'row'))

    fcs_table <- eventReactive(input$experiment_table_rows_selected, {
        row <- input$experiment_table_rows_selected
        exp_id <- exp_table()[row,1]
        fcs_table <- CytobankAPI::fcs_files.list(cyto_session(), exp_id)
        fcs_table[,c('id', 'filename')]
    })

    output$fcs_table <- DT::renderDataTable({
        fcs_table()
    }, selection = list(mode = 'multiple',
                        selected = c(1),
                        target = 'row'))

    output$test_fcs <- renderPrint({
        rows <- input$fcs_table_rows_selected
        fcsIDs <- unlist(fcs_table()[rows,1])
        fcsIDs
    })

    #downloads fcs files
    output$filepaths <- eventReactive(input$submit_download, {
        exp_row <- input$experiment_table_rows_selected
        fcs_rows <- input$fcs_table_rows_selected

        exp_ind <- unlist(exp_table()[exp_row,1])
        fcs_ind <- unlist(fcs_table()[fcs_rows,1])

        n <- 3
        #Create a Progress object
        progress <- shiny::Progress$new()
        progress$set(message = "Downloading: ", value = 0)
        # Close the progress when this reactive exits (even if there's an error)
        on.exit(progress$close())

        updateProgress <- function(detail = NULL) {
            progress$inc(amount = 1/n, detail = detail)
        }

        download_fcs_files(cyto_session = cyto_session(),
                           expID = exp_ind,
                           fcsFileIDs = fcs_ind,
                           lut = fcs_table(),
                           updateProgress = updateProgress)

    })

    #checks for fcs files in wd
    #if there are it reads them into flowframes
    #otherwise fcs_flowframes has character value ""
    fcs_flowframes <- reactivePoll(1000,
                                   session,
                                   checkFunc = function(){
                                       if(length(get_downloaded_exp_fcs_ids()[[1]] > 0)) {
                                           file.info(get_downloaded_exp_fcs_ids()[[1]])$mtime[1]
                                       }
                                       else{
                                           ""
                                       }
                                   },
                                   valueFunc = function(){
                                       obj <- downloaded_fcs_files_to_flowSet()
                                       #print(names(obj))
                                       obj
                                       #print("converting to flowset")
                                   }
    )

    exp_info <- callModule(exp_info, 'exp_info', cyto_session)

    api_exp <- reactive({
        return(list('mode' = 'api',
                    'fcs_flowframes' = fcs_flowframes(),
                    'exp_info' = exp_info(),
                    'cyto_session' = cyto_session()))
    })

    return(api_exp)


}

###############################################################################
## demo_mode module


demo_mode_ui <- function(id){
    ns <- NS(id)
    tagList(
        actionButton(ns('submit_demo'), 'Run demo')
    )
}

demo_mode <- function(input, output, session){
    demo_exp <- eventReactive(input$submit_demo, {
        load('./data/demo_data.Rdata')
        return(demo_data)
    })

    return(demo_exp)
}


###############################################################################
## reproduce_mode module

reproduce_mode_ui <- function(id){
    ns <- NS(id)
    tagList(
        selectInput(ns('select_domain'),
                    label = 'select_domain',
                    choices = c('silver',
                                'cellmass'),
                    selected = 'silver'),
        verbatimTextOutput(ns('selected_domain'))
    )
}

reproduce_mode <- function(input, output, session){
    output$selected_domain <- renderText({
        input$select_domain
    })
}

