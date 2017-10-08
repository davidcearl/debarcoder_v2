###############################################################################
## helper functions

# upload_fcs <- function(cyto_session, expID, new_exp_name,
#                        clone = FALSE, keep_fcs = c(-1), updateProgress = NULL) {
# 
#     init_name <- new_exp_name
#     if(clone == TRUE) {
#         if (is.function(updateProgress)) {
#             updateProgress(detail = "cloning experiment")
#         }
#         new.exp <- CytobankAPI::experiments.clone_selective(cyto_session, expID, init_name, fcs_files = keep_fcs)
#     } else {
#         if (is.function(updateProgress)) {
#             updateProgress(detail = paste("creating:", init_name))
#         }
#         purpose <- paste0(CytobankAPI::experiments.show(cyto_session, expID)[["purpose"]][[1]], ", debarcoded")
#         new.exp <- CytobankAPI::experiments.new(cyto_session, init_name, purpose)
#     }
# 
#     uploadlist <- list.files(file.path(expID, "debarcoded"))
# 
#     for(i in 1:length(uploadlist)) {
#         if (is.function(updateProgress)) {
#             updateProgress(detail = uploadlist[[i]])
#         }
#       CytobankAPI::fcs_files.upload(cyto_session, new.exp$id, file.path(expID, "debarcoded", uploadlist[i]))
#     }
# 
#     #neccessary to get around the multiple panel bug
#     if (is.function(updateProgress)) {
#         updateProgress(detail = "Wrapping up...")
#     }
#     #new.exp.clone <- experiments.clone_selective(cyto_session, new.exp$id, new_exp_name, clone_reagents = TRUE)
# 
#     #experiments.trash(cyto_session, new.exp$id)
# 
#     return(new.exp$id)
# }

upload_fcs <- function(cyto_session, expID, new_exp_name = "newexp",
                       opt = "clone", keep_fcs = c(-1), updateProgress = NULL, old_expID = "") {
  init_name <- new_exp_name
  new.exp <- list()
  if(opt == "clone") {
    if (is.function(updateProgress)) {
      updateProgress(detail = "cloning experiment")
    }
    new.exp$id <- CytobankAPI::experiments.clone_full(cyto_session, expID)
    # new.exp$id <- CytobankAPI::experiments.clone_selective(cyto_session, expID, init_name, fcs_files = keep_fcs)
    
    #?CytobankAPI::experiments.clone_selectiv
  } else if(opt == "oldexp") {
    print(old_expID)
    new.exp$id <- old_expID
  } else {
    if (is.function(updateProgress)) {
      updateProgress(detail = paste("creating:", init_name))
    }
    purpose <- paste0(CytobankAPI::experiments.show(cyto_session, expID)[["purpose"]][[1]], ", debarcoded")
    new.exp <- CytobankAPI::experiments.new(cyto_session, init_name, purpose)
  }
  
  uploadlist <- list.files(file.path(expID, "debarcoded"))
  
  #print(new.exp$id)
  for(i in 1:length(uploadlist)) {
    if (is.function(updateProgress)) {
      updateProgress(detail = uploadlist[[i]])
    }
    CytobankAPI::fcs_files.upload(cyto_session, new.exp$id, file.path(expID, "debarcoded", uploadlist[i]))
  }
  
  #neccessary to get around the multiple panel bug
  if (is.function(updateProgress)) {
    updateProgress(detail = "Wrapping up...")
  }
  #new.exp.clone <- experiments.clone_selective(cyto_session, new.exp$id, new_exp_name, clone_reagents = TRUE)
  
  #experiments.trash(cyto_session, new.exp$id)
  
  return(new.exp$id)
}

###############################################################################





upload_ui <- function(id) {
    ns <- NS(id)
    tagList(radioButtons(ns("cloneOpt"),
                         "Upload Option",
                         c("Clone Experiment" = "clone",
                           "Create New Experiment" = "newexp",
                           "Existing Experiment" = "oldexp")),
            conditionalPanel(paste0("input['",
                                    ns("cloneOpt"),
                                    "'] == 'oldexp' "),
                             DT::dataTableOutput(ns("experiment_table_upload"))#, ns("upload")
                             
            ),
            conditionalPanel(paste0("input['",
                                    ns("cloneOpt"),
                                    "'] == 'newexp' "),
                             textInput(ns("new_exp_name"), "Experiment Name")#, ns("upload")
            ),
            actionButton(ns("upload"), label = "Upload to Cytobank"))
}

###############################################################################

upload <- function(input, output, session, setup) {

    newexp <- reactiveValues(info = list())

    observeEvent(input$upload, {

        exp_id <- setup()[['exp_info']][['exp_id']]

        if(input$cloneOpt == "oldexp") {
          opt = 'oldexp'
          row <- input$experiment_table_upload_rows_selected
          old_exp_id <- unlist(setup()[['exp_info']][['exp_id']])
          #print(old_exp_id)
        } else if (input$cloneOpt == "newexp") {
          opt = 'newexp'
        } else {
          opt = 'clone'
        }
        # Create a Progress object
        progress <- shiny::Progress$new()
        progress$set(message = "Uploading: ", value = 0)
        # Close the progress when this reactive exits (even if there's an error)
        on.exit(progress$close())

        n <- length(list.files(file.path(exp_id, "debarcoded"))) + 1
        updateProgress <- function(detail = NULL) {
            progress$inc(amount = 1/n, detail = detail)
        }

        if(setup()[['mode']] == 'api'){
            newexp$info <- upload_fcs(setup()[['cyto_session']],
                                      exp_id,
                                      input$new_exp_name,
                                      opt = opt,
                                      updateProgress = updateProgress, 
                                      old_expID = old_exp_id)
        }

    })

      
    output$experiment_table_upload <- DT::renderDataTable({
      setup()[['exp_table']]
    }, selection = list(mode = 'single', selected = c(1), target = 'row'))
    
    
    return(newexp)
}
