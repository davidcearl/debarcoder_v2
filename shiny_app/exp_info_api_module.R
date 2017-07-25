
require(CytobankAPI)

###############################################################################

## helper functions

exp_id_from_downloaded_files <- function() {
    exp_id <- unlist(list.dirs())[unlist(lapply(list.dirs(), 
                                                grepl, 
                                                pattern = '^./[0-9]+$'))]
    exp_id <- substr(exp_id, 3, nchar(exp_id))
    return(exp_id)
}

#get exp population list from cytobank api
get_populations <- function(cyto_session, exp_id){
    return(populations.list(cyto_session, exp_id,  output = "default"))
}

#get exp gate list from cytobank api
get_gates <- function(cyto_session, exp_id){
    return(gates.list(cyto_session, exp_id, output = "default"))
}

get_compensations <- function(cyto_session, exp_id) {
    return(compensations.list(cyto_session, exp_id, output = 'default'))
}

#name lookup table for consistent naming?
# "Panel 1" is hardcoded
#needs exp_id from shiny ui
get_lut <- function(cyto_session, exp_id ) {
    scales <- scales.list(cyto_session, 
                          exp_id,  
                          output = "default")
    mypanel <- panels.list(cyto_session, 
                           exp_id, 
                           output = "default")[["Panel 1"]][["channels"]]
    scales_df <- as.data.frame(lapply(scales, function(X) unname(unlist(X))))
    mypanel_df <- as.data.frame(lapply(mypanel, function(X) unname(unlist(X))))
    lut <- merge.data.frame(mypanel_df, 
                            scales_df, 
                            by.x= "normalizedShortNameId", 
                            by.y ="normalizedShortNameId")
    return(lut)
}

## ui

exp_info_ui <- function(id) {
    ns <- NS(id)
    tagList(
        actionButton(ns('get_info_button'), 'Get exp info')
    )
}

## server

exp_info <- function(input, output, session, cyto_session) {
    #get experiment info from cytobank api
    #mayber add check_connection conditional
    exp_info <- eventReactive(input$get_info_button, {
        exp_id <- exp_id_from_downloaded_files()
        exp_pops <- get_populations(cyto_session(), exp_id)
        exp_comps <- get_compensations(cyto_session(), exp_id)
        exp_gates <- get_gates(cyto_session(), exp_id)
        exp_lut <- get_lut(cyto_session(), exp_id)
        return(list('exp_pops' = exp_pops,
                    'exp_comps' = exp_comps,
                    'exp_gates' = exp_gates,
                    'exp_lut' = exp_lut))
    })
    return(exp_info)
}
        
