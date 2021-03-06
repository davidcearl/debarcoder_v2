match_sample_tags <- function(cyto_session, new_expID, user_sampletags, extract = FALSE, updateProgress= NULL) {

    sample_tags <- CytobankAPI::sample_tags.download(cyto_session, new_expID)
    #st = sample tags
    st.df <- read.table(sample_tags, header = T, sep = '\t', stringsAsFactors = F)

    #useful when "well_XX" files are present
    well.ind <- grep("well_", st.df$FCS.Filename)
    well.ind <-well.ind[!well.ind %in% grep("NA", st.df$FCS.Filename)]

    st.df[well.ind, "Plate.Column"] <- as.numeric(stringr::str_extract(st.df$FCS.Filename[well.ind], "[[:digit:]]"))
    st.df[well.ind, "Plate.Row"] <- stringr::str_extract(st.df$FCS.Filename[well.ind], "[[:upper:]]")


    cb.wells <- paste0(st.df$Plate.Row, st.df$Plate.Column)
    user.wells <- paste0(user_sampletags$Row, user_sampletags$Column)
    cb.ind <- match(cb.wells, user.wells) #rows from user_sampletags which match wells

    #edits cytobank template with usertags, matching by indecies
    for (column in colnames(user_sampletags)[!colnames(user_sampletags) %in% c("Row", "Column")]) {
        st.df[, column]<- user_sampletags[cb.ind, column]
    }

    st.df[is.na(st.df)] <- "-" #neccessary for cytobank to correctly handle missing values
    colnames(st.df)<- gsub("\\.", " ", colnames(st.df)) #fixed column names


    return(st.df)
}


sample_tag_ui <- function(id){
    ns <- NS(id)
    tagList(
        selectInput(ns('select_platemap'),
                    label = 'platemap to use',
                    choices = c('Multiplexed Activity Metabolomics',
                                'upload custom')),
        uiOutput(ns('platemap')),
        br(),
        actionButton(ns('upload_sample_tags'), 'upload sample tags'),
        tableOutput(ns('contents'))
    )
}

sample_tag <- function(input, output, session, setup, upload_id){

    output$platemap <- renderUI({
        ns <- session$ns
        if(input$select_platemap == 'upload custom'){
            ui_element <- fluidRow(
                              fileInput(ns('upload_platemap'),
                                        label = 'Upload Platemap',
                                        accept = c('text/csv',
                                                   'text/comma-separated-values,text/plain',
                                                   '.csv')),
                              downloadButton(ns('downloadData'),
                                             'Download Template')
                          )
            return(ui_element)
        }
        #else {
        #    ui_element <- fluidRow(
        #                      actionButton('mam', "Multiplexed Activity Metabolomics")
        #                  )
        #}
    })

    output$downloadData <- downloadHandler(
        filename = "platemap_template.csv",
        content = function(file) {
            table <- read.csv("./data/platemap_template.csv")
            write.csv(table, file, row.names = FALSE, na = "")
        }
    )

    sample_tags_reactive <- eventReactive(input$upload_sample_tags, {

        new_expID <- upload_id()

        if(input$select_platemap == 'upload custom'){

            inFile <- input$upload_platemap

            if (is.null(inFile))
                return(NULL)

            userst.df <- read.csv(inFile$datapath, header=TRUE)

            matched.df <- match_sample_tags(setup()[['cyto_session']], new_expID, userst.df)

            write.table(matched.df, file='mysampletags.tsv', quote=FALSE, sep='\t', row.names = FALSE)

            st.upload <- CytobankAPI::sample_tags.upload(setup()[['cyto_session']], new_expID, "mysampletags.tsv")
            return(list(matched.df))
        }
        else {
            progress <- shiny::Progress$new()
            progress$set(message = "Working: ", value = 0)
            # Close the progress when this reactive exits (even if there's an error)
            on.exit(progress$close())

            n <- 3
            updateProgress <- function(detail = NULL) {
                progress$inc(amount = 1/n, detail = detail)
            }
            if (is.function(updateProgress)) {
                updateProgress(detail = "Downloading sample tags")
            }

            inFile <- "platemap_MAM.csv"

            if (is.null(inFile))
                return(NULL)

            userst.df <- read.csv(paste0('./data/', inFile), header=TRUE)
            userst.df$Timepoints <- stringr::str_pad(userst.df$Timepoints,2,pad = 0)

            matched.df <- match_sample_tags(setup()[['cyto_session']], new_expID, userst.df, updateProgress)

            write.table(matched.df, file='mysampletags.tsv', quote=FALSE, sep='\t', row.names = FALSE)
            if (is.function(updateProgress)) {
                updateProgress(detail = "Uploading Sample Tags")
            }
            st.upload <- CytobankAPI::sample_tags.upload(setup()[['cyto_session']], new_expID, "mysampletags.tsv")
            return(list(matched.df))
        }
    })

    output$contents<- renderTable({
        sample_tags_reactive()[[1]]
    })
}





