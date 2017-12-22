well_assignment_ui <- function(id) {
    ns <- NS(id)
    tagList(
        textInput(ns("well1"), label = 'Top Left Well', value = "A1", placeholder = "Format: A1"),
        selectInput(ns("nrow"), label = 'Number of barcoded rows',  c(1:8)),
        selectInput(ns("bc_row"), label = 'Which barcode corresponds to rows?', c("BC1", "BC2")),
        selectInput(ns("ncol"), label = 'Number of barcoded columns',  c(1:12)),
        selectInput(ns("bc_col"), label = 'Which barcode corresponds to columns?', c("BC1", "BC2"))
    )
}

###
# modulue reactive fuction
###

well_assignment <- function(input, output, session) {

    origin_well <- reactive({
      validate(
        origin_well <- grep("^[[:upper:]]{1}[[:digit:]]+$", input$well1, value = TRUE),  
        need(length(origin_well) != 0,
             "Please enter the starting well as a capital letter, followed by a number, ie: A1")
      )
      origin_well
    })
  
    observe({
        x <-  input$bc_layout
        
        row_ids <- c('A' = 1, 'B'= 2, 'C' = 3, 'D' = 4, 'E'= 5, 'F' = 6, 'G' = 7, 'H' = 8)
        col_ids <- c("1" = 1, "2" = 2, "3" = 3, "4" = 4, "5" = 5, "6" = 6, "7" = 7, "8" = 8, "9" = 9, "10" = 10, "11" = 11, "12" = 12)
        
        if (x == 'row') {
            x <- row_ids
        }
        else { 
            x <- col_ids
        }
        
        updateSelectInput(session, "bc_highest",
                          choices = x,
                          selected = head(x, 1)
        )
        updateSelectInput(session, "bc_lowest",
                          choices = x,
                          selected = tail(x, 1)
        )
    })
    
    bc_info <- reactive({
        info <- list("layout" = input$bc_layout, "high" = as.integer(input$bc_highest), "low" = as.integer(input$bc_lowest))
        info
    })
    return(bc_info)
}
