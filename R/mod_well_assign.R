well_assignment_ui <- function(id) {
    ns <- NS(id)
    tagList(
        selectInput(ns("bc_layout"), label = 'Barcoded Rows or Columns', c("column" = "col", "row" = 'row'), selected = 'column'),
        selectInput(ns("bc_highest"), label = 'Highest Level', c('A' = 1, 'B'= 2, 'C' =3, 'D' =4, 'E'= 5, 'F' =6, 'G' =7, 'H' = 8)),
        selectInput(ns("bc_lowest"), label = 'Lowest Level', c('A' = 1, 'B'= 2, 'C' =3, 'D' =4, 'E'= 5, 'F' =6, 'G' =7, 'H' = 8))
    )
}

###
# modulue reactive fuction
###

well_assignment <- function(input, output, session) {
    
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
