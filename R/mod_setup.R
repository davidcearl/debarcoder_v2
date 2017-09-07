###############################################################################

## ui
setup_ui <- function(id){
    ns <- NS(id)
    tagList(
        selectInput(ns('mode'),
                    label = 'Select mode',
                    choices = c('demo',
                                'api',
                                'reproduce'),
                    selected = 'demo'),
        uiOutput(ns('mode_ui'))

    )
}

## server
#caller expects return of list of flow_set and exp_info
setup <- function(input, output, session){
    demo_exp <- callModule(demo_mode, 'demo_mode')
    api_exp <- callModule(api_mode, 'api_mode')
    callModule(reproduce_mode, 'reproduce_mode')

    output$mode_ui <- renderUI({
        ns <- session$ns
        if(input$mode == 'demo'){
            return(demo_mode_ui(ns('demo_mode')))
        }
        else if(input$mode == 'api'){
            return(api_mode_ui(ns('api_mode')))
        }
        else if(input$mode == 'reproduce'){
            return(reproduce_mode_ui(ns('reproduce_mode')))
        }
    })

    exp <- reactive({
        ns <- session$ns
        if(input$mode == 'demo'){
            return(demo_exp())
        }
        else if(input$mode == 'api'){
            return(api_exp())
        }
    })

    return(exp)
}
