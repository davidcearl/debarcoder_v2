###############################################################################

## ui
debarcode_select_ui <- function(id){
    ns <- NS(id)
    
    tagList(
      fluidRow(
        # h1('My Text'),
        # textOutput(ns('bc_type')),
        tabsetPanel(
          tabPanel("Level-1", uiOutput(ns('debarcoder_ui_1'))),
          tabPanel("Level-2", uiOutput(ns('debarcoder_ui_2')))
        )
       # uiOutput(ns('mode_ui'))
      )
    )
}

## server
#caller expects return of list of flow_set and exp_info
debarcode_select <- function(input, output, session, fcb_dfs, bc_type, x){
 # print(bc_type)
  output$bc_type <- renderText(bc_type())
  debarcode_l1 <- callModule(run_debarcoder, 'run_debarcoder1', fcb_dfs, x,
                             channel = fcb_dfs()[[1]][['bc1_channel']])
  debarcode_l2 <- callModule(run_debarcoder, 'run_debarcoder2', fcb_dfs, x,
                             channel = fcb_dfs()[[1]][['bc2_channel']], bc_prev.df = debarcode_l1)
  
  output$debarcoder_ui_1 <- renderUI({
    ns <- session$ns
    return(run_debarcoder_ui(ns('run_debarcoder1')))
  })
  
  output$debarcoder_ui_2 <- renderUI({
    ns <- session$ns
    return(run_debarcoder_ui(ns('run_debarcoder2')))
  })
  
  return(debarcode_l2)
    # demo_exp <- callModule(demo_mode, 'demo_mode', x)
    # api_exp <- callModule(api_mode, 'api_mode', x)
    # callModule(reproduce_mode, 'reproduce_mode')

    # output$mode_ui <- renderUI({
    #     ns <- session$ns
    #     if(input$mode == 'demo'){
    #         return(demo_mode_ui(ns('demo_mode')))
    #     }
    #     else if(input$mode == 'api'){
    #         return(api_mode_ui(ns('api_mode')))
    #     }
    #     else if(input$mode == 'reproduce'){
    #         return(reproduce_mode_ui(ns('reproduce_mode')))
    #     }
    # })
    # 
    # exp <- reactive({
    #     ns <- session$ns
    #     if(input$mode == 'demo'){
    #         return(demo_exp())
    #     }
    #     else if(input$mode == 'api'){
    #         return(api_exp())
    #     }
    # })

    #return(bc_type)
}
