###############################################################################

## ui
debarcode_select_ui <- function(id){
    ns <- NS(id)
    
    tagList(
      fluidRow(
        # h1('My Text'),
        # textOutput(ns('bc_type')),
        uiOutput(ns('debarcoder_ui'))
       # uiOutput(ns('mode_ui'))
      )
    )
}

## server
#caller expects return of list of flow_set and exp_info
debarcode_select <- function(input, output, session, fcb_dfs, bc_type, x){
 # print(bc_type)
  output$bc_type <- renderText(bc_type())
  
  debarcode2 <- callModule(run_debarcoder2, 'run_debarcoder2', fcb_dfs, x)
  debarcode3 <- callModule(run_debarcoder3, 'run_debarcoder3', fcb_dfs, x)
  
  output$debarcoder_ui <- renderUI({
    ns <- session$ns
    if(bc_type() == '2dye'){
      #print('29, 2dye')
      return(run_debarcoder2_ui(ns('run_debarcoder2')))
    } else if(bc_type() == '3dye') {
      return(run_debarcoder3_ui(ns('run_debarcoder3')))
    }
  })
  
  mode <- reactive({
    ns <- session$ns
    #print(bc_type())
    if(bc_type() == '2dye'){
      return(debarcode2())
    } else if (bc_type() == '3dye') {
      return(debarcode3())
    }
  })
  
  return(mode)
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
