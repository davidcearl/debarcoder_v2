source("./R/loadfiles.R")


source("./R/launch.R")
source("./R/mod_assign_split.R")
source("./R/mod_debarcoder.R")
source("./R/mod_exp_info_api.R")
source("./R/mod_fcb_select.R")
source("./R/mod_mode.R")
source("./R/mod_sample_tag.R")
source("./R/mod_setup.R")
source("./R/mod_upload.R")
source("./R/mod_well_assign.R")
source("./R/util_debarcoder.R")
source("./R/util_gate_trans.R")
source("./R/util_regression.R")


###############################################################################
ui <- dashboardPage(dashboardHeader(disable = T),
                    dashboardSidebar(disable = T),
                    dashboardBody(box(width = 12, 
                      tabBox(width=12,id="tabBox_next_previous",
                     tabPanel("1: Select Mode and Import Experiment",
                              value = '1: Select Mode and Import Experiment',
                              sidebarLayout(
                                sidebarPanel(
                                  helpText("1. Before using this tool you need
                                           to have gated the population you
                                           want to debarcode on Cytobank and
                                           made a compensation matrix"),
                                  helpText("2. Enter your login and API token
                                           for Cytobank"),
                                  helpText("3. You can generate an API token by
                                           going to your account settings on
                                           Cytobank"),
                                  helpText("4. Don't try to type the token in,
                                           just copy and paste it in the token
                                           input box"),
                                  helpText("5. click submit and wait for the
                                           connection notification to appear"),
                                  helpText("6. Proceed to the next tab"),
                                  h5("Created by Benjamin Reisman and
                                     David Earl, 2017")
                                  ),
                                mainPanel(
                                  setup_ui('setup')
                                )
                                  )),
                     tabPanel("2: Select FCB Population", value = '2: Select FCB Population"',
                              sidebarLayout(
                                sidebarPanel(
                                  helpText("1. Make sure the right files are
                                           selected for the FCB and comps.
                                           Either the barcoded fcs file or the
                                           singly stained compensation controls
                                           can be used as barcoding controls."),
                                  helpText("2. Choose the FCB population and
                                           compensation to use"),
                                  helpText("3. Events with less than selected
                                           signal intensity will be removed"),
                                  helpText("This should be set to the intensity
                                           of your lowest barcode level"),
                                  helpText("4. Click done when you've finished
                                           setting the options and go to the
                                           next tab")
                                  ),
                                mainPanel(
                                  fcb_select_ui('fcb_select')
                                )
                                  )
                                  ),
                     tabPanel("3: Debarcode", value = '3: Debarcode',
                              sidebarLayout(
                                sidebarPanel(
                                  helpText("1. Specifiy the number of levels
                                           of BC1 with the slider"),
                                  helpText("2. Specifiy the uncertainty
                                           cutoff"),
                                  helpText("Events below the cutoff will not
                                           be assinged to a BC level"),
                                  helpText("3. Select the model to use
                                           for debarcoding"),
                                  helpText("4. Click run BC1 debarcoder"),
                                  helpText("This will take a few minutes to
                                           run. If the output plots look good
                                           scroll down to BC2 debarcoder.
                                           Otherwise you can adjust the
                                           options and rerun")
                                  ),
                                mainPanel(
                                  run_debarcoder_ui("run_db")
                                )
                                  )
                                  ),
                     tabPanel("4: Well Assignments", value = "4: Well Assignments",
                              sidebarLayout(
                                sidebarPanel(
                                  helpText("Please specify the correspondence
                                           between the barcode levels and
                                           the well positions")
                                  ),
                                mainPanel(
                                  assign_split_ui('assign_split')
                                )
                                  )
                     ),
                     tabPanel("5: Upload to Cytobank", value = "5: Upload to Cytobank",
                              sidebarLayout(
                                sidebarPanel(
                                  helpText("1. Choose whether to clone the
                                           original experiment and upload FCS
                                           files or whether to create a new
                                           experiment"),
                                  helpText("2. Chose an experiment name"),
                                  helpText("3. Press Upload")
                                  ),
                                mainPanel(
                                  upload_ui('upload')
                                )
                                )
                              ),
                     tabPanel("6: Edit Sample Tags", value = "6: Edit Sample Tags",
                              sidebarLayout(
                                sidebarPanel(
                                  helpText('foo')
                                ),
                                mainPanel(
                                  sample_tag_ui('sample_tag')
                                )
                              )
                     ),
                     tags$script("
                       $('body').mouseover(function() {
                                 list_tabs=[];
                                 $('#tabBox_next_previous li a').each(function(){
                                 list_tabs.push($(this).html())
                                 });
                                 Shiny.onInputChange('List_of_tab', list_tabs);})
                                 "
                      )
                      ),
                     uiOutput("Next_Previous")
                    ))
)


server <- function(input, output, session){

  # library(shiny)
  # library(DT)
  # library(CytobankAPI)
  # library(flowCore)
  # library(sn)
  # library(rlist)
  # library(MASS)
  # library(gridExtra)
  # library(mclust)
  # library(ggplot2)
  # library(viridis)
  # library(classInt)
  # library(pracma)
  # library(stringr)
  # library(mixsmsn)
  # library(shinythemes)
  Previous_Button=tags$div(actionButton("Prev_Tab",HTML('<div class="col-sm-4"><i class="fa fa-angle-double-left fa-2x"></i></div>
                                                                  ')))
  Next_Button=div(actionButton("Next_Tab",HTML('<div class="col-sm-4"><i class="fa fa-angle-double-right fa-2x"></i></div>')))
  
  output$Next_Previous=renderUI({
    tab_list=input$List_of_tab[-length(input$List_of_tab)]
    nb_tab=length(tab_list)
    
    print(tab_list)
    print(input$tabBox_next_previous)
    print(nb_tab)
    if (which(tab_list==input$tabBox_next_previous)==nb_tab)
      column(1,offset=1,Previous_Button)
    else if (which(tab_list==input$tabBox_next_previous)==1)
      column(1,offset = 10,Next_Button)
    else
      div(column(1,offset=1,Previous_Button),column(1,offset=8,Next_Button))
    
  })
  
  observeEvent(input$Prev_Tab,
               {

                 tab_list=input$List_of_tab
                 current_tab=which(tab_list==input$tabBox_next_previous)
                 updateTabsetPanel(session,"tabBox_next_previous",selected=tab_list[current_tab-1])
               }
  )
  observeEvent(input$Next_Tab,
               {
                 tab_list=input$List_of_tab
                 current_tab=which(tab_list==input$tabBox_next_previous)
                 print(current_tab)
                 updateTabsetPanel(session,"tabBox_next_previous",selected=tab_list[current_tab+1])
               }
  )
  source("./R/loadfiles.R")
  
  setup <- callModule(setup, 'setup')

  fcb <- callModule(fcb_select, 'fcb_select', setup)

  debarcoded <- callModule(run_debarcoder, 'run_db', fcb)

  callModule(assign_split, 'assign_split', setup, fcb, debarcoded)

  upload_id <- callModule(upload, 'upload', setup)

  #sample_tag <- callModule(sample_tag, 'sample_tag', setup, upload_id)

  #session$onSessionEnded(stopApp)
  

}

shinyApp(ui, server)

