library(shiny)
source('./shiny_app/setup_module.R')
source('./shiny_app/fcb_select_module.R')
source('./shiny_app/debarcoder_module.R')
source('./shiny_app/assign_split_module.R')
source('./shiny_app/upload_module.R')

###############################################################################
ui <- navbarPage(title = 'DebarcodeR',
                 tabPanel("1: Select Mode and Import Experiment", value = 'tab1',
                          sidebarLayout(
                              sidebarPanel(
                                  helpText("1. Before using this tool you need to have gated the population you want to debarcode on Cytobank and made a compensation matrix"),
                                  helpText("2. Enter your login and API token for Cytobank"),
                                  helpText("3. You can generate an API token by going to your account settings on Cytobank"),
                                  helpText("4. Don't try to type the token in, just copy and paste it in the token input box"),
                                  helpText("5. click submit and wait for the connection notification to appear"),
                                  helpText("6. Proceed to the next tab"),
                                  h5("Created by Benjamin Reisman and David Earl, 2017")
                              ),
                              mainPanel(
                                  setup_ui('setup')
                              )
                          )),
                 tabPanel("2: Select FCB Population", value = 'tab2',
                          sidebarLayout(
                              sidebarPanel(
                                  helpText("1. Make sure the right files are selected for the FCB and comps. Either the barcoded fcs file or the singly stained compensation controls can be used as barcoding controls."),
                                  helpText("2. Choose the FCB population and compensation to use"),
                                  helpText("3. Events with less than selected signal intensity will be removed"),
                                  helpText("This should be set to the intensity of your lowest barcode level"),
                                  helpText("4. Click done when you've finished setting the options and go to the next tab")
                              ),
                              mainPanel(
                                  fcb_select_ui('fcb_select')
                              )
                          )
                 ),
                 tabPanel("3: Debarcode", value = 'tab3',
                          sidebarLayout(
                              sidebarPanel(
                                  helpText("1. Specifiy the number of levels of BC1 with the slider"),
                                  helpText("2. Specifiy the uncertainty cutoff"),
                                  helpText("Events below the cutoff will not be assinged to a BC level"),
                                  helpText("3. Select the model to use for debarcoding"),
                                  helpText("4. Click run BC1 debarcoder"),
                                  helpText("This will take a few minutes to run.  If the output plots look good scroll down to BC2 debarcoder. Otherwise you can adjust the options and rerun")
                              ),
                              mainPanel(
                                  run_debarcoder_ui("run_db")
                              )
                          )
                 ),
                 tabPanel("4: Well Assignments", value = 'tab4',
                          sidebarLayout(
                              sidebarPanel(
                                  helpText("Please specify the correspondence between the barcode levels and the well positions")
                              ),
                              mainPanel(
                                      assign_split_ui('assign_split') 
                              )
                          )
                 ),
                 tabPanel("5: Upload to Cytobank", value = 'tab5',
                          sidebarLayout(
                              sidebarPanel(
                                  helpText("1. Choose whether to clone the original experiment and upload FCS files or whether to create a new experiment"),
                                  helpText("2. Chose an experiment name"),
                                  helpText("3. Press Upload")
                              ),
                              mainPanel(
                                  upload_ui('upload')
                              )
                          )
                 )
)

server <- function(input, output, session){
    
    setup <- callModule(setup, 'setup')
     
    fcb <- callModule(fcb_select, 'fcb_select', setup)
     
    debarcoded <- callModule(run_debarcoder, 'run_db', fcb)
    
    assign_split <- callModule(assign_split, 'assign_split', setup, fcb, debarcoded)
    
    upload <- callModule(upload, 'upload', setup)
    
    #session$onSessionEnded(stopApp)
}

shinyApp(ui, server)
