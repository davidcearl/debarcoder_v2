source("./R/loadfiles.R")

###############################################################################
ui <- navbarPage(title = 'DebarcodeR', id = "mainNavbarPage",
                 tabPanel(title = "1: Select Mode and Import Experiment",
                          id = 'tab1',
                          value = 'tab1',
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
                 tabPanel(title = "2: Select FCB Population",
                          id = 'tab2',
                          value = 'tab2',
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
                 tabPanel(title = "3: Debarcode",
                          id = 'tab3',
                          value = 'tab3',
                          sidebarLayout(
                            sidebarPanel(
                              helpText("1. Specificy the Morphology correction method.
                                       This corrects for the influence of cellular protein content
                                       on barcoding staining. A dedicated uptake control provides the best results, 
                                       but forward and side scatter can also be used to esimate cellular protien content"),
                              helpText("2. Select a modeling method. Gaussian mixture modeling uses parametric modeling to
                                       estimate the probability for each cell originating from each population, while 
                                       Jenks Natural Breaks Optimization non-parametrically splits the data into 
                                       a specified number of groups for optimial seperation. Choose gaussian mixture modeling if
                                       there is substantial overlap between one more populations.Otherwise, guassian mxiture modeling
                                       provides better results and the cost of increased runtime."),
                              helpText("3. Select the thresholds for assigning cells. If 
                                       gaussian mixture modeling was used, the ambiguity cutoff can also be set. 
                                       Lower thresholds are more stringent and will result in lower cell yields."),
                              helpText("4. Click run BC1 debarcoder"),
                              helpText("This will take a few minutes to
                                       run. If the output plots look good
                                       scroll down to BC2 debarcoder.
                                       Otherwise you can adjust the
                                       options and rerun"),
                              verbatimTextOutput('mytext')
                              ),
                            mainPanel(
                              debarcode_select_ui('debarcode_select')
                              #run_debarcoder_ui("run_db")
                            )
                              )
                              ),
                 tabPanel("4: Well Assignments",
                          value = 'tab4',
                          id = 'tab4',
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
                 tabPanel("5: Upload to Cytobank",
                          value = 'tab5',
                          id = 'tab5',
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
                 tabPanel("6: Edit Sample Tags", value = "
                          
                          
                          
                          tab6",
                          sidebarLayout(
                            sidebarPanel(
                              helpText('foo')
                            ),
                            mainPanel(
                              sample_tag_ui('sample_tag')
                            )
                          )
                 )
                 )

server <- function(input, output, session){
  #x=session syntax allows updateNavbarpage to automatically advance to next
  #tab as approriate
  setup <- callModule(setup, 'setup', x = session)
  
  
  fcb <- callModule(fcb_select, 'fcb_select', setup, x = session)
  
  #debarcoded <- callModule(run_debarcoder, 'run_db', fcb, x = session)
  debarcoded <- callModule(debarcode_select, 'debarcode_select',
                           fcb, reactive(input$bctype), x = session)
  
  callModule(assign_split, 'assign_split', setup, fcb, debarcoded, x = session)
  
  upload_id <- callModule(upload, 'upload', setup)
  
  #sample_tag <- callModule(sample_tag, 'sample_tag', setup, upload_id)
  
  session$onSessionEnded(stopApp)
}

shinyApp(ui, server)