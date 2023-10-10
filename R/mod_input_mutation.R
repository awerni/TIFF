mutationInputModeUI <- function(id){
  
  ns <- NS(id)
  
  list(
    fluidRow(
      column_3( 
        br(),
        textInput(
          inputId = ns("symbol"), 
          label = "gene symbol", 
          value = ""
        )
      ),
      column_5(
        br(),
        h3(textOutput(ns("genename"))))
    ),
    br(),
    tabsetPanel(
      id = ns("SelectType"),
      tabPanel(
        "Table",
        column_4(),
        column_8(
          br(),
          fluidRow(
            style = "text-align: right;",
            column_6(
              actionButton(ns("select_all"), "select all"),
              actionButton(ns("unselect_all"), "un-select all")
            ),
            column_6(
              fluidRow(
                additionalColumnsUI_sidebar(
                  ns("mutationDT")
                )
              )
            )
          )
        ),
        additionalColumnsUI_main(
          ns("mutationDT")
        )
      ),
      tabPanel(
        "Plot",
        restoreSelectionInputModeUI(
          id = ns("display")
        )
      )
    )
  )
}

mutationInputMode <- function(input, output, session, species, TissuePrefilter){
  
  GeneMutation <- reactive({
    symbol <- input$symbol
    if (is.null(symbol)) return()
    if (nchar(symbol) < 3) return()
    getGeneFromSymbol(symbol, species)
  })
  

  TissueMutation <- reactive({
    gene <- GeneMutation()$ensg
    prefilter <- TissuePrefilter()
    req(gene, prefilter)
    
    withProgress(
      expr = getInfoMutation(gene, prefilter),
      message = "Retrieving DNA status data",
      value = 0.3
    )
  })
  
  tableInfo <- additionalColumns(
    id = "mutationDT",
    Table = TissueMutation,
    defaultCols = c("tissuename", "tumortype", "aamutated"),
    rownames = FALSE,
    escape = FALSE,
    selection = list(
      mode = "multiple",
      target = "row"
    ),
    filter = "top",
    options = list(
      pageLength = 20,
      lengthMenu = c(10, 15, 20, 100)
    )
  )
  
  output$genename <- renderText(GeneMutation()$name)
  
  proxy <- tableInfo$proxy
  
  observeEvent(
    eventExpr = input$select_all, 
    handlerExpr = {
      proxy %>% DT::selectRows(tableInfo$AllRows())
    }
  )
  
  observeEvent(
    eventExpr = input$unselect_all, 
    handlerExpr = {
      proxy %>% DT::selectRows(list())
    }
  )
  
  
  Items <- callModule(
    module = restoreSelectionInputMode,
    id = "display",
    classStack = StashedData
  )
  
  StashedData <- reactive({
    TissueMutation()[tableInfo$AllRows(),]
  })
  
  observeEvent(
    eventExpr = input[["display-column"]], 
    handlerExpr = {
      updateSelectInput(
        session = session,
        inputId = "display-column",
        selected = "aamutated"
      )
    }, 
    once = TRUE, 
    ignoreInit = TRUE
  )
  
  Mutation_ti <- reactive({
    
    if(input$SelectType == "Table") {
      selRow <- tableInfo$SelectedRows()
      validate(need(selRow, "no row selected"))
      proxy %>% DT::selectRows(list())
      
      tissuenames <- TissueMutation() %>% 
        slice(selRow) %>% 
        pull(tissuename)
      
      list(
        tissuename = tissuenames, 
        source = NA
      )
    } else {
      Items()
    }
  })
  
  Mutation_ti
}
