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
        h3(textOutput(ns("genename")))
      ),
      column_4(
        br(), br(),  
        actionButton(
          inputId = ns("select_all"), 
          label = "select all"
        ),
        actionButton(
          inputId = ns("unselect_all"), 
          label = "un-select all"
        )
      )
    ),
    fluidRow_12(containerDT(ns("table")))
  )
}

mutationInputMode <- function(input, output, session, species, TissuePrefilter){
  GeneMutation <- reactive({
    symbol <- input$symbol
    if (is.null(symbol)) return()
    if (nchar(symbol) < 3) return()
    getGeneFromSymbol(symbol, species)
  })
  
  TableData <- reactive({
    if (is.null(GeneMutation())) return()
    
    withProgress(
      expr = getInfoMutation(GeneMutation()$ensg, TissuePrefilter()),
      message = 'Retrieving DNA status data', 
      value = 0.3
    )
  })
  
  output$genename <- renderText(GeneMutation()$name)
  
  output$table <- DT::renderDataTable(
    expr = {
      df <- TableData()
      validate(need(is.data.frame(df) && nrow(df) > 0, "no sequencing data available..."))
      df <- df %>% 
        mutate(
          dnamutation = gsub(";", "; ", dnamutation), 
          aamutation = gsub(";", "; ", aamutation),
          aamutated = as.factor(aamutated)
        )
      
      if ("tumortype" %in% names(df)){
        df <- df %>% mutate(tumortype = as.factor(tumortype))
      }
      
      df
    },
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
  
  proxy <- DT::dataTableProxy("table")
  
  observeEvent(
    eventExpr = input$select_all, 
    handlerExpr = {
      proxy %>% DT::selectRows(input$table_rows_all)
    }
  )
  
  observeEvent(
    eventExpr = input$unselect_all, 
    handlerExpr = {
      proxy %>% DT::selectRows(list())
    })
  
  Mutation_ti <- reactive({
    selRow <- input$table_rows_selected
    validate(need(selRow, "no row selected"))
    
    proxy %>% DT::selectRows(list())
    
    tissuenames <- TableData() %>% 
      slice(selRow) %>% 
      pull(tissuename)
    
    list(
      tissuename = tissuenames, 
      source = NA
    )
  })
  
  Mutation_ti
}
