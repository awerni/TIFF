annotationInputUI <- function(id){
  ns <- NS(id)
  tableId <- ns("table")
  
  list(
    fluidRow(
      column_2(
        br(),
        actionButton(
          inputId = ns("select_all"), 
          label = "select all"
        )
      ),
      column_2(
        br(),
        actionButton(
          inputId = ns("unselect_all"), 
          label = "un-select all"
        )
      ),
      column_3(),
      column_5(
        br(), 
        additionalColumnsUI_sidebar(tableId)
      )
    ),
    br(),
    additionalColumnsUI_main(tableId)
  )
}

annotationInput <- function(input, output, session, Table){
  tableInfo <- additionalColumns(
    id = "table",
    Table = Table,
    defaultCols = NULL,
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
  
  proxy <- tableInfo$proxy
  observeEvent(
    eventExpr = input$unselect_all, 
    handlerExpr = {
      proxy %>% DT::selectRows(list())
    }
  )
  
  observeEvent(
    eventExpr = input$select_all, 
    handlerExpr = {
      proxy %>% DT::selectRows(tableInfo$AllRows())
    }
  )
  
  Tissues <- reactive({
    selRow <- tableInfo$SelectedRows()
    proxy %>% DT::selectRows(list())
    
    tissuenames <- Table()[selRow, "tissuename"]
    
    list(
      tissuename = tissuenames, 
      source = NA
    )
  })
  
  Tissues
}
