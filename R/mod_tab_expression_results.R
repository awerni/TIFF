expressionResultsTabUI_main <- function(id){
  ns <- NS(id)
  tabLayoutUI_main(ns("tab"))
}

expressionResultsTabUI_sidebar <- function(id){
  ns <- NS(id)
  
  list(
    tabLayoutUI_sidebar(
      id = ns("tab"),
      defaults = list(
        left = "roc",
        middle = "point",
        right = "violin"
      )
    ),
    hr(),
    extendedTableDownloadUI(ns("download"), "expr")
  )
}

expressionResultsTab <- function(input, output, session, classLabel, classSelection, Results, TableData, species){
  ExpressionGene <- reactive({
    s <- SelectedRow()
    validate(need(s, "no gene selected"))
    
    df <- TableData()
    as.list(df[s, c("ensg", "symbol")])
  })
  
  TpmTissue <- reactive({
    ensg <- ExpressionGene()$ensg
    cs <- Results()$cs
    
    getTissueDataGeneExpressionById(ensg, cs)
  })
  
  PrettyTableData <- reactive({
    df <- TableData()
    req(df)
    
    df %>% mutate(location = getEnsemblLocationLink(location))
  })
  
  plotFun <- function(plotType){
    gene <- ExpressionGene() # validate selection first
    
    res <- Results()
    cl <- reactiveValuesToList(classLabel)
    
    generateExpressionPlot(
      df = TpmTissue(), 
      ca = makeClassAssignment(res$cs, cl), 
      plotType = plotType,
      title = paste0("\n", gene$symbol, " - ", gene$ensg)
    )
  }
  
  rowInfo <- callModule(
    module = tabLayout,
    id = "tab",
    plotFun = plotFun, 
    TableData = PrettyTableData,
    jsRowCallback = getEnsgRowCallback(species)
  )
  SelectedRow <- rowInfo$SelectedRow
  
  callModule(
    module = extendedTableDownload,
    id = "download",
    Table = PrettyTableData,
    Subject = TpmTissue,
    Item = ExpressionGene,
    classSelection = classSelection,
    classLabel = classLabel,
    by = "ensg",
    additional = "symbol"
  )
  
  rowInfo$AllRows
}
