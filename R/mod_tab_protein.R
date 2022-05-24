proteinTabUI <- function(id){
  ns <- NS(id)
  fluidRow(
    column_2(proteinTabUI_sidebar(ns)),
    column_10(tabLayoutUI_main(ns("tab")))
  )
}

proteinTabUI_sidebar <- function(ns){
  
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
    extendedTableDownloadUI(ns("download"), "score", "RPPA")
  )
}

proteinTab <- function(input, output, session, classSelection, classLabel, TissueAnnotationFocus){
  # Reactives -----------------------------------------------------------------
  ExpressionAntibody <- reactive({
    s <- SelectedRow()
    validate(need(s, "no antibody selected"))
    
    df <- TableData()
    as.list(df[s, c("antibody", "validation_status", "vendor")])
  })
  
  RPPAScoreTissue <- reactive({
    antibody <- ExpressionAntibody()$antibody
    cs <- reactiveValuesToList(classSelection)
    
    getTissueDataRppaById(antibody, cs)
  })
  
  AntibodyAnno <- reactive({
    getAntibodyInformation()
  })
  
  Result_protein <- reactive({
    validate(
      need(classSelection$class1, "goto input tab and select tissues for class1"),
      need(classSelection$class2, "goto input tab and select tissues for class2")
    )
    cs <- reactiveValuesToList(classSelection)
    df <- differentialProteinExpression(cs, p = session) 
    
    validateFunctionResult(df)
    
    ab_anno <- AntibodyAnno()
    
    df %>%
      left_join(ab_anno, by = "antibody") %>%
      select(antibody, validation_status, vendor, catalog_number, higher, logFC, P.Value, adj.p.val)
  })
  
  # Layout --------------------------------------------------------------------
  plotFun <- function(plotType) {
    antibody <- ExpressionAntibody()
    
    cs <- reactiveValuesToList(classSelection)
    cl <- reactiveValuesToList(classLabel)
    
    generateProteinExpressionPlot(
      df = RPPAScoreTissue(), 
      ca = makeClassAssignment(cs, cl), 
      plotType = plotType,
      title = paste0("\n", antibody$antibody, " - ", antibody$validation_status)
    )
  }
  
  TableData <- reactive({
    res <- Result_protein()
    validate(need(res, "protein expression has not been calculated"))
    
    cl <- reactiveValuesToList(classLabel)
    
    res %>% 
      mutate_at(c("logFC", "P.Value", "adj.p.val"), signif, 3) %>%
      mutate(higher = styleHigherCol(higher == "class1", cl)) %>%
      select(antibody, validation_status, vendor, catalog_number, higher, everything())
  })
  
  rowInfo <- callModule(
    module = tabLayout,
    id = "tab",
    plotFun = plotFun, 
    TableData = TableData
  )
  SelectedRow <- rowInfo$SelectedRow
  
  # Download ------------------------------------------------------------------
  callModule(
    module = extendedTableDownload,
    id = "download",
    Table = TableData,
    Subject = RPPAScoreTissue,
    Item = ExpressionAntibody,
    classSelection = classSelection,
    classLabel = classLabel,
    by = "antibody",
    additional = c("vendor", "validation_status")
  )
}
