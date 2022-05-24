copyNumberTabUI <- function(id){
  ns <- NS(id)
  
  fluidRow(
    column_2(copyNumberTabUI_sidebar(id)),
    column_10(tabLayoutUI_main(ns("tab")))
  )
}

copyNumberTabUI_sidebar <- function(id){
  ns <- NS(id)
  
  ret <- div(
    tabLayoutUI_sidebar(
      id = ns("tab"),
      defaults = list(
        left = "roc",
        middle = "point",
        right = "scatterplot"
      ),
      additionalChoices = c("scatter plot" = "scatterplot")
    ),
    hr(),
    extendedTableDownloadUI(ns("download"), "copynumber")
  )
  
  div(
    uiOutput(ns("sidebar_run")),
    hr(),
    ret
  )
}

copyNumberTab <- function(input, output, session, fm, classSelection, classLabel, gene_anno, species){
  # Sidebar -------------------------------------------------------------------
  ns <- session$ns
  output$sidebar_run <- renderUI({
    fmButtonBothClasses(ns("run"), fm, classSelection)
  })
  
  # Reactives -----------------------------------------------------------------
  CopynumberGene <- reactive({
    s <- SelectedRow()
    validate(need(s, "no gene selected"))
    
    df <- TableData()
    as.list(df[s, c("ensg", "symbol")])
  })
  
  RelativeCopynumberTissue <- reactive({
    ensg <- CopynumberGene()$ensg
    cs <- Results()$cs
    
    getTissueDataCopyNumberById(ensg, cs)
  })
  
  ExpressionCopynumberTissue <- reactive({
    ensg <- CopynumberGene()$ensg
    cs <- Results()$cs
    
    getTissueDataGeneExpressionById(ensg, cs)
  })
  
  # FutureManager -------------------------------------------------------------
  Results <- reactiveVal()
  Args <- reactive({
    list(
      cs = reactiveValuesToList(classSelection),
      anno = gene_anno()
    )
  })
  
  fm$registerRunObserver(
    inputId = ns("run"),
    label = "Copy Number",
    statusVar = Results,
    longFun = copyNumberLongFun,
    Args = Args
  )
  
  resultClassLabel <- registerFreezedClassLabel(output, classLabel, Results, fm, ns("run"))
  
  # Layout --------------------------------------------------------------------
  plotFun <- function(plotType){
    gene <- CopynumberGene() # validate selection first
    
    res <- Results()
    cl <- reactiveValuesToList(resultClassLabel)
    
    df <- RelativeCopynumberTissue()
    ca <- makeClassAssignment(res$cs, cl)
    title <- paste0("\n", gene$symbol, " - ", gene$ensg)
    
    if (plotType == "scatterplot"){
      generateCopynumberExpressionPlot(
        dfCN = df,
        dfExpr = ExpressionCopynumberTissue(),
        ca = ca,
        title = title
      )
    } else {
      generateCopynumberPlot(
        df = df,
        ca = ca,
        plotType = plotType,
        title = title
      )
    }
  }
  
  TableData <- reactive({
    res <- Results()
    fmValidate(res)
    
    cl <- reactiveValuesToList(resultClassLabel)
    
    res$df %>% 
      mutate(higher = styleHigherCol(log2FC > 0, cl)) %>%
      mutate_at(c("class1.median", "class2.median", "p.value", "log2FC", "adj.P.Value"), signif, 3) %>%
      rename_at(vars(contains('class1')), list(~ gsub("^class1", cl$class1_name, .))) %>%
      rename_at(vars(contains('class2')), list(~ gsub("^class2", cl$class2_name, .))) %>%
      mutate(location = getEnsemblLocationLink(location)) %>%
      relocate(higher, .after = location)
  })
  
  rowInfo <- callModule(
    module = tabLayout,
    id = "tab",
    plotFun = plotFun, 
    TableData = TableData,
    jsRowCallback = getEnsgRowCallback(species)
  )
  SelectedRow <- rowInfo$SelectedRow
  
  # Download ------------------------------------------------------------------
  callModule(
    module = extendedTableDownload,
    id = "download",
    Table = TableData,
    Subject = RelativeCopynumberTissue,
    Item = CopynumberGene,
    classSelection = classSelection,
    classLabel = resultClassLabel,
    by = "ensg",
    additional = "symbol"
  )
}

copyNumberLongFun <- function(task, cs, anno){
  df <- getCopyNumberAssociation(cs, anno, p = task)
  
  if (is.null(df)) return() # handle the task cancel
  if (is.fmError(df)) return(df)
  
  list(
    df = df,
    cs = cs
  )
}
