mutationPerGeneTabUI_main <- function(id){
  ns <- NS(id)
  tabLayoutUI_main(ns("tab"), FALSE)
}

mutationPerGeneTabUI_sidebar <- function(id){
  ns <- NS(id)
  
  list(
    tabLayoutUI_sidebar(
      id = ns("tab"),
      defaults = list(
        left = "bar",
        right = "pie"
      ),
      additionalChoices = c(
        "bar plot" = "bar", 
        "pie chart" = "pie"
      ),
      hidden = c("roc", "point", "violin", "box"),
      useMiddlePlot = FALSE
    ),
    hr(),
    extendedTableDownloadUI(ns("download"), "mutation")
  )
}

mutationPerGeneTab <- function(input, output, session, classSelection, classLabel, gene_anno, species){
  Result_mutation <- reactive({
    cs <- reactiveValuesToList(classSelection)
    validate(
      need(cs$class1, "goto input tab and select several tissues for class1"),
      need(cs$class2, "goto input tab and select several tissues for class2")
    )
    getMutationAssociation(cs, gene_anno(), p = session)
  })
  
  MutationGene <- reactive({
    s <- SelectedRow()
    validate(need(s, "no gene selected"))
    
    df <- TableData()
    as.list(df[s, c("ensg", "symbol")])
  })
  
  MutationTissue <- reactive({
    ensg <- MutationGene()$ensg
    cs <- reactiveValuesToList(classSelection)
    
    getTissueDataMutationById(ensg, cs)
  })
  
  plotFun <- function(plotType){
    cs <- reactiveValuesToList(classSelection)
    cl <- reactiveValuesToList(classLabel)
    ca <- makeClassAssignment(cs, cl)
    generateMutationPlot(
      df = MutationTissue(),
      plotType = plotType,
      ca = ca, 
      ensg = MutationGene()$ensg,
      geneSymbol = MutationGene()$symbol
    )
  }
  
  TableData <- reactive({
    rm <- Result_mutation()
    validateFunctionResult(rm)
    
    cl <- reactiveValuesToList(classLabel)
    
    rm %>% 
      rename_at(vars(contains('class1')), list(~ gsub("^class1", cl$class1_name, .))) %>%
      rename_at(vars(contains('class2')), list(~ gsub("^class2", cl$class2_name, .))) %>%
      mutate_at(c("p.value", "adj.P.Value"), signif, 3) %>%
      mutate(higher = styleHigherCol(higher == "class1", cl))
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
    Subject = MutationTissue,
    Item = MutationGene,
    classSelection = classSelection,
    classLabel = classLabel,
    by = "ensg",
    additional = "symbol"
  )
}
