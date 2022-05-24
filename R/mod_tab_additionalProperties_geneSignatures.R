additionalPropertiesGeneSignaturesTabUI_main <- function(id){
  ns <- NS(id)
  tabLayoutUI_main(ns("tab"))
}

additionalPropertiesGeneSignaturesTabUI_sidebar <- function(id){
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
    extendedTableDownloadUI(ns("download"), "score", "signature")
  )
}

additionalPropertiesGeneSignaturesTab <- function(input, output, session, classSelection, classLabel, TissueAnnotationFocus){
  # Reactives -----------------------------------------------------------------
  Result_signature <- reactive({
    validate(
      need(classSelection$class1, "goto input tab and select tissues for class1"),
      need(classSelection$class2, "goto input tab and select tissues for class2")
    )
    cs <- reactiveValuesToList(classSelection)
    df <- differentialGeneSignature(cs, p = session) 
    validateFunctionResult(df)
    df <- df %>% filter(!is.na(higher))
    validate(need(nrow(df) > 0, "not enough data available"))
    
    df %>%
      mutate_at(c("logFC", "P.Value", "adj.p.val"), signif, 3) %>% 
      select(signature, higher, logFC, P.Value, adj.p.val)
  })
  
  Signature <- reactive({
    s <- SelectedRow()
    validate(need(s, "no signature selected"))
    
    df <- TableData()
    as.list(df[s, "signature", drop = FALSE])
  })
  
  SignatureScoreTissue <- reactive({
    signature <- Signature()$signature
    cs <- reactiveValuesToList(classSelection)
    getTissueDataGeneSignatureById(signature, cs)
  })
  
  # Layout --------------------------------------------------------------------
  plotFun <- function(plotType) {
    signature <- Signature()
    
    cs <- reactiveValuesToList(classSelection)
    cl <- reactiveValuesToList(classLabel)
    
    generateGeneSignaturePlot(
      df = SignatureScoreTissue(), 
      ca = makeClassAssignment(cs, cl), 
      plotType = plotType,
      title = paste0("\n", signature$signature)
    )
  }
  
  TableData <- reactive({
    validate(
      need(Result_signature(), "gene signature test has not been calculated")
    )
    cl <- reactiveValuesToList(classLabel)
    Result_signature() %>% 
      mutate(higher = styleHigherCol(higher == "class1", cl)) %>%
      left_join(geneSignatures, by = "signature") %>%
      mutate(
        description = paste0('<a href="', hyperlink, '" target="_blank">', description, '</a>'), 
        hyperlink = NULL
      ) %>%
      select(signature, description, logFC, higher, everything())
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
    Subject = SignatureScoreTissue,
    Item = Signature,
    classSelection = classSelection,
    classLabel = classLabel,
    by = "signature"
  )
}
