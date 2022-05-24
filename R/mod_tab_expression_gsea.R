expressionGseaTabUI_main <- function(id){
  ns <- NS(id)
  div(
    fluidRow_12( 
      br(), 
      plotWrapperUI(ns("plot"))
    ),
    fluidRow_12(containerDT(ns("table")))
  )
}

expressionGseaTabUI_sidebar <- function(id){
  ns <- NS(id)
  
  list(
    radioButtons(
      inputId = ns("rank"), 
      label = "GSEA ranking", 
      choices = c(
        "direction + p-value" = "p.valueDir", 
        "logFC" = "logFC"
      ), 
      selected = "p.valueDir"
    ),
    radioButtons(
      inputId = ns("gene_set"), 
      label = "gene sets",
      choices = NA
    ),
    hr(),
    textInput(
      inputId = ns("filename"), 
      label = "Choose file name", 
      value = "GSEA-filename"
    ),
    downloadButton(
      outputId = ns("download"), 
      label = "Download"
    ),
    uiOutput(ns("indicator"))
  )
}

expressionGseaTab <- function(input, output, session, fm, classLabel, Results, species, msigDBLink){
  output$indicator <- renderUI({
    updateRadioButtons(
      session = session,
      inputId = "gene_set",
      choices = getGSEACollection(),
      selected = "hallmark"
    )
    
    NULL
  })
  
  ResultsGSEA <- reactiveVal()
  reRunRequired <- FALSE
  
  Args <- reactive({
    res <- Results()
    fmValidate(res)
    
    reRunRequired <<- TRUE
    list(
      geneSet = input$gene_set,
      res = res,
      rankType = input$rank,
      species = species
    )
  })
  
  ArgsTrigger <- reactive({
    hiddenId <- paste0("output_", session$ns("table"), "_hidden")
    list(
      args = Args(),
      isHidden = session$clientData[[hiddenId]]
    )
  })
  
  taskId <- "gsea_data" # some initial ID
  observeEvent(
    eventExpr = ArgsTrigger(),
    handlerExpr = {
      d <- ArgsTrigger()
      if (d$isHidden || !reRunRequired) return()
      reRunRequired <<- FALSE
      
      args <- d$args
      req(validateArgs(args))
      
      fm$cancel(taskId)
      
      # generate new ID in case of changing args
      # new process with the latest args will spin up
      # the previous process with outdated args will be canceled
      taskId <<- fmGenerateTaskId("gsea_data") 
      fm$showProgress(taskId, "GSEA", ResultsGSEA)
      fm$run(
        taskId = taskId,
        fun = gseaLongFun,
        args = args,
        statusVar = ResultsGSEA
      )
    }
  )
  
  Result_expressionGSEA <- reactive({
    validate(need(input$rank, "wait"))
    res <- ResultsGSEA()
    fmValidate(res)
    
    res$df %>%
      mutate_at(c("stat.mean", "p.val", "q.val"), signif, 3) %>%
      mutate(geneset = paste0('<a href="', msigDBLink, geneset, '" target="_blank">', geneset, '</a>')) %>%
      select(geneset, set.size, stat.mean, p.val, q.val)
  })
  
  GSEA_SelectedGeneSet <- reactive({
    s <- input$table_rows_selected
    validate(need(s, "no gene set selected"))
    
    gs <- unlist(Result_expressionGSEA()[s, "geneset"])
    gsub("</a>", "", gsub("<a href.*blank\\\">", "", gs)) # remove a href
  })
  
  GseaPlotExpr <- reactive({
    gs <- GSEA_SelectedGeneSet()
    validate(
      need(input$rank, "wait"),
      need(gs, "select gene set")
    )
    res <- ResultsGSEA()
    fmValidate(res)
    
    ensg_gs <- res$gsea_data[gs]
    cl <- reactiveValuesToList(classLabel)
    classNames <- c(cl$class1_name, cl$class2_name)
    generateGSEA_plot(res$res$df, ensg_gs, input$rank, classNames)
  })
  
  callModule(
    module = plotWrapper,
    id = "plot",
    PlotExpr = GseaPlotExpr
  )
  
  output$table <- DT::renderDataTable(
    expr = Result_expressionGSEA(),
    rownames = FALSE,
    filter = "top",
    options = list(
      pageLength = 10, 
      search = list(regex = TRUE)
    ),
    escape = FALSE,
    selection = list(
      mode = "single", 
      selected = 1, 
      target = "row"
    )
  )
  
  output$download <- downloadHandler(
    filename = function() { 
      paste0(input$filename, ".csv")
    },
    content = function(file) {
      res <- Result_expressionGSEA() %>% 
        mutate(geneset = gsub("<[^>]*>", "", geneset))
      
      readr::write_excel_csv(res, file, na = "")
    }
  )
}

gseaLongFun <- function(task, geneSet, res, rankType, species){
  if (!is.null(task)) {
    if (fmIsInterrupted(task)) return()
    fmUpdateProgress(task, progress = 0.05, msg = "fetching data for GSEA...")
  }
  
  geneSets <- if (geneSet == "hallmark") {
    res$gsea_data_hallmark
  } else {
    getGSEAdata(species, geneSet)
  }
  
  df <- getGSEA(
    diffExResult = res$df,
    geneSets = geneSets,
    rankType = rankType,
    p = task
  )
  
  list(
    gsea_data = geneSets,
    df = df,
    res = res
  )
}
