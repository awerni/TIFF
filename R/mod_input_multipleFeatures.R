multipleFeaturesInputModeUI <- function(id, defaultGene = "KRAS"){
  ns <- NS(id)
  
  choices <- c(
    "Gene expression" = "expr", 
    "Copy number alteration" = "copy"
  )
  
  list(
    fluidRow(
      column_3(
        tags$h4("X axis"),
        textInput(
          inputId = ns("gene_x"), 
          label = "gene symbol", 
          value = defaultGene
        ),
        selectInput(
          inputId = ns("source_x"), 
          label = "data type", 
          selected = "expr", 
          choices = choices
        )
      ),
      column_3( 
        tags$h4("Y axis"),
        textInput(
          inputId = ns("gene_y"), 
          label = "gene symbol", 
          value = "",
          placeholder = "the same as for X axis"
        ),
        selectInput(
          inputId = ns("source_y"), 
          label = "data type", 
          selected = "copy", 
          choices = choices
        )
      ),
      column_3(
        tags$h4("Mutation"),
        textInput(
          inputId = ns("gene_mut"), 
          label = "gene symbol", 
          value = ""
        ),
        shinyjs::hidden(radioButtons(
          inputId = ns("mutated"),
          label = "Show by mutation",
          choices = c(
            "all" = "all", 
            "mutated" = "mut", 
            "wild-type" = "wt"
          ),
          selected = "all"
        ))
      ),
      column_3()
    ),
    fluidRow_12(
      hr(style = "margin-top:5px; margin-bottom:5px;"),
      radioButtons(
        inputId = ns("plot_type"),
        label = NULL, 
        choices = c(
          "tissue selection" = "selection", 
          "data exploration" = "exploration"
        ),
        selected = "selection", 
        inline = TRUE
      )
    ),
    fluidRow_12(
      uiOutput(ns("container"))
    )
  )
}

multipleFeaturesInputMode <- function(input, output, session, species, TissuePrefilter){
  MultipleFeaturesGeneX <- reactive({
    symbol <- input$gene_x
    if (is.null(symbol)) return()
    if (nchar(symbol) < 3) return()
    getGeneFromSymbol(symbol, species)
  })
  MultipleFeaturesDataX <- reactive({
    gene <- MultipleFeaturesGeneX()
    if (is.null(gene)) return()
    
    withProgress(
      expr = getMultipleFeaturesAxisData(input$source_x, gene, "x", TissuePrefilter()),
      value = 0.3,
      message = "fetching X axis data"
    )
  })
  
  MultipleFeaturesGeneY <- reactive({
    symbol <- input$gene_y
    if (invalid(symbol)){
      # get gene from X axis
      symbol <- input$gene_x
    }
    if (invalid(symbol) || nchar(symbol) < 3) return()
    getGeneFromSymbol(symbol, species)
  })
  MultipleFeaturesDataY <- reactive({
    gene <- MultipleFeaturesGeneY()
    if (is.null(gene)) return()
    
    withProgress(
      expr = getMultipleFeaturesAxisData(input$source_y, gene, "y", TissuePrefilter()),
      value = 0.3,
      message = "fetching Y axis data"
    )
  })
  
  MultipleFeaturesGeneMut <- reactive({
    symbol <- input$gene_mut
    inheritsX <- invalid(symbol)
    
    if (inheritsX){
      # get gene from X axis
      symbol <- input$gene_x
    }
    
    isInvalid <- invalid(symbol) || nchar(symbol) < 3
    
    inputId <- "mutated"
    if (inheritsX || isInvalid){
      shinyjs::hide(id = inputId)
    } else {
      shinyjs::show(id = inputId)
    }
    
    if (!isInvalid){
      getGeneFromSymbol(symbol, species)
    }
  })
  MultipleFeaturesDataMut <- reactive({
    geneMut <- MultipleFeaturesGeneMut()
    if (is.null(geneMut)) return()
    dataX <- MultipleFeaturesDataX()
    if (is.null(dataX$d) || nrow(dataX$d) == 0) return()
    dataY <- MultipleFeaturesDataY()
    if (is.null(dataY$d) || nrow(dataY$d) == 0) return()
    
    fullData <- dataY$d
    duplicatedNames <- setdiff(intersect(names(fullData), names(dataX$d)), "tissuename")
    if (length(duplicatedNames) > 0){
      fullData <- fullData %>% select(-any_of(duplicatedNames))
    }
    fullData <- fullData %>% inner_join(dataX$d, by = "tissuename")
    
    getTissueDataMutationById(
      ensg = geneMut$ensg,
      tissueClasses = list(class1 = fullData$tissuename)
    ) %>%
      select(tissuename, aamutated) %>%
      right_join(fullData, by = "tissuename")
  })
  
  #FIXME: export as a plot function
  MultipleFeaturesPlot <- reactive({
    
    fullData <- MultipleFeaturesDataMut()
    
    dataX <- isolate(MultipleFeaturesDataX())
    dataY <- isolate(MultipleFeaturesDataY())
    validate(
      need(fullData, "data not available"),
      need(dataX, "data not available"),
      need(dataY, "data not available")
    )
    
    xlim <- range(fullData$xVar, na.rm = TRUE)
    ylim <- range(fullData$yVar, na.rm = TRUE)
    
    mutated <- input$mutated
    if (mutated != "all"){
      fullData <- filter(fullData, aamutated == mutated)
    }
    
    scaleX <- get(dataX$scale)
    scaleY <- get(dataY$scale)
    
    if (!nzchar(isolate(input$gene_mut))){
      colourAes <- TissuePrefilter()$db_col
      colourAes <- sym(colourAes)
      colourScale <- NULL
    } else {
      colourAes <- sym("aamutated")
      colourScale <- scale_colour_manual(
        na.value = "gray60",
        values = c(
          "mut" = "#F8766D",
          "wt" = "#00BFC4"
        )
      )
    }
    
    fullData <- fullData %>% mutate(colourLabel = !!colourAes)
    
    colourAes <- enquo(colourAes)
    mapping <- aes(
      x = xVar,
      y = yVar,
      colour = !!colourAes,
      tissuename = tissuename,
      colourLabel = colourLabel
    )
    
    p <- ggplot(
      data = fullData,
      mapping = mapping
    ) + 
      geom_point() + 
      scaleX() + 
      scaleY() +
      xlab(dataX$varLabel) + 
      ylab(dataY$varLabel) + 
      theme_bw() + 
      theme(legend.position = "bottom", text = element_text(size = 15)) +
      coord_cartesian(xlim = xlim, ylim = ylim)
    
    if (!is.null(colourScale)){
      p <- p + colourScale
    }
    
    p
  })
  
  MultipleFeaturesPlotCheck <- reactive({
    geneX <- MultipleFeaturesGeneX()
    geneY <- MultipleFeaturesGeneY()
    geneMut <- MultipleFeaturesGeneMut()
    sourceX <- input$source_x
    sourceY <- input$source_y
    
    !is.null(geneX) && !is.null(geneY) && !is.null(geneMut) &&
      !is.null(sourceX) && !is.null(sourceY) &&
      !is.null(input$mutated) &&
      !is.null(TissuePrefilter())
  })
  
  ns <- session$ns
  output$container <- renderUI({
    plotType <- input$plot_type
    validate(need(plotType, "select plot type"))
    
    if (plotType == "selection"){
      brushPlotUI(
        id = ns("plot"), 
        height = "600px", 
        direction = "xy"
      )
    } else if (plotType == "exploration"){
      plotWrapperUI(
        id = ns("tooltip_plot"), 
        height = "600px"
      )
    }
  })
  
  callModule(
    module = plotWrapper,
    id = "tooltip_plot",
    PlotExpr = MultipleFeaturesPlot,
    varDict = list(),
    PlotType = TRUE,
    tooltipCallback = function(x) {
      paste0(x$tissuename, "\n", x$colourLabel)
    }
  )
  
  MultipleFeatures_cl <- callModule(
    module = brushPlot,
    id = "plot",
    plotExpr = MultipleFeaturesPlot,
    checkExpr = MultipleFeaturesPlotCheck,
    defaultCutoff_x = function() { getCutoffValue(input$source_x) },
    defaultCutoff_y = function() { getCutoffValue(input$source_y) },
    dx = 0,
    dy = 0,
    availableChoices = "cutoff"
  )
  
  MultipleFeatures_cl
}

getCutoffValue <- function(dataSource, default = 0) {
  switch(
    EXPR = dataSource,
    expr = 10,
    avana = -1,
    copy = 1.5,
    default
  )
}
