expressionHeatmapTabUI_main <- function(id){
  ns <- NS(id)
  fluidRow_12(
    br(), 
    plotWrapperUI(ns("heatmap"))
  )
}

expressionHeatmapTabUI_sidebar <- function(id){
  ns <- NS(id)
  
  list(
    sliderInput(
      inputId = ns("ngenes"), 
      label = "number of diff. expressed genes", 
      min = 4, 
      max = 100, 
      step = 1, 
      value = 40
    ),
    checkboxInput(
      inputId = ns("flipaxes"), 
      label = "Flip Axes", 
      value = FALSE
    ),
    radioButtons(
      inputId = ns("scale"), 
      label = "Scale values", 
      choices = c(
        "none" = "none", 
        "genes" = "genes", 
        "tissues" = "tissues"
      ), 
      selected = "genes"
    ),
    selectInput(
      inputId = ns("label"), 
      label = "Sample text label", 
      choices = NULL,
      selectize = FALSE
    ),
    selectInput(
      inputId = ns("graphical_label"), 
      label = "Graphical label", 
      choices = NULL,
      selectize = FALSE
    ),
    uiOutput(ns("indicator"))
  )
}

expressionHeatmapTab <- function(input, output, session, classLabel, Results, TableData, AllRows, TissueAnnotationFocus){
  output$indicator <- renderUI({
    anno <- TissueAnnotationFocus()
    x <- colnames(anno)
    
    updateSelectInput(
      session = session,
      inputId = "label",
      choices = x,
      selected = isolate(input$label) %||% x[1]
    )
    
    updateSelectInput(
      session = session,
      inputId = "graphical_label",
      choices = setdiff(c("none", x), "tissuename"),
      selected = isolate(input$graphical_label) %||% "none"
    )
    
    NULL
  })
  
  HeatmapExpr <- reactive({
    res <- Results()
    fmValidate(res)
    
    nGenes <- input$ngenes
    labelCol <- input$label
    graphicalCol <- input$graphical_label
    req(
      TableData(),
      nGenes,
      labelCol,
      graphicalCol,
      input$scale,
      !is.null(input$flipaxes)
    )
    
    resultData <- TableData() %>% head(nGenes)
    
    sAnno <- res$TissueAnnotationFocus
    
    labelCol <- setdiff(labelCol, "none")
    graphicalCol <- setdiff(graphicalCol, "none")
    selcols <- c(labelCol, graphicalCol)
    if (length(selcols) != length(intersect(selcols, colnames(sAnno)))) return()
    
    cs <- res$cs
    if (!checkClassSelection(cs)) return()
    cl <- reactiveValuesToList(classLabel)

    generateExpressionHeatmap(
      genes = resultData[["ensg"]],
      sampleClasses = cs, 
      classLabel = cl,
      sampleLabel = labelCol,
      annoVariables = graphicalCol,
      heat_scale =  input$scale, 
      flipAxes = input$flipaxes, 
      sampleAnno = sAnno
    )
  })
  
  widthFun <- function() {
    res <- Results()
    fmValidate(res)
    
    cs <- res$cs
    if (is.null(cs)) return("auto")
    ti <- c(cs$class1, cs$class2)
    labels <- reactiveValuesToList(classLabel)
    classNames <- c(labels$class1_name, labels$class2_name)
    
    if (is.null(input$flipaxes) || is.null(input$ngenes)) return(300)
    
    wlabel <- 3 * max(nchar(classNames), na.rm = TRUE)
    if (input$flipaxes) {
      400 + length(ti) * 12 + wlabel # class legend reduces width
    } else {
      300 + input$ngenes * 12 + wlabel
    }
  }
  
  heightFun <- function() {
    res <- Results()
    fmValidate(res)
    
    cs <- res$cs
    if (is.null(cs)) return("auto")
    ti <- c(cs$class1, cs$class2)
    
    if (is.null(input$flipaxes) || is.null(input$ngenes)) return(300)
    
    if (input$flipaxes) {
      300 + input$ngenes * 12
    } else {
      300 + length(ti) * 12
    }
  }
  
  callModule(
    module = plotWrapper,
    id = "heatmap",
    PlotExpr = HeatmapExpr,
    width = widthFun,
    height = heightFun
  )
}
