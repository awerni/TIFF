tissueSamplePrefilterUI <- function(id) {
  ns <- NS(id)
  tissue_groups <- getTissueGroups()
  radioId <- ns("tissue_type")
  
  createCondition <- function(value){
    paste0("input['", radioId, "'] == '", value, "'")
  }
  
  tagList(
    radioButtons(
      inputId = radioId, 
      label = "Tissue type prefilter", 
      choices = c("tumor", "adjacent normal", "normal", "PDX"), 
      selected = "tumor"
    ),
    shinyjs::hidden(div(
      id = ns("panel"),
      conditionalPanel(
        condition = createCondition("tumor"),
        selectInput(
          inputId = ns("tumor_selector"), 
          label = "Tumor selection", 
          multiple = TRUE, 
          selectize = FALSE, 
          size = 6, 
          choices = tissue_groups$tumor, 
          selected = tissue_groups$tumor
        )
      ),
      conditionalPanel(
        condition = createCondition("adjacent normal"),
        selectInput(
          inputId = ns("adjacent_normal_selector"), 
          label = "Adjacent Tumor Selection", 
          multiple = TRUE, 
          selectize = FALSE, 
          size = 6,
          choices = tissue_groups$adjacent_normal, 
          selected = tissue_groups$adjacent_normal
        )
      ),
      conditionalPanel(
        condition = createCondition("normal"),
        selectInput(
          inputId = ns("normal_selector"), 
          label = "Organ selection", 
          multiple = TRUE, 
          selectize = FALSE, 
          size = 6, 
          choices = tissue_groups$normal, 
          selected = tissue_groups$normal
        )
      ),
      conditionalPanel(
        condition = createCondition("PDX"),
        selectInput(
          inputId = ns("pdx_selector"), 
          label = "PDX tumor selection", 
          multiple = TRUE, 
          selectize = FALSE, 
          size = 6, 
          choices = tissue_groups$pdx, 
          selected = tissue_groups$pdx
        )
      )
    )),
    uiOutput(ns("indicator"))
  )
}

tissueSamplePrefilter <- function(input, output, session) {
  output$indicator <- renderUI({
    shinyjs::show("panel")
    NULL
  })
  
  tissue_type <- reactive({
    validate(need(input$tissue_type, message = FALSE))
    input$tissue_type
  })
  
  tumor_selector <- reactive({
    validate(need(input$tumor_selector, message = FALSE))
    input$tumor_selector
  })
  
  adjacent_normal_selector <- reactive({
    validate(need(input$adjacent_normal_selector, message = FALSE))
    input$adjacent_normal_selector
  })
  
  normal_selector <- reactive({
    validate(need(input$normal_selector, message = FALSE))
    input$normal_selector
  })
  
  pdx_selector <- reactive({
    validate(need(input$pdx_selector, message = FALSE))
    input$pdx_selector
  })
  
  tissue_prefilter <- reactive({
    tissueType <- tissue_type()
    
    tissue <- switch(
      EXPR = tissueType,
      "tumor" = tumor_selector(),
      "adjacent normal" = adjacent_normal_selector(),
      "normal" = normal_selector(),
      "PDX" = pdx_selector()
    )
    
    makePrefilter(tissueType, tissue)

  })
  
  tissue_prefilter
}
