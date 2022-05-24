tissueAnnotationInputModeUI <- function(id){
  ns <- NS(id)
  annotationInputUI(ns("handler"))
}

tissueAnnotationInputMode <- function(input, output, session, TissuePrefilter, TissueAnnotation){
  Table <- reactive({
    withProgress(
      expr = {
        tp <- TissuePrefilter()
        ta <- TissueAnnotation() %>% filter(get(tp$db_col) %in% tp$tissue)
        
        if (tp$db_col == "organ") {
          ta <- ta %>% filter(tumortype == "normal" & is.na(tumortype_adjacent))
        }
        
        db_col <- sym(tp$db_col)
        ta %>% 
          filter(tissuepanel == tp$tissue_panel) %>%
          select(tp$db_cols) %>%
          arrange(!!db_col) %>%
          droplevels()
      },
      message = "Retrieving tissue annotation", 
      value = 0.3
    )
  })
  
  TissueAnnotation_ti <- callModule(
    module = annotationInput,
    id = "handler",
    Table = Table
  )
  
  TissueAnnotation_ti
}
