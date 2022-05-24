patientAnnotationInputModeUI <- function(id){
  ns <- NS(id)
  annotationInputUI(ns("handler"))
}

patientAnnotationInputMode <- function(input, output, session, TissuePrefilter, PatientAnnotationFuller){
  Table <- reactive({
    withProgress(
      expr = {
        tp <- TissuePrefilter()
        pa <- getPatientAnnotationFilter(tp, PatientAnnotationFuller())
        
        validate(need(pa, "no patient data available"))
        db_col <- sym(tp$db_col)
        pa %>% arrange(!!db_col)
      },
      message = "Retrieving tissue annotation", 
      value = 0.3
    )
  })
  
  PatientAnnotation_ti <- callModule(
    module = annotationInput,
    id = "handler",
    Table = Table
  )
  
  PatientAnnotation_ti
}
