derivedDataTabUI <- function(id){
  ns <- NS(id)
  
  tabsetPanel(
    id = ns("tabset"),
    columnTabPanel(
      title = "immune cells",
      value = "immune cells",
      derivedDataImmuneCellsTabUI_sidebar(ns("immune_cells")),
      derivedDataImmuneCellsTabUI_main(ns("immune_cells"))
    ),
    columnTabPanel(
      title = "metabolics",
      value = "metabolics",
      derivedDataMetabolicsTabUI_sidebar(ns("metabolics")),
      derivedDataMetabolicsTabUI_main(ns("metabolics"))
    ),
    columnTabPanel(
      title = "signaling",
      value = "signaling",
      derivedDataSignalingTabUI_sidebar(ns("signaling")),
      derivedDataSignalingTabUI_main(ns("signaling"))
    ),
    columnTabPanel(
      title = "clones",
      value = "clones",
      derivedDataClonesTabUI_sidebar(ns("clones")),
      derivedDataClonesTabUI_main(ns("clones"))
    )
  )
}


derivedDataTab <- function(input, output, session, classSelection, classLabel){
  callModule(
    module = derivedDataImmuneCellsTab,
    id = "immune_cells",
    classSelection = classSelection, 
    classLabel = classLabel
  )
  
  callModule(
    module = derivedDataMetabolicsTab,
    id = "metabolics",
    classSelection = classSelection, 
    classLabel = classLabel
  )
  
  callModule(
    module = derivedDataSignalingTab,
    id = "signaling",
    classSelection = classSelection, 
    classLabel = classLabel
  )
  
  callModule(
    module = derivedDataClonesTab,
    id = "clones",
    classSelection = classSelection, 
    classLabel = classLabel
  )
}
