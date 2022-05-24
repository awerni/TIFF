#' Get TIFF settings
#' 
#' Retrieves settings from an yml file. By default, the file from the package 
#' is used. If you want to use your custom file, please set R_CONFIG_FILE 
#' environmental variable.
#' 
#' If the `obj` argument is not NULL, it will be returned from the function.
#' 
#' @param obj NULL or the settings object
#' @export
getSettings <- function(obj = NULL){
  defaultPath <- system.file("shiny/config.yml", package = "TIFF")
  path <- Sys.getenv("R_CONFIG_FILE", unset = defaultPath)
  
  if (is.null(obj)){
    config::get(file = path)
  } else {
    obj
  }
}

#' Check the integrity of the settings object
#'
#' @param settings settings object
#'
#' @return TRUE if settings object is valid.
#' @export
#'
#' @examples
#' 
#' config <- getSettings()
#' initialSettingsCheck(config)
#' 
initialSettingsCheck <- function(settings){
  # 1) required fields check
  reqFields <- c("db", "links", "species")
  missingFields <- setdiff(reqFields, names(settings))
  if (length(missingFields) > 0){
    message("Missing setting fields!")
    message(paste(missingFields, collapse = ", "))
    
    return(structure(FALSE, reason = "opts"))
  }
  
  # 2) DB settings check
  reqFields <- c("host", "name", "user", "password", "port")
  missingFields <- setdiff(reqFields, names(settings[["db"]]))
  if (length(missingFields) > 0){
    message("Missing DB setting fields!")
    message(paste(missingFields, collapse = ", "))
    
    return(structure(FALSE, reason = "opts"))
  }
  
  # 3) Links check
  if (is.null(settings[["links"]][["msigDBLink"]])){
    message("Missing link setting fields!")
    message("msigDBLink")
    
    return(structure(FALSE, reason = "opts"))
  }
  
  if (is.null(settings[["links"]][["docuLink"]])){
    message("Missing link setting fields!")
    message("docuLink")
    
    return(structure(FALSE, reason = "opts"))
  }
  
  TRUE
}

#' Utility function for making the text callback function.
#'
#' @param scoreTemplate glue template which will render the text.
#' @param sigx indicating the number of decimal places for rx.
#' @param sigy indicating the number of decimal places for ry.
#'
#' @return
#' @export
#'
#' @examples
#' 
#' cl <- makeTextCallback( "Cxpression score from {ry[2]} to {ry[1]}")
#' cl(100, c(0,10), c(0, 100))
#' cl(100, c(0,10), c(0, NA))
#' 
makeTextCallback <- function(scoreTemplate, sigx = 4, sigy = 4) {
  
  function(n, rx, ry){
    
    label <-  glue::glue("{get_label(n)} selected.")
    
    if(!anyNA(ry)) {
      ry <- signif(ry, sigx)
      rx <- signif(rx, sigy)
      label <- paste(label, scoreTemplate)
      label <- glue::glue(label)
    }
    
    return(label)
  }
  
}
