library(shiny)
library(shinyBS)
library(dplyr)
library(XIFF)
library(TIFF)
library(logger)

log_threshold(TRACE)

# use R_CONFIG_FILE environmental variable to use your own config (based on config.yml)
settings <- getSettings() 

# check on the app restart, will return an error: "App failed to start" in the browser
if (!initialSettingsCheck(settings)) stop("Initial check failed!")

strategy <- getFmStrategy(settings)
FutureManager::plan(strategy, substitute = FALSE)
setDbOptions(settings)

geneSignatures <- getAvailableGeneSignatures()
hallmarkGeneSets <- getAvailableHallmarkGeneSets()
