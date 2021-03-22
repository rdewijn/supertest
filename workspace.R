library(shiny)
library(tercen)
library(dplyr)
library(tidyr)


############################################

source("ui.R")
source("server.R")
options("tercen.workflowId"= "fa717b14ad3479b68291b372be004818")
options("tercen.stepId"= "91ed2bb9-5421-4b84-83b3-e2a193862498")

runApp(shinyApp(ui, server))  
