library(shiny)
library(tercen)
library(dplyr)
library(tidyr)


############################################

source("ui.R")
source("server.R")
options("tercen.workflowId"= "3a85f3a46cd4009b31373508220064f2")
options("tercen.stepId"= "86799058-7709-4ec3-adc5-41a4ac6e3c11")

runApp(shinyApp(ui, server))  
