library(shiny)
library(shinydashboard)
library(shinyjs)

ui = shinyUI
( 
  fluidPage
  ( 
    
    ui <- dashboardPage
    (
      dashboardHeader(title = "supertest"),
      dashboardSidebar(
        sidebarMenu(
          menuItem("home", tabName = "home", icon = icon("home")),
          menuItem("results table", tabName = "results", icon = icon("th")),
          actionLink("done", "run", icon = icon("check-circle"))
        )
      )
      ,
      dashboardBody(
        tabItems(

          
          tabItem(tabName = "home",
                  # Boxes need to be put in a row (or column)
                  
                  fluidRow(
                    box(title = "Supertest",
                        textOutput("grpTxt")),
                    
                    box(title = "General Settings",
                        sliderInput("seqhom", "Minimal Sequence Homolgy:",0,1, 0.9,),
                        sliderInput("p0", "Minimal confidence:", 0,1,0.25),
                        sliderInput("thr", "Specificity threshold:", 0, 0.08, 0.018),
                        sliderInput("pss", "Minimum peptide set size:", 0, 6, 4)
                    ),
                    fluidRow(
                      
                      box(
                        title = "In Vivo In Vitro DB settings",
                        sliderInput("iviv", "Confidence:", 0,1, 0.9)
                      ),
                      box(
                        title = "PhosphoNET DB settings",
                        sliderInput("pn.breakpoints", "Rank breakpoints (define mid range): ", 1,25,value = c(4,12)),
                        sliderInput("pn.low", "Low rank confidence", 0,1,0.9),
                        sliderInput("pn.mid", "Mid rank confidence", 0,1,0.7),
                        sliderInput("pn.high", "Low rank confidence", 0,1,0.25)
                      )
                    )
                  )
          ),
          tabItem(tabName = "results",
                  fluidRow(
                    box(title = "Activity up",
                        textOutput("uptext"),
                        tableOutput("kup")),
                    box(title = "Activity down",
                        textOutput("downtext"),
                        tableOutput("kdn"))
                    
                  )
          )
          
          
        )
        
      )
    ),
    shinyjs::useShinyjs()
  )
)
