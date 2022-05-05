## Author: Jason Yeung
## jy357@bu.edu
## BF591 Final Project

library(shiny)
library(dplyr)
library(ggplot2)
library(colourpicker)

ui <- fluidPage(
  titlePanel("BF591 Final Project"),
  
  mainPanel(
    tabsetPanel(type = "tabs",
                tabPanel("Samples",
                         sidebarLayout(
                           sidebarPanel(
                             fileInput("csv_file", paste0("Sample File (.csv)")),
                             submitButton("Submit")
                             ),
                           mainPanel(
                             tabsetPanel(
                               tabPanel("Summary"),
                               tabPanel("Table"),
                               tabPanel("Plots")
                             )
                           )
                         )
                         )
                )
    )
)

server <- function(input, output, session) {
  # TODO
}


shinyApp(ui=ui, server=server)