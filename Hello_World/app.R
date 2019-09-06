# Hello World
## Program1

# ui <- "Hello World!"
# 
# server <- function(input, output, session) {
#   
# }
# 
# shinyApp(ui, server)


## Program2

# ui <- fluidPage( h1("Hello World!")) ## Using h1/ heading tag from HTML
# server <- function(input, output, session) {
#   
# }
# 
# shinyApp(ui, server)

## Program3

ui <- fluidPage( h1("Hello World!"), ## Using h1/ heading tag from HTML
                 sidebarLayout(sidebarPanel(),
                               mainPanel())
                 ) 
server <- function(input, output, session) {

}

shinyApp(ui, server)

## Program4

ui <- fluidPage( 
  h1("Hello World!"), ## Using h1/ heading tag from HTML
                 sidebarLayout(
                   sidebarPanel(
                     selectInput("dataset", "Choose a dataset:",
                                 choices = ls ("package:datasets"),
                                 selected = "pressure")
                   ),
                   mainPanel()
                   )
) 
server <- function(input, output, session) {
  
}

shinyApp(ui, server)

