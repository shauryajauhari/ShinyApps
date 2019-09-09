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

# ui <- fluidPage( h1("Hello World!"), ## Using h1/ heading tag from HTML
#                  sidebarLayout(sidebarPanel(),
#                                mainPanel())
#                  ) 
# server <- function(input, output, session) {
# 
# }
# 
# shinyApp(ui, server)

## Program4

# ui <- fluidPage( 
#   h1("Hello World!"), ## Using h1/ heading tag from HTML
#                  sidebarLayout(
#                    sidebarPanel(
#                      selectInput("dataset", "Choose a dataset:",
#                                  choices = ls ("package:datasets"),
#                                  selected = "pressure")
#                    ),
#                    mainPanel()
#                    )
# ) 
# server <- function(input, output, session) {
#   
# }
# 
# shinyApp(ui, server)


ui <- fluidPage(
  titlePanel(title = "Welcome to Shaurya's Blog."),
  h4("This application is for getting familiar with the basic concepts of 
     Shiny. To begin with, let's explore some datasets bundled with R."),
  sidebarLayout(
    sidebarPanel(),
    mainPanel()
    )
) 
server <- function(input, output, session) {
  plot()
  
}

shinyApp(ui, server)


