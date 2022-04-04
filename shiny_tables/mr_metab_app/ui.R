# UI
ui <- fluidPage(
  
  # Application title
  titlePanel(
    h1("Metabolomic MR associations (V1) for Open Targets Genetics (5% FDR)", align = "center", style = "font-size: 30px")),
  
  DT::dataTableOutput("table")
  #  reactableOutput("table")
  
)