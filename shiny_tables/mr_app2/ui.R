# UI
ui <- fluidPage(
  
  # Application title
  titlePanel(
    h1("Proteomic pan-MR with cis-coloc associations for Open Targets Genetics (bxy_pval < 0.0005 or coloc_h4_h3 > 1)", align = "center", style = "font-size: 30px")),
  
  DT::dataTableOutput("table")
  #  reactableOutput("table")
  
)