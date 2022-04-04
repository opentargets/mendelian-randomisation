
library(dplyr)

# Server
server <- function(input, output) {
  
  df <- readRDS("5_combined_filtered_annotated_metab_mr_file.rds")
  
  df <- df %>% select(Data, trait, outcome, n_initial, n_cases, outcome_trait, outcome_trait_category, outcome_trait_efo2,
                      cis_trans2, cis_ensid, cis_gene, 
                      nsnp, 
                      bxy, bxy_lci, bxy_uci, or_bxy, or_lci_bxy, or_uci_bxy, bxy_pval, bxy_pval_mantissa, bxy_pval_exponent,
                      or_mr_egger_int, or_lci_mr_egger_int, or_uci_mr_egger_int, mr_egger_int_p, mr_egger_int_p_mantissa, mr_egger_int_p_exponent,
                      or_mr_wm, or_lci_mr_wm, or_uci_mr_wm, mr_wm_p, mr_wm_p_mantissa, mr_wm_p_exponent)
  
  output$table <- DT::renderDataTable({
    DT::datatable(df,
                  filter = list(position = 'top'), ## turn on per-column search box
                  options = list(pageLength = 10,  ## number of rows to output for each page
                                 scrollX = TRUE,   ## enable scrolling on X axis
                                 scrollY = TRUE,   ## enable scrolling on Y axis
                                 autoWidth = TRUE, ## use smart column width handling
                                 server = TRUE,   ## use client-side processing
                                 columnDefs = list(list(targets = '_all', className = 'dt-center'))
                  )
    ) %>%
      DT::formatSignif(
        which(sapply(df, function(x) class(x) == "numeric") == TRUE), 
        digits = 5
      )
  })
  
  
  # output$table <- renderReactable({
  #   reactable(df, minRows = 10, searchable = TRUE, bordered = TRUE,
  # highlight = TRUE,  defaultColDef = reactable::colDef(
  #   align = "center",
  #   minWidth = 170,
  #   headerStyle = list(background = "#c3c3c9")))
  # })
}