library(DT)
library(dplyr)

server <- function(input, output) { # adapted from here https://stackoverflow.com/questions/51425442/parent-child-rows-in-r-shiny-package
  
#   # load data (need to do this just once every time the core datasets are updated)
#   
#   df <- readRDS("~/mr_prot_snp_filtered_dataset_v1_v2.rds")
#   dfsnp <- df %>% select(exp_out_gsmr_coloc, snp, cis_trans, bzx, bzx_se, bzx_pval, bzy, bzy_se, bzy_pval, fpred_max_label_index, fpred_max_label_tag, pav_cis, hgnc_v2g)
# 
#   dfmr <- readRDS("~/mr_prot_filtered_dataset_v1_v2.rds")
# 
#   ## merge the row details
# 
#   expout <- unique(dfmr$exp_out_gsmr_coloc)
# 
#   dflst <- parallel::mclapply(seq(expout), function (x) {
# 
#     # print(expout[x])
# 
#     dfsnp1 <- dfsnp[which(dfsnp$exp_out_gsmr_coloc==expout[x]),]
#     dfmr1 <- dfmr[which(dfmr$exp_out_gsmr_coloc==expout[x]),]
# 
#     # print("success")
# 
#     subdats <- lapply(list(dfsnp1), purrr::transpose)
# 
#     Dat <- cbind(" " = "&oplus;", dfmr1, "_details" = I(subdats))
# 
#   }, mc.cores = parallel::detectCores())
# 
# 
# Dat <- data.table::rbindlist(dflst)
# 
# saveRDS(Dat, "~/mr_app2/mrdf_nested_for_app.rds")
#   
  

  Dat <- readRDS("mrdf_nested_for_app.rds") %>% select(1, "Data", "exposure", "protein_trait", "hgnc_protein", "ensid", "outcome", "outcome_trait", "outcome_trait_efo", "nsnp", "n_cases", "n_initial", "cis_trans_mr", "coloc", "varid_left", "coloc_h4", "coloc_h4_h3",  "bxy", "bxy_se", "bxy_pval", "bxy_pval_mantissa", "bxy_pval_exponent", "mr_egger", "mr_egger_se", "mr_egger_p", "mr_egger_p_mantissa", "mr_egger_p_exponent", "mr_wm", "mr_wm_se", "mr_wm_p", "mr_wm_p_mantissa", "mr_wm_p_exponent", "pav_cismr", "fpred_max_label_index_mr", "fpred_max_label_tag_mr", "v2g_mr", "hgnc_v2g_mr", "_details")
  
Dat <- Dat[with(Dat, order(bxy_pval_exponent, bxy_pval_mantissa)),]
  
  ## the callback
  callback = JS(
    "table.column(1).nodes().to$().css({cursor: 'pointer'});",
    "",
    "// make the table header of the nested table",
    "var format = function(d, childId){",
    "  if(d != null){",
    "    var html = ", 
    "      '<table class=\"display compact hover\" id=\"' + childId + '\"><thead><tr>';",
    "    for (var key in d[d.length-1][0]) {",
    "      html += '<th>' + key + '</th>';",
    "    }",
    "    html += '</tr></thead></table>'",
    "    return html;",
    "  } else {",
    "    return '';",
    "  }",
    "};",
    "",
    "// row callback to style the rows of the child tables",
    "var rowCallback = function(row, dat, displayNum, index){",
    "  if($(row).hasClass('odd')){",
    "    $(row).css('background-color', 'papayawhip');",
    "    $(row).hover(function(){",
    "      $(this).css('background-color', '#E6FF99');",
    "    }, function() {",
    "      $(this).css('background-color', 'papayawhip');",
    "    });",
    "  } else {",
    "    $(row).css('background-color', 'lemonchiffon');",
    "    $(row).hover(function(){",
    "      $(this).css('background-color', '#DDFF75');",
    "    }, function() {",
    "      $(this).css('background-color', 'lemonchiffon');",
    "    });",
    "  }",
    "};",
    "",
    "// header callback to style the header of the child tables",
    "var headerCallback = function(thead, data, start, end, display){",
    "  $('th', thead).css({",
    "    'border-top': '3px solid indigo',", 
    "    'color': 'indigo',",
    "    'background-color': '#fadadd'",
    "  });",
    "};",
    "",
    "// make the datatable",
    "var format_datatable = function(d, childId){",
    "  var dataset = [];",
    "  var n = d.length - 1;",
    "  for(var i = 0; i < d[n].length; i++){",
    "    var datarow = $.map(d[n][i], function (value, index) {",
    "      return [value];",
    "    });",
    "    dataset.push(datarow);",
    "  }",
    "  var id = 'table#' + childId;",
    "  if (Object.keys(d[n][0]).indexOf('_details') === -1) {",
    "    var subtable = $(id).DataTable({",
    "                 'data': dataset,",
    "                 'autoWidth': true,",
    "                 'deferRender': true,",
    "                 'info': false,",
    "                 'lengthChange': false,",
    "                 'ordering': d[n].length > 1,",
    "                 'order': [],",
    "                 'paging': false,",
    "                 'scrollX': false,",
    "                 'scrollY': false,",
    "                 'searching': false,",
    "                 'sortClasses': false,",
    "                 'rowCallback': rowCallback,",
    "                 'headerCallback': headerCallback,",
    "                 'columnDefs': [{targets: '_all', className: 'dt-center'}]",
    "               });",
    "  } else {",
    "    var subtable = $(id).DataTable({",
    "            'data': dataset,",
    "            'autoWidth': true,",
    "            'deferRender': true,",
    "            'info': false,",
    "            'lengthChange': false,",
    "            'ordering': d[n].length > 1,",
    "            'order': [],",
    "            'paging': false,",
    "            'scrollX': false,",
    "            'scrollY': false,",
    "            'searching': false,",
    "            'sortClasses': false,",
    "            'rowCallback': rowCallback,",
    "            'headerCallback': headerCallback,",
    "            'columnDefs': [", 
    "              {targets: -1, visible: false},", 
    "              {targets: 0, orderable: false, className: 'details-control'},", 
    "              {targets: '_all', className: 'dt-center'}",
    "             ]",
    "          }).column(0).nodes().to$().css({cursor: 'pointer'});",
    "  }",
    "};",
    "",
    "// display the child table on click",
    "table.on('click', 'td.details-control', function(){",
    "  var tbl = $(this).closest('table'),",
    "      tblId = tbl.attr('id'),",
    "      td = $(this),",
    "      row = $(tbl).DataTable().row(td.closest('tr')),",
    "      rowIdx = row.index();",
    "  if(row.child.isShown()){",
    "    row.child.hide();",
    "    td.html('&oplus;');",
    "  } else {",
    "    var childId = tblId + '-child-' + rowIdx;",
    "    row.child(format(row.data(), childId)).show();",
    "    td.html('&CircleMinus;');",
    "    format_datatable(row.data(), childId);",
    "  }",
    "});")
  
  ## datatable
  output$table <- DT::renderDataTable({
    datatable(Dat, 
              callback = callback, 
              escape = -2,
              filter = list(position = 'top'),
              
              options = list(
              columnDefs = list(
                list(visible = FALSE, targets = ncol(Dat)),
                list(orderable = FALSE, className = 'details-control', targets = 1),
                list(className = "dt-center", targets = "_all"),
                list(targets = '_all', className = 'dt-center')
              ),
              pageLength = 10,  ## number of rows to output for each page
              scrollX = TRUE,   ## enable scrolling on X axis
              scrollY = TRUE,   ## enable scrolling on Y axis
              autoWidth = TRUE, ## use smart column width handling
              server = TRUE   ## use client-side processing
            ))
    
  })
}