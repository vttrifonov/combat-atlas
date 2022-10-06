#install.packages('renv')

#renv::init(bare=TRUE)

renv::install(c(
  'data.table',
  'magrittr',
  'markdown',
  'rmarkdown',
  'Matrix',
  'bioc::limma',
  'reticulate',
  'languageserver',
  'httpgd',
  'ggplot2',
  'glue',
  'reshape2',
  'gridExtra',
  'bioc::fgsea',
  'DT',
  'shiny',
  'bioc::GSVA'
))
