library('ggplot2')
library('iNEXT')
library("data.table")


data_r = fread("prepare_for_inext_input.txt",sep = '\t', header=TRUE)
otutable=list(as.matrix(data_r))
otu <- iNEXT(otutable, q=0, datatype="incidence_raw",endpoint=1200000,knots=2000)
otu$iNextEst
