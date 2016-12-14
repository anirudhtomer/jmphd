library(ggplot2)

##############################################
# Load the data set
##############################################
amctx = read.csv2(file.choose(), header = T)
#The zis number has one to one correspondance with amctx
amctx$zis = NULL

##############################################
# Check if any column has missing data: yuhoooooo none of them
##############################################
apply(amctx, MARGIN = 2, FUN= function(x){any(is.na(x))})

##############################################
# Data type cleaning before conversion to long
##############################################
amctx$amctx = as.factor(amctx$amctx)
amctx$gl_loss = factor(amctx$gl_loss, labels = c("no", "yes"))
amctx$gl_death = factor(amctx$gl_death, labels = c("no", "yes"))
amctx$gl_failure = factor(amctx$gl_failure, labels = c("no", "yes"))

##############################################
# PCR is the only output of interest
##############################################
amctx_pcr = amctx[amctx$measure=="pcr",]
