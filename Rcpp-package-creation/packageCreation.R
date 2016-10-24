install.packages("curl")
install.packages("crayon")

install.packages("devtools")
library("devtools")

devtools::install_github("klutometis/roxygen")
library(roxygen2)

setwd("C:/Users/838035/Dropbox/PhD/src/testing")
create("JMbayesAncillary")

# Now add all files in the R folder of the package


#Now we generate documentation
setwd("C:/Users/838035/Dropbox/PhD/src/testing/JMbayesAncillary")
document()

setwd("..")
install("cats")