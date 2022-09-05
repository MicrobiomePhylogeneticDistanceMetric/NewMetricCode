idata <-read.csv("Wavelet Transforms in R/Raw Data/ORAL RAW DATA WAVELETS/oral-4W4H3.csv",header = TRUE, stringsAsFactors = FALSE)
ldata <- idata[,-1]
jdata <-as.matrix(ldata)
kdata <- log(1000000* (jdata + (1-min(jdata))))
mdata<- t(kdata)

Bray_Curtis <- function(X){
  input <- data.matrix(X)
  bc_result <- bcdist(input, rmzero = FALSE)
  bc.matrix <- data.matrix(bc_result)
  write.csv(bc.matrix, file = "BCDist4B4RHigh3.csv", row.names = FALSE)
}

Bray_Curtis(mdata)