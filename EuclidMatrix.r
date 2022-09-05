idata <- read.csv("oral-rawdata.csv",header = TRUE, stringsAsFactors = FALSE)
ldata <- idata[,-1]
jdata <-as.matrix(ldata)
kdata <- log(1000000* (jdata + (1-min(jdata))))
mdata<- t(kdata)

euclidean <- function(a, b) sqrt(sum((a - b)^2))
eh <- function(X) {
  d <- matrix(0, nrow = nrow(X), ncol = nrow(X))
  for ( i in 1:(nrow(X) - 1) ) {
    for ( j in (i + 1):nrow(X) ) {
      d[i,j] <- d[j,i] <- euclidean(X[i,], X[j,])
    }
  }
  d
}

euclid<-eh(mdata)
write.csv(euclid, file = "RawDataEuclid" , row.names = FALSE)
