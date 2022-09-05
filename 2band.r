idata <-read.csv('Wavelet Transforms in R/Raw Data/stool-rawdata.csv',header = TRUE, stringsAsFactors = FALSE)
idata <- idata[,-1]
cdata <- as.matrix(idata)

two_band <- function(X){
  n=dim(X)[1]
  m=dim(X)[2]
  h0=c(0.6830127,1.1830127,0.3169873,-0.1830127)
  h1=c(-0.1830127, -0.3169873, 1.1830127,-0.6830127)
  nn=2*ceiling(n/2)
  
  vv1=matrix(0,nn,m)
  ww1=matrix(0,nn,m)
  for (k in 1:m){
    ss=X[,k]
    if(n%%2){ss=c(ss,rep(0,2-n%%2))}
    ss=c(ss,ss[1:2])
    v1=rep(0,nn/2)
    w1=rep(0,nn/2)
    j=1
    for (i in 1:(nn/2)){
      v1[i]=sum(h0*ss[j:(j+3)])
      w1[i]=sum(h1*ss[j:(j+3)])
      j=j+2
    }
    v1=c(v1[nn/2],v1)
    w1=c(w1[nn/2],w1)
    
    for (i in 1:(nn/2)){
      vv1[2*i-1,k]=h0[3]*v1[i]+h0[1]*v1[i+1]
      vv1[2*i,k]=h0[4]*v1[i]+h0[2]*v1[i+1]
      ww1[2*i-1,k]=h1[3]*w1[i]+h1[1]*w1[i+1]
      ww1[2*i,k]=h1[4]*w1[i]+h1[2]*w1[i+1]
    }
  }
  vv1=vv1[1:n,]
  ww1=ww1[1:n,]
  write.csv(vv1, file = "OralW2L.csv", row.names = FALSE)
  write.csv(ww1, file = "OralW2H.csv", row.names = FALSE)
}

two_band(cdata)
