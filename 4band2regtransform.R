idata <-read.csv('Wavelet Transforms in R/Raw Data/stool-rawdata.csv',header = TRUE, stringsAsFactors = FALSE)
idata <- idata[,-1]
cdata <- as.matrix(idata)

four_band2reg <- function(X){
  n=dim(X)[1]
  m=dim(X)[2]
  h0=c(-0.067371764,0.094195111,0.40580489,0.567371764,0.567371764,0.40580489,0.094195111,-0.067371764)
  h1=c(-0.094195111,0.067371764,0.567371764 ,0.40580489,-0.40580489,-0.567371764,-0.067371764,0.094195111)
  h2=c(-0.094195111,-0.067371764,0.567371764,-0.40580489,-0.40580489,0.567371764,-0.067371764,-0.094195111)
  h3=c(-0.067371764,-0.094195111,0.40580489,-0.567371764,0.567371764,-0.40580489,0.094195111,0.067371764)
  nn=4*ceiling(n/4)
  
  vv1=matrix(0,nn,m)
  ww1=matrix(0,nn,m)
  ww2=matrix(0,nn,m)
  ww3=matrix(0,nn,m)
  for (k in 1:m){
    ss=X[,k]
    if(n%%4){ss=c(ss,rep(0,4-n%%4))}
    ss=c(ss,ss[1:4])
    v1=rep(0,nn/4)
    w1=rep(0,nn/4)
    w2=rep(0,nn/4)
    w3=rep(0,nn/4)
    j=1
    for (i in 1:(nn/4)){
      v1[i]=sum(h0*ss[j:(j+6)])
      w1[i]=sum(h1*ss[j:(j+6)])
      w2[i]=sum(h2*ss[j:(j+6)])
      w3[i]=sum(h3*ss[j:(j+6)])
      j=j+4
    }
    v1=c(v1[nn/4],v1)
    w1=c(w1[nn/4],w1)
    w2=c(w2[nn/4],w2)
    w3=c(w3[nn/4],w3)
    
    for (i in 1:(nn/4)){
      vv1[4*i-3,k]=h0[3]*v1[i]+h0[1]*v1[i+1]
      vv1[4*i-2,k]=h0[4]*v1[i]+h0[2]*v1[i+1]
      vv1[4*i-1,k]=h0[5]*v1[i]+h0[3]*v1[i+1]
      vv1[4*i,k]=h0[6]*v1[i]+h0[4]*v1[i+1]
      ww1[4*i-3,k]=h0[3]*w1[i]+h0[1]*w1[i+1]
      ww1[4*i-2,k]=h0[4]*w1[i]+h0[2]*w1[i+1]
      ww1[4*i-1,k]=h0[5]*w1[i]+h0[3]*w1[i+1]
      ww1[4*i,k]=h0[6]*w1[i]+h0[4]*w1[i+1]
      ww2[4*i-3,k]=h0[3]*w2[i]+h0[1]*w2[i+1]
      ww2[4*i-2,k]=h0[4]*w2[i]+h0[2]*w2[i+1]
      ww2[4*i-1,k]=h0[5]*w2[i]+h0[3]*w2[i+1]
      ww2[4*i,k]=h0[6]*w2[i]+h0[4]*w2[i+1]
      ww3[4*i-3,k]=h0[3]*w3[i]+h0[1]*w3[i+1]
      ww3[4*i-2,k]=h0[4]*w3[i]+h0[2]*w3[i+1]
      ww3[4*i-1,k]=h0[5]*w3[i]+h0[3]*w3[i+1]
      ww3[4*i,k]=h0[6]*w3[i]+h0[4]*w3[i+1]
    }
  }
  vv1=vv1[1:n,]
  ww1=ww1[1:n,]
  ww2=ww2[1:n,]
  ww3=ww3[1:n,]
  write.csv(vv1, file = "Oral4W2RLow.csv", row.names = FALSE)
  write.csv(ww1, file = "Oral4W2RHigh1.csv", row.names = FALSE)
  write.csv(ww2, file = "Oral4W2RHigh2.csv", row.names = FALSE)
  write.csv(ww3, file = "Oral4W2RHigh3.csv", row.names = FALSE)
}

four_band2reg(cdata)