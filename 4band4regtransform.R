idata <-read.csv('Wavelet Transforms in R/Raw Data/stool-rawdata.csv',header = TRUE, stringsAsFactors = FALSE)
idata <- idata[,-1]
cdata <- as.matrix(idata)

four_band4reg <- function(X){
  n=dim(X)[1]
  m=dim(X)[2]
  rnames <- rownames(X)
  cnames <- colnames(X)
  h0=c(0.0857130200,0.1931394393,0.3491805097,0.5616494215,0.4955029828,0.4145647737,0.2190308939,-0.1145361261,-0.0952930728,-0.1306948909,-0.0827496793,0.0719795354,0.0140770701,0.0229906779,0.0145382757,-0.0190928308)
  h1=c(-0.1045086525,0.1183282069,-0.1011065044,-0.0115563891,0.6005913823,-0.2550401616,-0.4264277361,-0.0827398180,0.0722022649,0.2684936992,0.1691549718,-0.4437039320,0.0849964877,0.1388163056,0.0877812188,-0.1152813433)
  h2=c(0.2560950163,-0.2048089157,-0.2503433230,-0.2484277272,0.4477496752,0.0010274000,-0.0621881917,0.5562313118,-0.2245618041,-0.3300536827,-0.2088643503,0.2202951830,0.0207171125,0.0338351983,0.0213958651,-0.0280987676)
  h3=c(0.1839986022,-0.6622893130,0.6880085746,-0.1379502447,0.0446493766,-0.0823301969,-0.0923899104,-0.0233349758,0.0290655661,0.0702950474,0.0443561794,-0.0918374833,0.0128845052,0.0210429802,0.0133066389,-0.0174753464)
  nn=4*ceiling(n/4)
  
  vv1=matrix(0,nn,m)
  rownames(vv1) <- rnames
  colnames(vv1) <- cnames
  ww1=matrix(0,nn,m)
  rownames(ww1) <- rnames
  colnames(ww1) <- cnames
  ww2=matrix(0,nn,m)
  rownames(ww2) <- rnames
  colnames(ww2) <- cnames
  ww3=matrix(0,nn,m)
  rownames(ww3) <- rnames
  colnames(ww3) <- cnames
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
  write.csv(vv1, file = "Oral4W4RLow", row.names = FALSE)
  write.csv(ww1, file = "Oral4W4RHigh1", row.names = FALSE)
  write.csv(ww2, file = "Oral4W4RHigh2", row.names = FALSE)
  write.csv(ww3, file = "Oral4W4RHigh3", row.names = FALSE)
}

four_band4reg(cdata)