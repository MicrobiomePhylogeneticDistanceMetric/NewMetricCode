idata <-read.csv('Wavelet Transforms in R/Raw Data/stool-rawdata.csv',header = TRUE, stringsAsFactors = FALSE)
idata <- idata[,-1]
cdata <- as.matrix(idata)

three_band <- function(X){
  n=dim(X)[1]
  m=dim(X)[2]
  rnames <- rownames(X)
  cnames <- colnames(X)
  h0=c(0.33838609728386,0.53083618701374,0.72328627674361,0.23896417190576,0.04651408217589,-0.14593600755399)
  h1=c(-0.11737701613483,0.54433105395181,-0.01870574735313,-0.69911956479289,-0.13608276348796,0.42695403781698)
  h2=c(0.40363686892892,-0.62853936105471,0.46060475252131,-0.40363686892892,-0.07856742013185,0.24650202866523)
  nn=3*ceiling(n/3)
  
  vv1=matrix(0,nn,m)
  rownames(vv1) <- rnames
  colnames(vv1) <- cnames
  ww1=matrix(0,nn,m)
  rownames(ww1) <- rnames
  colnames(ww1) <- cnames
  ww2=matrix(0,nn,m)
  rownames(ww2) <- rnames
  colnames(ww2) <- cnames
  for (k in 1:m){
    ss=X[,k]
    if(n%%3){ss=c(ss,rep(0,3-n%%3))}
    ss=c(ss,ss[1:3])
    v1=rep(0,nn/3)
    w1=rep(0,nn/3)
    w2=rep(0,nn/3)
    j=1
    for (i in 1:(nn/3)){
      v1[i]=sum(h0*ss[j:(j+5)])
      w1[i]=sum(h1*ss[j:(j+5)])
      w2[i]=sum(h2*ss[j:(j+5)])
      j=j+3
    }
    v1=c(v1[nn/3],v1)
    w1=c(w1[nn/3],w1)
    w2=c(w2[nn/3],w2)
    
    for (i in 1:(nn/3)){
      vv1[3*i-2,k]=h0[4]*v1[i]+h0[1]*v1[i+1]
      vv1[3*i-1,k]=h0[5]*v1[i]+h0[2]*v1[i+1]
      vv1[3*i,k]=h0[6]*v1[i]+h0[3]*v1[i+1]
      ww1[3*i-2,k]=h1[4]*w1[i]+h1[1]*w1[i+1]
      ww1[3*i-1,k]=h1[5]*w1[i]+h1[2]*w1[i+1]
      ww1[3*i,k]=h1[6]*w1[i]+h1[3]*w1[i+1]
      ww2[3*i-2,k]=h2[4]*w2[i]+h2[1]*w2[i+1]
      ww2[3*i-1,k]=h2[5]*w2[i]+h2[2]*w2[i+1]
      ww2[3*i,k]=h2[6]*w2[i]+h2[3]*w2[i+1]
    }
  }
  vv1=vv1[1:n,]
  ww1=ww1[1:n,]
  ww2=ww2[1:n,]
  write.csv(vv1, file = "Oral3WLow", row.names = FALSE)
  write.csv(ww1, file = "Oral3WHigh1", row.names = FALSE)
  write.csv(ww2, file = "Oral3WHigh2", row.names = FALSE)
}

three_band(cdata)