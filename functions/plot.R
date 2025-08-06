# Polygons/Boxes functions
polygons = function(x, XLAB='', YLAB='', MAIN='', DATE=FALSE, DATES=NULL, YLIM=FALSE, YLIMS=NULL ) {
  
  nTimes = dim(x)[2]
  if (DATE) {
    days = DATES
  } else {
    days = colnames(x)
  }
  low_counts = apply( x, 2, quantile, 0.025, na.rm=T)
  upp_counts = apply( x, 2, quantile, 0.975, na.rm=T)
  pnt_counts = apply( x, 2, mean, na.rm=T)
  MN = min( low_counts, na.rm=T )
  MX = max( upp_counts,na.rm=T )
  if (YLIM) {
    plot(1:nTimes, NA*pnt_counts, type='l', main=MAIN, xlab=XLAB, ylab=YLAB, ylim=YLIMS,xaxt='n',cex.main=1.5,cex.axis=1.25,cex.lab=1.25)
  } else {
    plot(1:nTimes, NA*pnt_counts, type='l', main=MAIN, xlab=XLAB, ylab=YLAB, ylim=c(MN,MX),xaxt='n',cex.main=1.5,cex.axis=1.25,cex.lab=1.25)
  }
  
  
  
  
  tmp1 = c( 1:nTimes, rev(1:nTimes) )
  tmp2 = c( low_counts, rev(upp_counts) )
  polygon( tmp1, tmp2, col=rgb(211/255,211/255,211/255,alpha=0.33), border=NA)
  lines(pnt_counts,col='darkgray',lwd=2)
  plot_at = round(seq(1,nTimes,length.out=4))
  axis(1,plot_at,days[plot_at],cex.axis=1.25)
  
}
boxes    = function(x, XLAB='', YLAB='', MAIN='', plus=0.475, outliers=TRUE, custom_max=list(FALSE,0,0) ) {
  
  if (plus>0.475) {
    plus = 0.475
  }
  if (plus<0) {
    plus = 0.475
  }
  
  if ( outliers==TRUE ) {
    MN = min(x,na.rm=T)
    MX = max(x,na.rm=T)
  } else {
    MN = min( apply(x,2,quantile,0.025,na.rm=T) )
    MX = max( apply(x,2,quantile,0.975,na.rm=T) )
  }
  
  if (custom_max[[1]]) {
    MN = custom_max[[2]]
    MX = custom_max[[3]]
  }
  
  nCol = dim(x)[2]
  XLABELS = colnames(x)
  plot_at = round( seq(1,nCol,length.out=min(nCol,6)) )
  plot(NA, NA, ylim=c(MN,MX),xlim=c(1-plus,nCol+plus), xlab=XLAB, ylab=YLAB, main=MAIN, xaxt='n', cex.main=1.5,cex.axis=1.25,cex.lab=1.25)
  axis( 1, (1:nCol)[plot_at], XLABELS[plot_at],cex.axis=1.25 )
  
  for (i in 1:nCol) {
    
    q025 = quantile( x[,i], 0.025, na.rm=T )
    q975 = quantile( x[,i], 0.975, na.rm=T )
    q25  = quantile( x[,i], 0.25, na.rm=T )
    q75  = quantile( x[,i], 0.75, na.rm=T )
    
    tmp1 = c( i-plus, i+plus, i+plus, i-plus )
    tmp2 = c( q25, q25, q75, q75 )
    polygon(tmp1,tmp2,col='lightgray')
    
    segments( i-plus, mean(x[,i],na.rm=T), i+plus, mean(x[,i],na.rm=T) , lty=2, col=2, lwd=0.5)
    segments( i-plus, median(x[,i],na.rm=T), i+plus, median(x[,i],na.rm=T) , lty=1, lwd=1.0)
    segments( i-.5*plus, q025, i+.5*plus, q025, lty=1)
    segments( i-.5*plus, q975, i+.5*plus, q975, lty=1)
    segments( i, q75, i, q975, lty=2 )
    segments( i, q25, i, q025, lty=2 )
    
    if (outliers) {
      tmp1 = x[x[,i]<q025,i]
      tmp2 = x[x[,i]>q975,i]
      tmp  = c(tmp1,tmp2)
      if (length(tmp)>0) {
        points( rep(i,length(tmp)),tmp, pch=8 ,cex=0.2)
      }
    }
    
  }
}
