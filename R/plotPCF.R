plotPCF <- function(PCF,n=15,persp=FALSE,legend=TRUE,...)
{
devl=dev.list()
if(is.null(devl))
{
old.par <- par(no.readonly = TRUE)    
on.exit(par(old.par))
}
else
{
dev.off(devl[length(devl)])
X11()
old.par <- par(no.readonly = TRUE)    
on.exit(par(old.par))
}

 k=PCF$pcf
 K=PCF

 titl <- expression(hat(g)* group("(",list(u,v),")") )


  colo <- colorRampPalette(c("red", "white", "blue"))
  M <- max(abs(range(k)))
  M <- pretty(c(-M,M),n=n)
  n <- length(M)
  COL <- colo(n)
  if (isTRUE(persp))
    {
      mask <- matrix(0,ncol=length(K$times),nrow=length(K$dist))
      for(i in 1:length(K$dist)){ for(j in 1:length(K$times)){mask[i,j] <- COL[findInterval(x=k[i,j],vec=M)]}}
      COL <- mask[1:(length(K$dist)-1),1:(length(K$times)-1)]
      
      if(isTRUE(legend))
        {
          par(cex.lab=2,cex.axis=1.5,font=2,lwd=1,mar=c(0,0,3,0))
          par(fig=c(0,0.825,0,1))
          persp(x=K$dist, y=K$times, z=k, xlab="u",ylab="v", zlab="",expand=1, col=COL, ...)
          title(titl,cex.main=2)
          par(fig=c(0.825,1,0,1))
          mini <- findInterval(x=min(k,na.rm=TRUE),vec=M)
          maxi <- findInterval(x=max(k,na.rm=TRUE),vec=M)
          legend("right",fill=colo(n)[maxi:mini],legend=M[maxi:mini],horiz=F,bty="n")
        }
      else
        {
          par(cex.lab=2,cex.axis=1.5,font=2,lwd=1)
          persp(x=K$dist, y=K$times, z=k, xlab="u",ylab="v", zlab="", expand=1, col=COL, ...)
          title(titl,cex.main=2)
        }
    }
  else
    {
      if(isTRUE(legend))
        {
          par(cex.lab=1.5,cex.axis=1.5,font=2,lwd=1,plt=c(0,1,0,1),mar=c(0.5,0.5,2.5,.5),las=1)
          par(fig=c(0.1,0.825,0.1,1))
          contour(K$dist, K$times, k, labcex=1.5,levels=M,drawlabels=F,col=colo(n),zlim=range(M),axes=F)
          box(lwd=2)
          at <- axTicks(1)
          axis(1,at=at[1:(length(at)-1)],labels=at[1:(length(at)-1)])
          axis(1,at=at[length(at)],labels="u",cex.axis=2)
          at <- axTicks(2)
          axis(2,at=at[1:(length(at)-1)],labels=at[1:(length(at)-1)])
          axis(2,at=at[length(at)],labels="v",cex.axis=2)
          title(titl,cex.main=2)
	    par(fig=c(0,1,0.1,1))
          mini <- findInterval(x=min(k,na.rm=TRUE),vec=M)
          maxi <- findInterval(x=max(k,na.rm=TRUE),vec=M)
          legend("right",fill=colo(n)[maxi:mini],legend=M[maxi:mini],horiz=F,bty="n")
        }
      else
        {
          par(cex.lab=2,cex.axis=1.5,font=2,lwd=2,las=1)
          contour(K$dist, K$times, k, labcex=1.5,levels=M,drawlabels=T,col=colo(n),zlim=range(M),axes=F)
          box(lwd=2)
          at <- axTicks(1)
          axis(1,at=at[1:(length(at)-1)],labels=at[1:(length(at)-1)])
          axis(1,at=at[length(at)],labels="u",cex.axis=2)
          at <- axTicks(2)
          axis(2,at=at[1:(length(at)-1)],labels=at[1:(length(at)-1)])
          axis(2,at=at[length(at)],labels="v",cex.axis=2)
          title(titl,cex.main=2)
        }
    }
par(old.par)
}
  
