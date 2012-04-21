plotK <- function(K,n=15,L=TRUE,persp=FALSE,legend=TRUE,...)
{
  if (isTRUE(L))
    {
      k <- K$Khat-K$Ktheo
      if (isTRUE(K$infectious))	
        titl <- expression(hat(K)[ST] * group("(",list(u,v),")") - pi*u^2*v)
      else
        titl <- expression(hat(K)[ST] * group("(",list(u,v),")") - 2*pi*u^2*v)
    }
  else
    {
      k <- K$Khat
      titl <- expression(hat(K)[I] * group("(",list(u,v),")") )
    }


  colo <- colorRampPalette(c("red", "white", "blue"))
  M <- max(abs(range(k)))
  M <- pretty(c(-M,M),n=n)
  n <- length(M)
  COL <- colo(n)
  if (isTRUE(persp))
    {
      mask <- matrix(0,ncol=length(K$times),nrow=length(K$dist))
      for(i in 1:length(K$dist)){ for(j in 1:length(K$times)){mask[i,j] <- COL[.iplace(x=k[i,j],X=M,xinc=diff(range(M))/(n-1))]}}
      COL <- mask[1:(length(K$dist)-1),1:(length(K$times)-1)]
      
      if(isTRUE(legend))
        {
          par(cex.lab=2,cex.axis=1.5,font=2,lwd=1,mar=c(0,0,3,0))
          par(fig=c(0,0.825,0,1))
          persp(x=K$dist, y=K$times, z=k, xlab="u",ylab="v", zlab="",expand=1, col=COL, ...)
          title(titl,cex.main=2)
          par(fig=c(0.825,1,0,1))
          mini <- .iplace(x=min(k,na.rm=T),X=M,xinc=diff(range(k))/(n-1))
          maxi <- .iplace(x=max(k,na.rm=T),X=M,xinc=diff(range(k))/(n-1))
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
          par(cex.lab=1.5,cex.axis=1.5,font=2,lwd=1,plt=c(0,1,0,1),mar=c(0.5,0.5,2.5,0.5),las=1)
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
          mini <- .iplace(x=min(k,na.rm=T),X=M,xinc=diff(range(k))/(n-1))
          maxi <- .iplace(x=max(k,na.rm=T),X=M,xinc=diff(range(k))/(n-1))
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
}
  
