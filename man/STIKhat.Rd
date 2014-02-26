\name{STIKhat} \alias{STIKhat}
\title{Estimation of the Space-Time Inhomogeneous K-function}

\description{Compute an estimate of the Space-Time
Inhomogeneous K-function.}

\usage{STIKhat(xyt, s.region, t.region, dist, times, lambda, 
correction="isotropic", infectious=FALSE) }

\arguments{
  \item{xyt}{coordinates and times (x,y,t) of the point pattern.}
    \item{s.region}{two-column matrix specifying polygonal region containing
  all data locations.
   If \code{s.region} is missing, the bounding box of \code{xyt[,1:2]} is considered.}
  \item{t.region}{vector containing the minimum and maximum values of
  the time interval.
   If \code{t.region} is missing, the range of \code{xyt[,3]} is considered.}
  \item{dist}{vector of distances u at which \eqn{K(u,v)}{K(u,v)} is computed.
   If missing, the maximum of \code{dist} is given by min(Sx,Sy)/4, where Sx and Sy
   represent the maximum width and height of \code{s.region}.
  
  }
   \item{times}{vector of times v at which \eqn{K(u,v)}{K(u,v)} is computed.}
    \item{lambda}{vector of values of the space-time intensity function
  evaluated at the points (x,y,t) in SxT.
  If \code{lambda} is missing, the estimate of the space-time
K-function is computed as for the homogeneous case (Diggle et
al., 1995), i.e. considering \eqn{n/|S \times T|}{n/|SxT|} as
an estimate of the space-time intensity.}
\item{correction}{A character vector specifying the edge correction(s) to be applied 
among \code{"isotropic"}, \code{"border"}, \code{"modified.border"}, \code{"translate"}
 and \code{"none"} (see Details). The default is \code{"isotropic"}.}
\item{infectious}{logical value. If \code{TRUE}, only future events
  are considered and the isotropic edge correction method is used. See Details.}
 }

\details{Gabriel (2014) proposes the following unbiased estimator for the STIK-function, based
on data giving the locations of events \eqn{x_i: i=1,...n}{xi:
i=1...,n} on a spatio-temporal region \eqn{S \times
T}{SxT}, where S is an arbitrary polygon and T is a time interval:
\deqn{\widehat {K}(u,v) = \sum_{i=1}^{n}
  \sum_{j\neq i} \frac{1}{w_{ij}} \frac{1}{\lambda(x_i) \lambda(x_j)}{\bf 1}_{\lbrace \|s_i - s_j\| \leq u
  \ ; \ |t_i - t_j| \leq v \rbrace}}{K(u,v)=sum_{i=1,...,n} sum_{j != i} 1/wij 
  1/(lambda(x_i)lambda(x_j)) 1{uij <= u ; tj-ti <= v}}
where \eqn{\lambda(x_i)}{lambda(xi)} is the intensity at \eqn{x_i = (s_i,t_i)}{xi=(si,ti)} and
\eqn{w_{ij}}{wij} is an edge correction factor to deal with spatial-temporal edge effects. 
The edge correction methods implemented are:

\code{isotropic}: \eqn{w_{ij} = |S \times T| w_{ij}^{(t)} w_{ij}^{(s)}}{wij = |S x T| wij(s) wij(t)}, 
where the temporal edge correction factor \eqn{w_{ij}^{(t)} = 1}{wij(t)=1} if both ends of the interval 
of length \eqn{2 |t_i - t_j|}{2|ti-tj|} centred at \eqn{t_i}{ti} lie within T and \eqn{w_{ij}^{(t)}=1/2}{wij(t)=1/2} 
otherwise and \eqn{w_{ij}^{(s)}}{wij(s)} is the proportion of the circumference of a circle centred 
at the location \eqn{s_i}{si} with radius \eqn{\|s_i -s_j\|}{||si-sj||} lying in S (also called Ripley's
edge correction factor).

\code{border}: \eqn{w_{ij} = \frac{\sum_{j=1}^n {\bf 1}\lbrace d(s_j,S) > u \ ; \ d(t_j,T) >v\rbrace / 
\lambda(x_j)}{{\bf 1}_{\lbrace d(s_i,S) > u \ ; \ d(t_i,T) >v \rbrace}}}{wij=(sum_j 1{d(sj,S)>u ; d(tj,T)>v}/
lambda(xj)) / 1{d(si,S)>u ; d(ti,T)>v}}, where \eqn{d(s_i,S)}{d(si,S)} denotes the distance 
between \eqn{s_i}{si} and the boundary of S and \eqn{d(t_i,T)}{d(ti,T)} the distance between 
\eqn{t_i}{ti} and the boundary of T.

\code{modified.border}: \eqn{w_{ij} = \frac{|S_{\ominus u}|\times|T_{\ominus v}|}
{{\bf 1}_{\lbrace d(s_i,S) > u \ ; \ d(t_i,T) >v \rbrace}}}{wij = |S(-u) x T(-v)| / 1{d(si,S)>u ; d(ti,T)>v}}, where
\eqn{S_{\ominus u}}{S(-u)} and \eqn{T_{\ominus v}}{T(-v)} are the eroded spatial and temporal region 
respectively, obtained by trimming off a margin of width u and v from the border of the original region.

\code{translate}: \eqn{w_{ij} =|S \cap S_{s_i-s_j}| \times |T \cap T_{t_i-t_j}|}{|S intersect S(si-sj)
 x T intersect T(ti-tj)|}, where \eqn{S_{s_i-s_j}}{S(si-sj)} and \eqn{T_{t_i-t_j}}{T(ti-tj)} 
 are the translated spatial and temporal regions.

\code{none}: No edge correction is performed and \eqn{w_{ij}=|S \times T|}{|S x T|}.

If parameter \code{infectious = TRUE}, ony future 
events are considered and the estimator is, using an isotropic edge correction factor 
(Gabriel and Diggle, 2009): 
\deqn{\widehat K(u,v) =  \frac{1}{|S \times T|} \frac{n}{n_v} \sum_{i=1}^{n_v}
  \sum_{j=1; j > i}^{n_v} \frac{1}{w_{ij}} \frac{1}{\lambda(x_i) \lambda(x_j)}{\bf 1}_{\left\lbrace u_{ij} \leq u
  \right\rbrace} {\bf 1}_{\left\lbrace t_j - t_i \leq v \right\rbrace}}{K(u,v) = 1/|SxT| n/nv sum_{i=1,...,nv}
  sum_{j=1,...,nv; j > i} 1/wij 1/(lambda(x_i)lambda(x_j)) 1{uij <=
  u} 1{tj - ti <= v}}
In this equation, the points \eqn{x_i=(s_i,t_i)}{xi=(si,ti)}
 are ordered so that \eqn{t_i < t_{i+1}}{ti <
t(i+1)}, with ties due to round-off error broken by randomly
unrounding if necessary. To deal with temporal edge-effects,
for each v, \eqn{n_v}{nv} denotes the number of events for
which \eqn{t_i \leq T_1 -v}{ti <= T1-v}, with \eqn{T=[T_0,T_1]}{T=[T0,T1]}. 
To deal with spatial edge-effects, we use Ripley's method.

If \code{lambda} is missing in argument, \code{STIKhat}
computes an estimate of the space-time (homogeneous)
K-function: \deqn{\widehat{K}(u,v) = \frac{|S \times
T|}{n_v(n-1)} \sum_{i=1}^{n_v}
  \sum_{j=1; j > i}^{n_v} \frac{1}{w_{ij}} {\bf 1}_{\lbrace u_{ij} \leq u \rbrace} {\bf 1}_{ \lbrace t_j - t_i \leq v
  \rbrace}}{K(u,v) = |SxT|/(nv(n-1)) sum_{i=1,...,nv} sum_{j=1,...,nv; j>i} 1/wij 1{uij <=
  u} 1{tj - ti <= v}}
}

 \value{ A list containing:

 \item{Khat}{ndist x ntimes matrix containing values of \eqn{\hat
 K_{ST}(u,v)}{K(u,v)}.}
\item{Ktheo}{ndist x ntimes matrix containing theoretical
values for a Poisson process; \eqn{\pi u^2 v}{pi u^2 v} for K
and \eqn{2 \pi u^2 v}{2 pi u^2 v}) for K^*.}

 \item{dist, times, infectious}{parameters passed in argument.}
 
  \item{correction}{the name(s) of the edge correction method(s) passed in argument.}
 }

\references{
Gabriel E. (2014) Estimating second-order characteristics of 
inhomogeneous spatio-temporal point processes: influence of 
edge correction methods and intensity estimates.
Methodology and computing in Applied Probabillity, 16(1).

Gabriel E., Diggle P. (2009) Second-order analysis
of inhomogeneous spatio-temporal point process data. Statistica
Neerlandica, 63, 43--51.

Gabriel E., Rowlingson B., Diggle P. (2013)
stpp: an R package for plotting, simulating and analyzing 
Spatio-Temporal Point Patterns. Journal of Statistical Software, 
53(2), 1--29.

Baddeley A, Turner R (2000) Practical maximum pseudolikelihood for spatial point patterns.
Aust NZ J Stat 42, 283--322.

Baddeley A., Moller J. and Waagepetersen R. (2000). Non- and
semi-parametric estimation of interaction in inhomogeneous
point patterns. Statistica Neerlandica, 54, 329--350.

Diggle P. , Chedwynd A., Haggkvist R. and Morris S. (1995).
Second-order analysis of space-time clustering. Statistical
Methods in Medical Research, 4, 124--136. }

\author{
Edith Gabriel <edith.gabriel@univ-avignon.fr> }


\examples{\dontrun{
## First example
#################

data(fmd)
data(northcumbria)
FMD<-as.3dpoints(fmd[,1]/1000,fmd[,2]/1000,fmd[,3])
Northcumbria=northcumbria/1000

# estimation of the temporal intensity
Mt<-density(FMD[,3],n=1000)
mut<-Mt$y[findInterval(FMD[,3],Mt$x)]*dim(FMD)[1]

# estimation of the spatial intensity
h<-mse2d(as.points(FMD[,1:2]), Northcumbria, nsmse=50, range=4)
h<-h$h[which.min(h$mse)]
Ms<-kernel2d(as.points(FMD[,1:2]), Northcumbria, h, nx=5000, ny=5000)
atx<-findInterval(x=FMD[,1],vec=Ms$x)
aty<-findInterval(x=FMD[,2],vec=Ms$y)
mhat<-NULL
for(i in 1:length(atx)) mhat<-c(mhat,Ms$z[atx[i],aty[i]])

# estimation of the STIK function
u <- seq(0,10,by=1)
v <- seq(0,15,by=1)
stik1 <- STIKhat(xyt=FMD, s.region=northcumbria/1000,t.region=c(1,200), 
lambda=mhat*mut/dim(FMD)[1], dist=u, times=v, infectious=TRUE)

# plotting the estimation
plotK(stik1)
plotK(stik1,type="persp",theta=-65,phi=35)

## Second example
##################

xyt=rpp(lambda=200)
stik2=STIKhat(xyt$xyt,dist=seq(0,0.16,by=0.02),
times=seq(0,0.16,by=0.02),correction=c("border","translate"))

plotK(stik2,type="contour",legend=T,which="translate")
}}
