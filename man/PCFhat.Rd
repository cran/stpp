\name{PCFhat} \alias{PCFhat}
\title{Estimation of the space-time pair correlation function}

\description{Compute an estimate of the space-time pair
correlation function.}

\usage{PCFhat(xyt, s.region, t.region, dist, times, lambda,
ks="box", hs, kt="box", ht, correction = "isotropic") }

\arguments{
  \item{xyt}{coordinates and times (x,y,t) of the point pattern.}
    \item{s.region}{two-column matrix specifying polygonal region containing
  all data locations.
   If \code{s.region} is missing, the bounding box of \code{xyt[,1:2]} is considered.}
  \item{t.region}{vector containing the minimum and maximum values of
  the time interval.
   If \code{t.region} is missing, the range of \code{xyt[,3]} is considered.}
  \item{dist}{vector of distances u at which \eqn{g(u,v)}{g(u,v)} is computed.}
   \item{times}{vector of times v at which \eqn{g(u,v)}{g(u,v)} is computed.}
    \item{lambda}{vector of values of the space-time intensity function
  evaluated at the points (x,y,t) in SxT.
  If \code{lambda} is missing, the estimate of the space-time
pair correlation function is computed considering \eqn{n/|S
\times T|}{n/|SxT|} as an estimate of the space-time
intensity.}

\item{ks}{Kernel function for the spatial distances. Default is
the \code{"box"} kernel. Can also be \code{"epanech"} for the
Epanechnikov kernel or \code{"gaussian"} or \code{"biweight"}.}

\item{hs}{Bandwidth of the kernel function \code{ks}.}

\item{kt}{Kernel function for the temporal distances. Default
is the \code{"box"} kernel. Can also be \code{"epanech"} for
the Epanechnikov kernel or \code{"gaussian"} or
\code{"biweight"}.}

\item{ht}{Bandwidth of the kernel function \code{kt}.}

\item{correction}{A character vector specifying the edge correction(s) to be applied 
among \code{"isotropic"}, \code{"border"}, \code{"modified.border"}, \code{"translate"}
 and \code{"none"} (see Details). The default is \code{"isotropic"}.}
 }

\details{An approximately unbiased estimator for the space-time
pair correlation function, based on data giving the locations
of events \eqn{x_i: i=1,...n}{xi: i=1...,n} on a
spatio-temporal region \eqn{S \times T}{SxT}, where S is an
arbitrary polygon and T a time interval: 
\deqn{\widehat g(u,v) = \frac{1}{4 \pi u} \sum_{i=1}^{n}
  \sum_{j \neq i} \frac{1}{w_{ij}} 
  \frac{k_{s}(u-\|s_i-s_j \|)k_{t}(v-|t_i-t_j|)}{\lambda(x_i) \lambda(x_j)}}{g(u,v) = 1/|SxT| 1/(4 pi u) sum_{i=1,...,n}
  sum_{j=1,...,n; j != i} 1/wij  ks(u -
  ||si-sj||)kt(v-|ti-tj|)/(lambda(xi)lambda(xj))}
where \eqn{\lambda(x_i)}{lambda(xi)} is the intensity at \eqn{x_i = (s_i,t_i)}{xi=(si,ti)} and
\eqn{w_{ij}}{wij} is an edge correction factor to deal with spatial-temporal edge effects. 
The edge correction methods implemented are:

\code{isotropic}: \eqn{w_{ij} = |S \times T| w_{ij}^{(t)} w_{ij}^{(s)}}{wij = |S x T| wij(s) wij(t)}, 
where the temporal edge correction factor \eqn{w_{ij}^{(t)} = 1}{wij(t)=1} if both ends of the interval 
of length \eqn{2 |t_i - t_j|}{2|ti-tj|} centred at \eqn{t_i}{ti} lie within T and \eqn{w_{ij}^{(t)}=1/2}{wij(t)=1/2} 
otherwise and \eqn{w_{ij}^{(s)}}{wij(s)} is the proportion of the circumference of a circle centred 
at the location \eqn{s_i}{si} with radius \eqn{\|s_i -s_j\|}{||si-sj||} lying in S (also called Ripley's
edge correction factor).

\code{border}: \eqn{w_{ij} = \frac{\sum_{j=1}^n {\bf 1}_{\lbrace d(s_j,S) > u \ ; \ d(t_j,T) >v\rbrace} / 
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


\eqn{k_s()}{ks()} and \eqn{k_t()}{kt()} denotes kernel
functions with bandwidth \eqn{h_s}{hs} and \eqn{h_t}{ht}.
Experience with pair correlation function estimation recommends
box kernels (the default), see Illian et al. (2008).
Epanechnikov, Gaussian and biweight kernels are also
implemented. Whatever the kernel function, if the bandwidth is
missing, a value is obtain from the function \code{dpik} of the
package KernSmooth. Note that the bandwidths play an important
role and their choice is crucial in the quality of the
estimators as they heavily influence their variance. }

 \value{A list containing:

 \item{pcf}{ndist x ntimes matrix containing values of \eqn{\hat
 g(u,v)}{g(u,v)}.}

 \item{dist, times}{parameters passed in argument.}

 \item{kernel}{a vector of names and bandwidths of the spatial and temporal kernels.}
 
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

Illian JB, Penttinen A, Stoyan H and Stoyan, D.
(2008). Statistical Analysis and Modelling of Spatial Point
Patterns. John Wiley and Sons, London.
 }

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

# estimation of the pair correlation function
g1 <- PCFhat(xyt=FMD, dist=1:20, times=1:20, lambda=mhat*mut/dim(FMD)[1],
 s.region=northcumbria/1000,t.region=c(1,200))

# plotting the estimation 

plotPCF(g1)
plotPCF(g1,type="persp",theta=-65,phi=35) 

## Second example
###################

xyt=rpp(lambda=200)
g2=PCFhat(xyt$xyt,dist=seq(0,0.16,by=0.02),
times=seq(0,0.16,by=0.02),correction=c("border","translate"))

plotPCF(g2,type="contour",which="border")

}}
