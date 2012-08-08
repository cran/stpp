\name{PCFhat} \alias{PCFhat}
\title{Estimation of the space-time pair correlation function}

\description{Compute an estimate of the space-time pair
correlation function.}

\usage{PCFhat(xyt, s.region, t.region, dist, times, lambda,
ks="box", hs, kt="box", ht, correction = TRUE) }

\arguments{
  \item{xyt}{coordinates and times (x,y,t) of the point pattern.}
    \item{s.region}{two-column matrix specifying polygonal region containing
  all data locations.
   If \code{s.region} is missing, the bounding box of \code{xyt[,1:2]} is considered.}
  \item{t.region}{vector containing the minimum and maximum values of
  the time interval.
   If \code{t.region} is missing, the range of \code{xyt[,3]} is considered.}
  \item{dist}{vector of distances u at which \eqn{K_{ST}(u,v)}{K(u,v)} is computed.}
   \item{times}{vector of times v at which \eqn{K_{ST}(u,v)}{K(u,v)} is computed.}
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


\item{correction}{logical value. If \code{TRUE}, spatial
(Ripley's) and temporal edge corrections are used.}
 }

\details{An approximately unbiased estimator for the space-time
pair correlation function, based on data giving the locations
of events \eqn{x_i: i=1,...n}{xi: i=1...,n} on a
spatio-temporal region \eqn{S \times T}{SxT}, where S is an
arbitrary polygon and \eqn{T=[T_0,T_1]}{T=[T0,T1]}: \deqn{
  \widehat g(u,v) =  \frac{1}{|S \times T|} \sum_{i=1}^{n}
  \sum_{j \neq j}^{n} \frac{1}{w_{ij} v_{ij}} \frac{k_{s}(u - \|s_i-s_j\|)k_{t}(v-|t_i-t_j|)}{\lambda(x_i) \lambda(x_j)}}{g(u,v) = 1/|SxT| sum_{i=1,...,n}
  sum_{j=1,...,n; j \neq j} 1/(wij*vij)  ks(u -
  ||si-sj||)kt(v-|ti-tj|)/(lambda(xi)lambda(xj))}
To deal with spatial edge-effects, we use Ripley's method, in
which \eqn{w_{ij}}{wij} is the proportion of the circle
centered on \eqn{s_i}{si} and passing through \eqn{s_j}{sj},
i.e. of radius \eqn{u_{ij}= \| s_i - s_j\|}{uij=||si-sj||},
that lies inside S. To deal with temporal edge effects,
\eqn{v_{ij}}{vij} is equal to 1 if both ends of the
  interval of length \eqn{2 |t_i - t_j|}{2|ti-tj|} centred at \eqn{t_i}{ti} lie within
  T and 1/2 otherwise.

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
 }

\references{ Illian JB, Penttinen A, Stoyan H and Stoyan, D.
(2008). Statistical Analysis and Modelling of Spatial Point
Patterns. John Wiley and Sons, London.
 }

\author{
Edith Gabriel <edith.gabriel@univ-avignon.fr> }


\examples{\dontrun{
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
g <- PCFhat(xyt=FMD, dist=1:20, times=1:20, lambda=mhat*mut/dim(FMD)[1],
 s.region=northcumbria/1000,t.region=c(1,200))

# plotting the estimation plotPCF(g)
plotPCF(g,persp=TRUE,theta=-65,phi=35) }}
