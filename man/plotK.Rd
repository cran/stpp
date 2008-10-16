\name{plotK} \alias{plotK}
\title{Plot the estimation of the Space-Time Inhomogeneous K-function}

\description{Contour plot or perspective plot the Space-Time
Inhomogeneous K-function estimate.}

\usage{plotK(K,n=15,L=TRUE,persp=FALSE,legend=TRUE,...) }

\arguments{
  \item{K}{result of the \code{STIkhat} function.}
    \item{n}{number of contour levels desired.}
  \item{L}{logical indicating whether \eqn{K_{ST}(u,v)}{K(u,v)} or \eqn{L(u,v)=K_{ST}(u,v)-\pi u^2 v}{K(u,v)-pi u^2 v} must be plotted.}
   \item{persp}{logical indicating whether a contour plot (default) or a perspective plot must be plotted.}
    \item{legend}{logical indicating whether a legend must be added to the plot.}
\item{...}{additional arguments to \code{persp} if \code{persp=TRUE}, such as \code{theta} and \code{phi}.}
 }

\seealso{
 \code{\link{contour}}, \code{\link{persp}} and \code{\link{STIKhat}} for an example.
 }

\author{
Edith Gabriel <edith.gabriel@univ-avignon.fr> }


