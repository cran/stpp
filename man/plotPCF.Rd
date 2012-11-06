\name{plotPCF} \alias{plotPCF}
\title{Plot the estimation of the Space-Time Inhomogeneous pair correlation function}

\description{Contour plot or perspective plot the Space-Time
Inhomogeneous pair correlation function estimate.}

\usage{plotPCF(PCF,n=15,persp=FALSE,legend=TRUE,...) }

\arguments{
  \item{PCF}{result of the \code{PCFhat} function.}
    \item{n}{number of contour levels desired.}
   \item{persp}{logical indicating whether a contour plot (default) or a perspective plot must be plotted.}
    \item{legend}{logical indicating whether a legend must be added to the plot.}
\item{...}{additional arguments to \code{persp} if \code{persp=TRUE}, such as \code{theta} and \code{phi}.}
 }

\seealso{
 \code{\link{contour}}, \code{\link{persp}} and \code{\link{PCFhat}} for an example.
 }

\author{
Edith Gabriel <edith.gabriel@univ-avignon.fr> }


