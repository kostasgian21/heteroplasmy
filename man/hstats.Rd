% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hstats.R
\name{hstats}
\alias{hstats}
\title{Calculate various statistics for a heteroplasmy set h}
\usage{
hstats(h, h0 = F, usepopn = F)
}
\arguments{
\item{h}{A vector containig heteroplasmy measurements. Every observation should be in \code{[0,1]}.}

\item{h0}{Logical parameter. A particular h0 value  Default is to treat h0 as a fit parameter}

\item{usepopn}{Logical parameter. Use of population or sample statistics (T and F, respectively)}
}
\value{
The maximum likelihood for the input data according to the Kimura distribution (using bootstrapping)
}
\description{
Can enforce an initial \emph{h0} or leave as a free parameter. Can can use population or sample statistics. \emph{analyticVar} offers a simplified version of this function to compute the standard error of the variance.
}
\examples{
 X.1 = rnorm(50,0.5,0.1)
hstats(X.1)
}
\keyword{statistic}
\keyword{summary}
