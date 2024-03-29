% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/kimura_neg_loglik.R
\name{kimura_neg_loglik}
\alias{kimura_neg_loglik}
\title{Calculate negative log-likelihood for heteroplasmy measurements}
\usage{
kimura_neg_loglik(theta, h, h0 = F)
}
\arguments{
\item{theta}{Kimura parameters p (or h0 here) and b.}

\item{h}{A vector containig heteroplasmy measurements. Every observation should be in \code{[0,1]}.}

\item{h0}{Logical parameter. A particular h0 value  Default is to treat h0 as a fit parameter}
}
\value{
The negative log likelihood  for the input.
}
\description{
Calculate negative log-likelihood for a given set of heteroplasmy measurements
\emph{h} and parameters theta = {logit(p), logit(b)} (we write \emph{h0} for \emph{p}). We can do this by
enforcing a particular \emph{h0} value (passed as an argument) or treating \emph{h0} as a fit
parameter (default). The logit transform is used to ensure \emph{h0} and \emph{b} remain in
the \code{[0,1]} interval regardless of what real-valued argument the numerical optimiser
attempts.
}
\examples{
 X.1 = rnorm(50,0.5,0.1)
kimura_neg_loglik(c(0.5,0.91),X.1)
}
\keyword{kimura}
\keyword{likelihood}
\keyword{log}
\keyword{negative}
