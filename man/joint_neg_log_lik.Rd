% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/joint_neg_log_lik.R
\name{joint_neg_log_lik}
\alias{joint_neg_log_lik}
\title{Joint negative log-likelihood function for several heteroplasmy measurements}
\usage{
joint_neg_log_lik(theta, hlist, use.h0s = F, h0s = -1)
}
\arguments{
\item{theta}{Kimura parameters p (or h0 here) and b.}

\item{hlist}{TBDD list of different sets of heteroplasmy measurements}

\item{use.h0s}{Logical parameter. TBD}

\item{h0s}{TBDD}
}
\value{
The negative log likelihood  for the list of inputs.
}
\description{
Joint negative log-likelihood function for several families' heteroplasmy
measurements.
theta = \code{[ b, h0.1, h0.2, ... ]} (use h values if use.h0s=F, otherwise initial heteroplasmies are enforced via h0s).
}
\examples{
 X.1 = rnorm(50,0.5,0.1)
joint_neg_log_lik(c(0.5,0.91),X.1)
}
\keyword{joint}
\keyword{kimura}
\keyword{likelihood}
\keyword{log}
\keyword{negative}
