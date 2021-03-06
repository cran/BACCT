% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BAC_binom.R
\name{BAC_binom}
\alias{BAC_binom}
\title{Bayesian Augmented Control for Binary Responses}
\usage{
BAC_binom(yh, nh, n1, n2, y1.range = 0:n1, y2.range = 0:n2, n.chain = 5,
  tau.alpha = 0.001, tau.beta = 0.001, prior.type = "nonmixture",
  criterion.type = c("diff", "prob"), prob.threshold, sim.mode = c("full",
  "express"))
}
\arguments{
\item{yh, nh}{Vector of the numbers of events (subjects) in the historical
trial(s). Must be of equal length.}

\item{n1, n2}{Number of subjects in the control or treatment arm of the current
trial.}

\item{y1.range, y2.range}{Number of events in control or treatment arm of the
current trial. See "Details".}

\item{n.chain}{Controls the number of posterior samples. Each chain contains
20,000 samples.}

\item{tau.alpha, tau.beta}{Hyperparameters of the inverse gamma distribution
controling the extent of borrowing.}

\item{prior.type}{Type of prior on control groups. Currenly, only the
inverse-gamma prior is implemented.}

\item{criterion.type}{Type of posterior quantities to be monitored. See
"Details."}

\item{prob.threshold}{For \code{criterion.type="prob"} only. See "Details".}

\item{sim.mode}{Simulation duration reduces greatly in \code{"express"}
mode, if treatment and control arms are independent. See "Details".}
}
\value{
An object of class "BAC".
}
\description{
Calling JAGS to implement BAC for binary responses
}
\details{
There are two types of posterior quantities for
  \code{criterion.type} argument. With \code{"diff"} option, the quantity
  computed is \eqn{p_{T} - p_{C}}; with \code{"prob,"} such quantity is
  \eqn{pr(p_{T} - p_{C}>\Delta)}, where \eqn{\Delta} is specified by
  \code{prob.threshold} argument.

  By default, \code{y1.range} and \code{y2.range} cover all possible outcomes
  and should be left unspecified in most cases. However, when \code{n1}
  and/or \code{n2} is fairly large, it is acceptable to use a reduced range
  that covers the outcomes that are most likely (e.g., within 95\% CI) to be
  observed. This may help shorten the time to run MCMC.

  Another way that can greatly shorten the MCMC running time is to specify
  \code{"express"} mode in \code{sim.mode} argument. Express mode reduces the
  number of simulations from \code{length(y1.range)*length(y2.range)} to
  \code{length(y1.range)+length(y2.range)}. Express mode is proper when the
  treatment arm rate is independent of control arm rate.
}
\examples{
\dontrun{
library(BACCT)
#borrow from 3 historical trials#
yh = c(11,300,52);nh = c(45,877,128)
#specify current trial sample sizes#
n1 = 20;n2 = 30

#Difference criterion type in full simulation mode#
obj1 = BAC_binom(yh=yh,nh=nh,n1=n1,n2=n2,n.chain=5,
criterion.type="diff",sim.mode="full")

#Probability criterion type in express simulation mode#
obj2 = BAC_binom(yh=yh,nh=nh,n1=n1,n2=n2,n.chain=5,
criterion.type="prob",prob.threshold=0.1,sim.mode="express")

#S3 method for class "BAC"
summary(obj1)
}
}
\author{
Hongtao Zhang
}

