% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mediateP.R
\name{mediate_binaY_contM}
\alias{mediate_binaY_contM}
\title{Mediation analysis function for binary outcome and continuous mediator}
\usage{
mediate_binaY_contM(
  data,
  outcome = "Y",
  mediator = "M",
  exposure = "X",
  covariateY = c("X1", "X2"),
  covariateM = c("X1", "X2"),
  x0 = 0,
  x1 = 1,
  cY = c(0, 0),
  cM = c(0, 0)
)
}
\arguments{
\item{data}{A dataset.}

\item{outcome}{The outcome variable.}

\item{mediator}{The mediator variable.}

\item{exposure}{The exposure variable.}

\item{covariateY}{A vector of confounders in the outcome regression.}

\item{covariateM}{A vector of confounders in the mediator regression.}

\item{x0}{The baseline exposure level.}

\item{x1}{The new exposure level.}

\item{cY}{conditional levels of covariateY}

\item{cM}{conditional levels of covariateM}
}
\value{
A list containing NIE, NDE and MP point and interval estimates.
}
\description{
Mediation analysis function for binary outcome and continuous mediator
}
