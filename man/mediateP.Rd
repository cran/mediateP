\name{mediate}
\alias{mediate}
%- R functions for calculating the point and interval estimates for the NIE, TE and MP, based on the product method.
\title{
Mediation analysis based on the product method
}
\description{
Function for calculating the point and interval estimates for the NIE, TE and MP, based on the product method.
}
\usage{
mediate(data, outcome, mediator, exposure, binary.outcome,
        binary.mediator, covariate.outcome, covariate.mediator,
        x0, x1, c.outcome, c.mediator, boot, R)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{(Required) The name of the dataset.}
  \item{outcome}{(Required) Name of the outcome variable, which should be either a continuous or binary datatype.}
  \item{mediator}{(Required) Name of the mediator variable, which should be either a continuous or binary datatype.}
  \item{exposure}{(Required) Name of the exposure variable, which should be either a continuous or binary datatype.}
  \item{binary.outcome}{(Required) If the outcome is binary, set to 1. If the outcome is continuous, set to 0.}
  \item{binary.mediator}{(Required) If the mediator is binary, set to 1. If the mediator is continuous, set to 0.}
  \item{covariate.outcome}{A vector of names showing the confounding variables used in the outcome model. The default value is NULL, which represents no confounding variables.  We only accepted continuous and binary confounding variables, if one confounding variable is categorical, please set it to a series of binary variables in advance.}
  \item{covariate.mediator}{A vector of names showing the confounding variables used in the mediator model. The default value is NULL, which represents no confounding variables. We only accepted continuous and binary confounding variables, if one confounding variable is categorical, please set it to a series of binary variables in advance.}
  \item{x0}{The baseline exposure level (i.e., $x^*$). The default value is 0.}
  \item{x1}{The new exposure level (i.e., $x$). The default value is 1.}
  \item{c.outcome}{A vector of numbers representing the conditional level of the confounding variables in the outcome model. The default value is a vector of 0.}
  \item{c.mediator}{A vector of numbers representing the conditional level of the confounding variables in the mediator model. The default value is a vector of 0.}
  \item{boot}{If a percentile bootstrap confidence interval needed to be added, set to 1. Otherwise, set to 0. The default value is 0.}
  \item{R}{(Required if using bootstrap approach) The number of replications when apply the percentile bootstrap method to calculate the confidence interval. The default value is 2,000.}
}

\value{
A list presents the mediation analysis results, which includes two components. The first component `res' presents the point estimator, standard error, and the 95\% confidence interval for the NIE, TE, and MP; the second component `class' clarifies this is a `mediate' object.
}

\author{
Chao Cheng (c.cheng@yale.edu)
}

\references{Chao Cheng, Donna Spiegelman, and Fan Li. (2020+) An evaluation of the product method for estimation and inference in mediation analysis. Submitted to Epidemiologic Methods. Manuscript is available at <https://github.com/chaochengstat/chaochengstat.github.io/blob/master/manuscript_Epidemiologic_Methods.pdf>.}

\examples{
## Generate a data set including two confounders (C1 and C2), one binary
## mediator (M), one binary exposure (X), one binary outcome (Y1), and
## one continuous outcome (Y2)
C1=rnorm(2000)>0
C2=rnorm(2000)
X = rnorm(2000)
M= as.numeric(runif(2000)< 1/(1+exp(0-0.9*X+0.1*C1)))
Y1= as.numeric(runif(2000)< 1/(1+exp(-(-2+0.5*X+0.5*M+0.2*C1+0.2*C2))))
Y2= -2+0.9*X-0.6*M+0.2*C1+0.2*C2+rnorm(2000,mean=0,sd=1)
mydata=as.data.frame(cbind(Y1,Y2,M,X,C1,C2))

## Example 1: Mediation analysis for X in change from 0 to 1 on Outcome Y1.
## In the outcome model, we adjusted for C1 and C2. In the mediator model,
## we only adjusted for C1. We calculated the NIE, TE and MP for exposure in
## change from 0 to 1, conditional on C1=0 and C2=1.
result=mediate(data=mydata, outcome="Y1", mediator="M", exposure="X",
                 binary.outcome=1, binary.mediator=1,
                 covariate.outcome=c("C1","C2"),
                 covariate.mediator=c("C1"),
                 x0=0, x1=1, c.outcome=c(0,1), c.mediator=c(0),
                 boot=0)

print(result)

## Example 2: Mediation analysis for X in change from 0 to 1 on Outcome Y2.
## In the outcome model, we adjusted for C1 and C2. In the mediator model,
## we did not adjust for any confounders. We calculated the NIE, TE and MP
## for exposure in change from 0 to 1, conditional on C1=0 and C2=1.

result=mediate(data=mydata, outcome="Y2", mediator="M", exposure="X",
                 binary.outcome=0, binary.mediator=1,
                 covariate.outcome=c("C1","C2"),
                 covariate.mediator=NULL,
                 x0=0, x1=1, c.outcome=c(0,1), c.mediator=NULL,
                 boot=0)

print(result)

## Example 3: Mediation analysis without any confounding adjustment
## We conduct mediation analysis for X in change from 0 to 1 on Outcome
## Y1, but here no confounders in the outcome and mediator models were
## considered.

result=mediate(data=mydata, outcome="Y1", mediator="M", exposure="X",
                 binary.outcome=1, binary.mediator=1,
                 covariate.outcome=NULL,
                 covariate.mediator=NULL,
                 x0=0, x1=1, c.outcome=NULL, c.mediator=NULL,
                 boot=0)

print(result)

}

