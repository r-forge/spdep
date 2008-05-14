\name{stsls}
\alias{stsls}
\alias{print.stsls}
\alias{print.summary.stsls}
\alias{summary.stsls}
\alias{residuals.stsls}
\alias{coef.stsls}
\alias{deviance.stsls}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Generalized spatial two stage least squares}
\description{
  The function fits a spatial lag model by two stage least squares, with the option of adjusting the results for heteroskedasticity.
}
\usage{
stsls(formula, data = list(), listw, zero.policy = FALSE,
 na.action = na.fail, robust = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{a symbolic description of the model to be fit. The details 
of model specification are given for \code{lm()}}
  \item{data}{an optional data frame containing the variables in the model. 
By default the variables are taken from the environment which the function 
is called.}
  \item{listw}{a \code{listw} object created for example by \code{nb2listw}}
  \item{zero.policy}{if TRUE assign zero to the lagged value of zones without 
neighbours, if FALSE (default) assign NA - causing \code{lagsarlm()} to terminate with an error}
  \item{na.action}{a function (default \code{na.fail}), can also be \code{na.omit} or \code{na.exclude} with consequences for residuals and fitted values - in these cases the weights list will be subsetted to remove NAs in the data. It may be necessary to set zero.policy to TRUE because this subsetting may create no-neighbour observations. Note that only weights lists created without using the glist argument to \code{nb2listw} may be subsetted.}
  \item{robust}{default FALSE, if TRUE, apply a heteroskedasticity correction to the coefficient estimatess and their covariances}
}
\details{
  The fitting implementation fits a spatial lag model:

\deqn{y = \rho W y + X \beta + \varepsilon}{y = rho W y + X beta + e}

by using spatially lagged X variables as instruments for the spatially lagged dependent variable.
}
\value{
  an object of class "stsls" containing:
  \item{coefficients}{coefficient estimates}
  \item{var}{coefficient covariance matrix}
  \item{sse}{sum of squared errors}
  \item{residuals}{model residuals}
  \item{df}{degrees of freedom}
}
\references{Kelejian, H.H. and I.R. Prucha (1998). A generalized spatial two
stage least squares procedure for estimating a spatial autoregressive
model with autoregressive disturbances. \emph{Journal of Real Estate
Finance and Economics} 17, 99-121.}
\author{Luc Anselin and Roger Bivand}

\seealso{\code{\link{lagsarlm}}}
\examples{
data(oldcol)
COL.lag.eig <- lagsarlm(CRIME ~ INC + HOVAL, data=COL.OLD, nb2listw(COL.nb))
summary(COL.lag.eig, correlation=TRUE)
COL.lag.stsls <- stsls(CRIME ~ INC + HOVAL, data=COL.OLD, nb2listw(COL.nb))
summary(COL.lag.stsls, correlation=TRUE)
COL.lag.stslsR <- stsls(CRIME ~ INC + HOVAL, data=COL.OLD, nb2listw(COL.nb),
robust=TRUE)
summary(COL.lag.stslsR, correlation=TRUE)
data(boston)
gp2a <- stsls(log(CMEDV) ~ CRIM + ZN + INDUS + CHAS + I(NOX^2) + I(RM^2) +
  AGE + log(DIS) + log(RAD) + TAX + PTRATIO + B + log(LSTAT),
 data=boston.c, nb2listw(boston.soi))
summary(gp2a)
}
\keyword{spatial}