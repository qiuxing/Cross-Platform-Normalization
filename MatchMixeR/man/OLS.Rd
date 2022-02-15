\name{OLS}
\alias{OLS}

\title{A cross-platform normalization procedure for matched samples
  based on OLS regression.}

\description{ This function takes in two matched gene expression
  matrices X (first platform) and Y (second platform). It runs OLS
  regression for each gene and uses the resulting simple linear model to
  transform X into the corresponding gene expressions on the second
  platform.
}

\usage{OLS(Xmat, Ymat)}

\arguments{
  \item{X}{gene expression level matrix on platform X (to be transformed)}

  \item{Y}{gene expression level matrix on platform Y ( not to be transformed)}
}

\details{
  This function fits the OLS regression model to transform X to Y.
}

\value{8 lists including betahat, betamat, Yhat}

\author{Xing Qiu}

\examples{
fit_OLS <- OLS(gpl96,gpl570)
}

\keyword{MatchMixeR}
\keyword{flmer}
