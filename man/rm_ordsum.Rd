% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/main.R
\name{rm_ordsum}
\alias{rm_ordsum}
\title{Fit and format an ordinal logistic regression using polr from the {MASS} package. The parallel regression assumption can
be tested using the Brant test in the Brant package and visually. Only logistic ordinal regression is supported currently.}
\usage{
rm_ordsum(
  data,
  covs,
  response,
  reflevel = "NULL",
  caption = NULL,
  showN = T,
  mv = FALSE,
  excludeLevels = NULL,
  testPO = TRUE,
  digits = 2,
  CIwidth = 0.95
)
}
\arguments{
\item{data}{dataframe containing data [REQUIRED]}

\item{covs}{character vector with the names of columns to include in table [REQUIRED]}

\item{response}{ordinal outcome variable [REQUIRED]}

\item{reflevel}{manual specification of the reference level, must match level exactly}

\item{caption}{Table caption}

\item{showN}{logical, should sample sizes be shown for each lvel, defaults to TRUE}

\item{mv}{logical indicating whether you want to summarise univariate analyses or a single model only mv=FALSE currently supported}

\item{excludeLevels}{a named list of levels to exclude from factor variables. Currently, this has only been implemented for the response variable.}

\item{testPO}{logical, should the proportional odds (parallel regression) assumption be tested with the Brant test, defaults to TRUE}

\item{digits}{number of digits to display, defaults to}

\item{CIwidth}{level of significance for computing the confidence intervals, default is 0.95}
}
\value{
A formatted table displaying the odds ratio associated with each covariate
}
\description{
Fit and format an ordinal logistic regression using polr from the {MASS} package. The parallel regression assumption can
be tested using the Brant test in the Brant package and visually. Only logistic ordinal regression is supported currently.
}
\keyword{Brant}
\keyword{ordinal}
\keyword{regression,}
\keyword{test}