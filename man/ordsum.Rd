% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/main.R
\name{ordsum}
\alias{ordsum}
\title{Fit and format an ordinal logistic regression using polr from the {MASS} package. The parallel regression assumption can
be tested using the Brant test (modfied from the the Brant package).}
\usage{
ordsum(
  data,
  covs,
  response,
  reflevel,
  markup = FALSE,
  sanitize = TRUE,
  nicenames = TRUE,
  excludeLevels,
  testPO = TRUE,
  showN = TRUE,
  digits = 1,
  CIwidth = 0.95
)
}
\arguments{
\item{data}{dataframe containing data [REQUIRED]}

\item{covs}{character vector with the names of columns to include in table [REQUIRED]}

\item{response}{ordinal outcome variable [REQUIRED]}

\item{reflevel}{manual specification of the reference level, must match level exactly}

\item{markup}{boolean indicating if you want latex markup}

\item{sanitize}{boolean indicating if you want to sanitize all strings to not break LaTeX}

\item{nicenames}{booling indicating if you want to replace . and _ in strings with a space}

\item{excludeLevels}{a named list of levels to exclude from the response variable}

\item{testPO}{logical, should the proportional odds (parallel regression) assumption be tested with the Brant test, defaults to TRUE, values greater than alpha are desirable.}

\item{showN}{logical, should sample sizes be shown for each lvel, defaults to TRUE}

\item{digits}{number of digits to display, defaults to}

\item{CIwidth}{level of significance for computing the confidence intervals, default is 0.95}
}
\value{
A formatted table displaying the odds ratio associated with each covariate
}
\description{
Fit and format an ordinal logistic regression using polr from the {MASS} package. The parallel regression assumption can
be tested using the Brant test (modfied from the the Brant package).
}
\keyword{Brant}
\keyword{ordinal}
\keyword{regression,}
\keyword{test}
