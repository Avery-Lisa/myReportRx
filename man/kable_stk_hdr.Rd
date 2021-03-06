% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/main.R
\name{kable_stk_hdr}
\alias{kable_stk_hdr}
\title{Stack columns in a table for clearer viewing}
\usage{
kable_stk_hdr(
  data,
  head_col,
  to_col,
  caption = NULL,
  indent = TRUE,
  hdr_prefix = "",
  hdr_suffix = "",
  tblOnly = FALSE
)
}
\arguments{
\item{data}{dataframe}

\item{head_col}{character value specifying the column name with the headers}

\item{to_col}{character value specifying the column name to add the headers into}

\item{caption}{table caption}

\item{indent}{should the original values in the to_col be indented}

\item{hdr_prefix}{character value that will prefix headers}

\item{hdr_suffix}{character value that will suffix headers}

\item{tblOnly}{boolean indicating if the table should be formatted for prining or returned as a data frame}
}
\description{
Stack columns in a table for clearer viewing
}
