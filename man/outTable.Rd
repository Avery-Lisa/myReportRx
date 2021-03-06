% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/main.R
\name{outTable}
\alias{outTable}
\title{The output function for the print methods
Table output defaults to kable, but the kableExtra package doesn't work well with Word.
To export nice tables to Word use options('doc_type'='doc')}
\usage{
outTable(
  tab,
  to_indent = numeric(0),
  to_bold = numeric(0),
  caption = NULL,
  chunk_label,
  ...
)
}
\arguments{
\item{tab}{a table to format}

\item{to_indent}{numeric vector the length of nrow(tab) indicating which rows to indent}

\item{to_bold}{numeric vector the length of nrow(tab) indicating which rows to bold}

\item{caption}{table caption}

\item{chunk_label}{only used if out_fmt = doc to allow cross-referencing}

\item{...}{other variables passed to covsum and the table output function}
}
\description{
The output function for the print methods
Table output defaults to kable, but the kableExtra package doesn't work well with Word.
To export nice tables to Word use options('doc_type'='doc')
}
