% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/showDendro.R
\name{showDendro}
\alias{showDendro}
\title{dendrogram from unsupervised clustering}
\usage{
showDendro(M, pheno, groupPal, topVar = 10 * 1000L, verbose = TRUE,
  ttl = "", ...)
}
\arguments{
\item{M}{(matrix) beta values}

\item{pheno}{(data.frame) row order should match column order of M.
should contain column names for groupBy. Note that groupBy columns 
should be factors}

\item{groupPal}{(list) keys are columns in pheno by which data are to
be grouped. Values are RColorBrewer palette names}

\item{topVar}{(integer) top most-variable probes to show}

\item{verbose}{(logical) print messages}

\item{...}{parameters for squash::dendromat}
}
\value{
No value. Side effect of printing a dendrogram)
}
\description{
dendrogram from unsupervised clustering
}

