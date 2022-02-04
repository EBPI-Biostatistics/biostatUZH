
\name{combineSubgroups}
\alias{combineSubgroups}
\title{Combines summary data across subgroups}

\description{Combines summary data for a continuous outcome given in subgroups}

\usage{
combineSubgroups(n, means, variances, treatment)
}

\arguments{
  \item{n}{vector with sample sizes in each subgroup.}
  \item{means}{vector with sample means}
  \item{variances}{vector with sample variances}
  \item{treatment}{Indicates which entries of each vector belong to the treatment group}
}

\value{
  Gives overall difference in means with a standard error.
}

\author{Leonhard Held}

\examples{
combineSubgroups(n=c(10,20,30,40), means=c(12,11,10,9),
variances=c(3,4,3,4), treatment=c(2,4)) 
}

