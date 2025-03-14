% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions.R
\name{get_ANOVA_index}
\alias{get_ANOVA_index}
\title{Discriminative power weights (AN) function}
\usage{
get_ANOVA_index(x, y)
}
\arguments{
\item{x:}{matrix with continuous covariates}

\item{y:}{response vector as factor variable}
}
\value{
DP weight AN
}
\description{
Discriminative power weights (AN) function
}
\details{
Caclulates the weight of AN
}
\references{
Fuetterer, C., Nalenz, M., and Augustin, T. (2021).
Discriminative Power Lasso – Incorporating Discriminative Power of Genes into
Regularization-Based Variable Selection. Technical Report. Available under:
\url{https://epub.ub.uni-muenchen.de/91666/1/DPL_TR_2022_03.pdf}
}
