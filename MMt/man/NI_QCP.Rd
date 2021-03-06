% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/NI_QCP.R
\name{NI_QCP}
\alias{NI_QCP}
\title{Cauchy mixture data generation}
\usage{
NI_QCP(y, m, general_function, threshold = 0.001, weight_omit = 0.01)
}
\arguments{
\item{y}{mixture data}

\item{m}{initial guess of mixture components}

\item{general_function}{general pdf function of components comes from}

\item{threshold}{relative change ratio of loss}

\item{weight_omit}{threshold for too small weight, default is 0.01}
}
\value{
number of components, location paramters, scale parameters, and weight paramters
}
\description{
Cauchy mixture data generation
}
