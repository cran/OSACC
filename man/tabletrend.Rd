\name{tabletrend}
\alias{tabletrend}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Cochran-Armitage test for trend
}
\description{
named for William Cochran and Peter Armitage, is used in categorical data analysis when the aim is to assess for the presence of an association between a variable with two categories and a variable with k categories}
\usage{
tabletrend(x, transpose = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{the array to be used. 
}
  \item{transpose}{"Cochran-Armitage test for trend must be used with a (R,2) table. Use transpose argument",call.=FALSE
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
\item{estimate }{the value of the statistic.}
\item{dim }{dimension of  the table X}
\item{p.value.uni}{the p-value for the test from one.sided}
\item{p.value.bi }{the p-value for the test from two.sided}

}
\references{Agresti, Alan (2002). Categorical Data Analysis (Second Edition). Wiley. ISBN 0-471-36093-7.
}
\author{A code is from Eric Lecoutre on R -help Archive http://tolstoy.newcastle.edu.au/R/help/05/07/9442.html
}

