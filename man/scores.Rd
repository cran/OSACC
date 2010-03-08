\name{scores}
\alias{scores}
\title{Scores}
\description{ It is for Cochran-Armitage-trend  test
}
\usage{
scores(x, MARGIN = 1, method = "table", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ the array to be used }
  \item{MARGIN}{ 1 indicates rows, 2 indicates columns 
}
  \item{method}{Methods for ranks are default, the other option are ridit,modridit
}
\item{\dots}{... optional arguments
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
\item{out  }{a score from the table}
}
\references{Agresti, Alan (2002). Categorical Data Analysis (Second Edition). Wiley. ISBN 0-471-36093-7.
}
\author{A code is from Eric Lecoutre on R -help Archive http://tolstoy.newcastle.edu.au/R/help/05/07/9442.html
}


