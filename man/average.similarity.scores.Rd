\name{average.similarity.scores}
\alias{average.similarity.scores}

\title{labelLevelsGOSetsHelpfile}
\description{
This function displays an object of class \code{character} specifying the names of the two levels compared in the experiment. It takes as input an object of class \code{\linkS4class{GOsims}}.
}
\usage{average.similarity.scores(object)}
\arguments{
 \item{object}{an object of class \code{\linkS4class{GOsims}}.}
 }

\author{
	Toma Tebaldi, Erik Dassi, Galena Kostoska
}

\seealso{
\code{\linkS4class{GOsims}}
\code{\link{GOComparison}}
}

\examples{
data(tRanslatomeSampleData)
average.similarity.scores(CCComparison)
}

\keyword{GOsims}
\keyword{average.similarity.scores}
