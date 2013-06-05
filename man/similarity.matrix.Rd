\name{similarity.matrix}
\alias{similarity.matrix}

\title{similaritymatrixGOSetsHelpfile}
\description{
This function displays an object of class \code{character} specifying the names of the two levels compared in the experiment. It takes as input an object of class \code{\linkS4class{GOsims}}.
}
\usage{similarity.matrix(object)}
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
similarity.matrix(CCComparison)
}

\keyword{GOsims}
\keyword{similarity.matrix}