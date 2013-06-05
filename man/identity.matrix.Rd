\name{identity.matrix}
\alias{identity.matrix}

\title{labelLevelsGOSetsHelpfile}
\description{
This function displays an object of class \code{character} specifying the names of the two levels compared in the experiment. It takes as input an object of class \code{\linkS4class{GOsims}}.
}
\usage{identity.matrix(object)}
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
identity.matrix(CCComparison)
}

\keyword{GOsims}
\keyword{identity.matrix}
