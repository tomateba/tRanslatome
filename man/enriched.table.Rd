\name{enriched.table}
\alias{enriched.table}

\title{enrichedtableenrichedSetsHelpfile}
\description{
This function displays an object of class \code{character} specifying the names of the two levels compared in the experiment. It takes as input an object of class \code{\linkS4class{EnrichedSets}}.
}
\usage{enriched.table(object)}
\arguments{
 \item{object}{an object of class \code{\linkS4class{EnrichedSets}}.}
 }

\author{
	Toma Tebaldi, Erik Dassi, Galena Kostoska
}

\seealso{
\code{\linkS4class{EnrichedSets}}
\code{\link{GOEnrichment}}
}

\examples{
data(tRanslatomeSampleData)
enriched.table(CCEnrichment)
}

\keyword{EnrichedSets}
\keyword{enriched.table}
