\name{label.level.enriched}
\alias{label.level.enriched}

\title{labelLevelsEnrichedSetsHelpfile}
\description{
This function displays an object of class \code{character} specifying the names of the two levels compared in the experiment. It takes as input an object of class \code{\linkS4class{EnrichedSets}} or \code{\linkS4class{GOsets}}.
}
\usage{label.level.enriched(object)}
\arguments{
 \item{object}{an object of class \code{\linkS4class{EnrichedSets}} or \code{\linkS4class{GOsets}}.}
 }

\author{
	Toma Tebaldi, Erik Dassi, Galena Kostoska
}

\seealso{
\code{\linkS4class{EnrichedSets}}
\code{\link{RegulatoryEnrichment}}
\code{\linkS4class{GOsets}}
\code{\link{GOEnrichment}}
}

\examples{
data(tRanslatomeSampleData)
label.level.enriched(CCEnrichment)
}

\keyword{EnrichedSets}
\keyword{GOsets}
\keyword{label.level.enriched}
