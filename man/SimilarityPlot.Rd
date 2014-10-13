\name{SimilarityPlot}
\alias{SimilarityPlot}

\title{SimilarityPlotHelpfile}
\description{
This function displays in a barplot, for each GO ontology, the average semantic similarity value between GO terms showing enrichment at the first or at the second level of analysis.
}
\usage{SimilarityPlot(object, outputformat="on screen")}
\arguments{
   \item{object}{an object of class \code{\linkS4class{GOsims}}}
   \item{outputformat}{a character string specifying if the plot is saved in jpeg (\code{jpeg}), postscript (\code{postscript}), pdf (\code{pdf}) format, or it is simply displayed on the screen(\code{on screen}). By default this value is \code{on screen}.}
   }
\author{
	Toma Tebaldi, Erik Dassi, Galena Kostoska
}

\seealso{
\code{\link{GOComparison}}
\code{\linkS4class{GOsets}}
\code{\linkS4class{GOsims}}
}

\examples{
data(tRanslatomeSampleData)
SimilarityPlot(CCComparison, outputformat="on screen")
}

\keyword{GOsims}
\keyword{GOcomparison}
\keyword{GOsets}
\keyword{GOenrichment}

