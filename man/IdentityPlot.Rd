\name{IdentityPlot}
\alias{IdentityPlot}

\title{IdentityPlotHelpfile}
\description{
This function displays in a barplot, for each GO ontology, the number of GO terms showing enrichment at both levels or only at one level.
}
\usage{IdentityPlot(object, outputformat="on screen")}
\arguments{
   \item{object}{an object of class \code{\linkS4class{GOsims}}}
   \item{outputformat}{a character string specifying if the plot is saved in jpeg (\code{jpeg}), postscript (\code{postscript}), pdf (\code{pdf}) format, or it is simply displayed on the screen(\code{on screen}). By default this value is \code{on screen}.}
   }
\author{
	Toma Tebaldi, Erik Dassi, Galena Kostoska
}

\seealso{
\code{\link{GOComparison}}
\code{\linkS4class{GOsims}}
\code{\linkS4class{GOsets}}
\code{\link{GOEnrichment}}
}

\examples{
data(tRanslatomeSampleData)
IdentityPlot(CCComparison, outputformat="on screen")
}

\keyword{GOsims}
\keyword{GOcomparison}
\keyword{GOsets}
\keyword{GOenrichment}

