\name{Heatmap}
\alias{Heatmap}

\title{HeatmapEnrichedHelpfile}
\description{
Heatmap is a function that plots the top enriched regulatory element or GO terms for the first and second level of analysis in a heatmap. The number of terms to be displayed can be set.
}
\usage{Heatmap(object, outputformat="on screen", n.nodes.1stlevel="5",
    n.nodes.2ndlevel="5", mult.cor=TRUE, ...)}
\arguments{
  \item{object}{an object of class \code{\linkS4class{GOsets}}.}
  \item{outputformat}{a character string specifying if the plot is saved in jpeg (\code{jpeg}), postscript (\code{postscript}), pdf (\code{pdf}) format, or it is simply displayed on the screen(\code{on screen}). By default this value is \code{on screen}.}
  \item{n.nodes.1stlevel}{a numeric value specifying the number of top enriched GO terms, from the first level, that will be represented on the plot. By default the value is set to \code{5}.}
  \item{n.nodes.2ndlevel}{a numeric value specifying the number of top enriched GO terms, from the second level, that will be represented on the plot. By default the value is set to \code{5}.}
  \item{mult.cor}{a boolean variable specifying whether the displayed significance values are multiple test corrected or the original p-values obtained from the selected enrichment method. By default it is set to \code{TRUE}.}
  \item{...}{if the class of object is \code{GOsets}, it accepts a character string named \code{ontology} selecting the GO ontology of interest, either \code{CC},\code{BP} or \code{MF}.}
}
\author{
	Toma Tebaldi, Erik Dassi, Galena Kostoska
}

\seealso{
\code{\linkS4class{EnrichedSets}}
\code{\link{RegulatoryEnrichment}}
\code{\linkS4class{GOsets}}
\code{\link{GOEnrichment}}
\code{\linkS4class{GOsims}}
}

\examples{
data(tRanslatomeSampleData)
Heatmap(CCEnrichment, ontology="MF")
}

\keyword{DEGs}
\keyword{GO}
\keyword{Regulatory}
\keyword{enrichment analysis}
\keyword{heatmap}