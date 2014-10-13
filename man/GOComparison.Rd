\name{GOComparison}
\alias{GOComparison}

\title{GOComparisonHelpfile}
\description{
GOComparison is a function which takes as input an object of class \code{\linkS4class{GOsets}}, containing the results of a GO enrichment analysis applied to both levels, and returns an object of class \code{\linkS4class{GOsims}}, containing a variety of comparisons among the enriched GO terms.
}
\usage{GOComparison(object)}
\arguments{
 \item{object}{an object of class \code{\linkS4class{GOsets}}}
}
\value{
An object of class \code{\linkS4class{GOsims}}}

\references{
Wang et al.(2007) A new method to measure the semantic similarity of go terms Bioinformatics , 23:0 1274-81, May 2007.
}
\author{
	Toma Tebaldi, Erik Dassi, Galena Kostoska
}

\seealso{
\code{\link{GOEnrichment}}
\code{\linkS4class{GOsets}}
\code{\linkS4class{GOsims}}
}

\examples{
data(tRanslatomeSampleData)
GOComparison(CCEnrichment)
}

\keyword{DEGs}
\keyword{GO}
\keyword{enrichment analysis}
\keyword{Fisher test}