\name{RegulatoryEnrichment}
\alias{RegulatoryEnrichment}

\title{RegulatoryEnrichmentHelpfile}
\description{
RegulatoryEnrichment is a function which, given as input an object of class \code{\linkS4class{DEGs}}, identifies overrepresented post-transcriptional regulators (RNA-binding proteins, microRNA, etc) controlling differentially expressed genes. The analysis is applied to a dataset of experimentally determined post-transcriptional interactions (i.e. regulator-UTR interaction) extracted from AURA (http://aura.science.unitn.it). Moreover, the function can identify enriched regulators for separate classes of genes of interest: only up-regulated genes, only down-regulated genes or both of them together. The output of the function is an object of class \code{\linkS4class{EnrichedSets}}, containing the results of the enrichment analysis. 
}
\usage{RegulatoryEnrichment(object, classOfDEGs="both",
    significance.threshold = 0.05, mult.cor=TRUE)}
\arguments{
 \item{object}{an object of class \code{\linkS4class{DEGs}}}
 \item{classOfDEGs}{a character string specifying the class of genes for which we want to detect enriched regulators: \code{up} for considering only up-regulated genes, \code{down} for considering only down-regulated genes, \code{both} for considering all DEGs, independently from the direction of their changes. The default is set to \code{both}.}
 \item{significance.threshold}{a numeric value specifying the significance threshold upon which the regulators are considered significantly over-represented. By default it is se to \code{0.05}.}
 \item{mult.cor}{a boolean variable specifying whether the significance threshold is applied to the multiple test corrected or to the original p-values obtained from the selected enrichment method. By default it is set to \code{TRUE}.}
}
\value{
	An object of class \code{\linkS4class{EnrichedSets}}
}
\references{
	Dassi E et al(2012). AURA: Atlas of UTR Regulatory Activity. Bioinformatics. 28(1):142-4. 
}
\author{
	Erik Dassi, Toma Tebaldi
}

\seealso{
\code{\linkS4class{TranslatomeDataset}}
\code{\link{computeDEGs}}
\code{\linkS4class{DEGs}}
}

\examples{
data(tRanslatomeSampleData)
RegulatoryEnrichment(limma.DEGs, significance.threshold = 0.05)
}

\keyword{Regulatory}
\keyword{Enrichment}
\keyword{DEGs}
\keyword{GO}
\keyword{enrichment analysis}