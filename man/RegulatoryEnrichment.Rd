\name{RegulatoryEnrichment}
\alias{RegulatoryEnrichment}

\title{RegulatoryEnrichmentHelpfile}
\description{
RegulatoryEnrichment is a function which, given as input an object of class \code{\linkS4class{DEGs}}, identifies overrepresented post-transcriptional regulators (RNA-binding proteins, microRNA, etc) controlling differentially expressed genes. The analysis is by default applied to a dataset of experimentally determined post-transcriptional interactions (i.e. regulator-UTR interaction) extracted from AURA (http://aura.science.unitn.it). However, the user can specify a custom dataset onto which the analysis can be performed (see arguments for details). Moreover, the function can identify enriched regulators for separate classes of genes of interest: only up-regulated genes, only down-regulated genes or both of them together. The method works by exploiting two lists: one containing all genes regulated by each of the post-transcriptional regulators, and the other containing the number of regulated and non-regulated genes for each of these post-transcriptional regulators in the backgroung gene set (usually the whole genome). By means of these two lists it is possible to compute a Fisher enrichment p-value indicating whether a significant group of genes in the DEGs list is likely to be regulated by one or more of these post-transcriptional regulators. The output of the function is an object of class \code{\linkS4class{EnrichedSets}}, containing the results of the enrichment analysis. 
}
\usage{RegulatoryEnrichment(object, classOfDEGs="both",
    significance.threshold = 0.05, mult.cor=TRUE, regulated.identities=NULL, regulated.counts=NULL)}
\arguments{
 \item{object}{an object of class \code{\linkS4class{DEGs}}}
 \item{classOfDEGs}{a character string specifying the class of genes for which we want to detect enriched regulators: \code{up} for considering only up-regulated genes, \code{down} for considering only down-regulated genes, \code{both} for considering all DEGs, independently from the direction of their changes. The default is set to \code{both}.}
 \item{significance.threshold}{a numeric value specifying the significance threshold upon which the regulators are considered significantly over-represented. By default it is se to \code{0.05}.}
 \item{mult.cor}{a boolean variable specifying whether the significance threshold is applied to the multiple test corrected or to the original p-values obtained from the selected enrichment method. By default it is set to \code{TRUE}.}
 \item{regulated.identities}{a matrix containing two columns (RegulatoryElement, RegulatedGenes) specyifing, for each row, a regulatory element name and the comma-separated list of genes it regulates. The user can use the regulatory.elements.regulated table in the tRanslatomeSampleData dataset as a template. By default this argument is NULL, which implies the dataset obtained from AURA will be used.}
 \item{regulated.counts}{a matrix containing three columns (RegulatoryElement, RegulatedGenes, NonRegulatedGenes) specyifing, for each row, a regulatory element name, the number of genes it regulates in the background gene set and the number of genes it does not regulate in the background gene set. The user can use the regulatory.elements.counts table in the tRanslatomeSampleData dataset as a template. By default this argument is NULL, which implies the dataset obtained from AURA will be used.}
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