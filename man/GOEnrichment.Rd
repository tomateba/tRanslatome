\name{GOEnrichment}
\alias{GOEnrichment}

\title{GOEnrichmentHelpfile}
\description{
GOEnrichment is a function which, given as input an object of class \code{\linkS4class{DEGs}}, identifies overrepresented GO terms among differentially expressed genes. The analysis can be applied to all the GO ontologies or restricted to GO terms specifically belonging to one ontology: molecular function, cellular component or biological process. Moreover the function can identify enriched GO terms for separate classes of genes of interest: only up-regulated genes, only down-regulated genes or both of them together. The output of the function is an object of class \code{\linkS4class{GOsets}}, containing the results of the enrichment analysis. 
}
\usage{GOEnrichment(object, ontology="all", classOfDEGs="both",
    test.method="classic", test.threshold = 0.05, mult.cor=TRUE)}
\arguments{
 \item{object}{an object of class \code{\linkS4class{DEGs}}}
 \item{ontology}{a character string specifying the GO ontology of interest: \code{CC} for Cellular Component, \code{BP} for Biological Process, \code{MF} for Molecular Function or \code{all} for all the three ontologies. The default is set to \code{all}.}
 \item{classOfDEGs}{a character string specifying the class of genes for which we want to detect enriched GO terms: \code{up} for considering only up-regulated genes, \code{down} for considering only down-regulated genes, \code{both} for considering all DEGs, independently from the direction of their changes. The default is set to \code{both}.}
 \item{test.method}{a character string specifying the statistical method to calculate the enrichment. By default it is set to \code{classic} (enrichment is measured with the classic Fisher exact test), but it can also be set to \code{elim}, \code{weight}, \code{weight01} or \code{parentchild}. All these methods are implemented in the \code{topGO} Bioconductor package}
 \item{test.threshold}{a numeric value specifying the significance threshold upon which the GO terms are considered significantly over-represented. By default it is se to \code{0.05}.}
 \item{mult.cor}{a boolean variable specifying whether the significance threshold is applied to the multiple test corrected or to the original p-values obtained from the selected enrichment method. By default it is set to \code{TRUE}.}
}
\value{
An object of class \code{\linkS4class{GOsets}}
}
\references{
Ashburner M et al(2000). Gene ontology: tool for the unification of
biology. Nat. Genet. May 2000;25(1):25-9.

Alexa A, Rahnenfuhrer J, Lengauer T. Improved scoring of functional
groups from gene expression data by decorrelating go graph
structure. Bioinformatics 2006, 22(13):1600-7
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
GOEnrichment(limma.DEGs,ontology="CC",classOfDEGs="up",
    test.method="classic", test.threshold = 0.05,mult.cor = TRUE)
}

\keyword{DEGs}
\keyword{GO}
\keyword{enrichment analysis}