\name{GOsets-class}
\docType{class}
\alias{GOsets-class}
\alias{Radar,GOsets-method}
\alias{Heatmap,GOsets-method}
\alias{GOComparison,GOsets-method}
\alias{show,GOsets-method}

\title{Class \code{GOsets}}
\description{
A class generated from the function \code{GOEnrichment()} containing the overrepresented Gene Ontology terms from the lists of differentially expressed genes coming from the analysis of two different expression levels. 
}
\section{Slots}{
  \describe{
    \item{\code{enriched.table}:}{Object of class \code{data.frame} containing all the GO enrichment analysis. The data frame contains one row for each GO term, while columns are organized as follows:
	\enumerate{
	\item \code{ontology}, GO ontology of the term, either \code{BP},\code{CC} or \code{MF}.
	\item \code{level}, level of analysis on which the enrichment is calculated. The names of the two levels are taken from the object \code{\linkS4class{DEGs}} given as an input to the function \code{GOEnrichment}. 
	\item \code{GO.ID}, Gene Ontology accession ID of the term.  
	\item \code{term}, Gene Ontology name of the term. 
	\item \code{annotated}, number of genes in the human genome associated to the GO term. 
	\item \code{significant}, number of genes in the list of DEGs associated to the GO term. 
	\item \code{expected}, number of genes in the list of DEGs expected to be associated to the GO term by chance. 
	\item \code{expected}, number of genes in the list of DEGs expected to be associated to the GO term by chance. 
	\item \code{pv.fisher}, enrichment p-value calculated with the Fisher test. The result is dependent on the method chosen to calculate enrichment(\code{classic},\code{elim}, \code{weight},\code{weight01} or \code{parentchild}).
	\item \code{pv.fisher.BH}, the Benjamini-Hochberg corrected enrichment p-value calculated with the Fisher test according to the chosen enrichment method.

}}    
	\item{\code{label.level.enriched}:}{Object of class \code{character} specifying the names of the two levels compared in the experiment.}

}
}
\section{Accessors}{
  \describe{
		\item{enriched.table}{\code{signature(object = "GOsets")}:  displays the slot \code{enriched.table}.}
		\item{label.level.enriched}{\code{signature(object = "GOsets")}: displays the slot \code{label.level.enriched}, slot of class \code{character} specifying the names of the two levels compared in the experiment.}
}
}
\section{Methods}{
  \describe{
    \item{GOComparison}{\code{signature(object = "GOsets")}: returns an object of class \code{\linkS4class{GOsims}}, containing identity and semantic similarity comparisons among the enriched GO terms.}
    \item{Heatmap}{\code{signature(object = "GOsets")}: enables the generation of a heatmap of the top enriched GO terms for the first and second level of analysis.}
    \item{Radar}{\code{signature(object = "GOsets")}: enables the generation of a radar plot of the top enriched GO terms for the first and second level of analysis.}
    \item{show}{\code{signature(object = "GOsets")}: displays all the two slots of the class.}
	 }
}
\references{
Ashburner M et al(2000). Gene ontology: tool for the unification of
biology. Nat. Genet. May 2000;25(1):25-9.

Alexa A, Rahnenfuhrer J, Lengauer T. Improved scoring of functional
groups from gene expression data by decorrelating go graph
structure. Bioinformatics 2006, 22(13):1600-7.
}
\author{
	Toma Tebaldi, Erik Dassi, Galena Kostoska
}

\seealso{
\code{\linkS4class{EnrichedSets}}
\code{\link{GOComparison}}
\code{\linkS4class{DEGs}}
}
\examples{
showClass("GOsets")
}
\keyword{classes}
\keyword{GO}
\keyword{enrichment analysis}
