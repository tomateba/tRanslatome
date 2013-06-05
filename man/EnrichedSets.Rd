\name{EnrichedSets-class}
\docType{class}
\alias{EnrichedSets-class}
\alias{enriched.table,EnrichedSets-method}
\alias{label.level.enriched,EnrichedSets-method}
\alias{Radar,EnrichedSets-method}
\alias{Heatmap,EnrichedSets-method}
\alias{show,EnrichedSets-method}

\title{Class \code{EnrichedSets}}
\description{
A class generated from the function \code{RegulatoryEnrichment()} containing the overrepresented post-transcriptional regulatory factors (RBPs, miRNAs, etc) from the lists of differentially expressed genes coming from the analysis of two different expression levels. 
}
\section{Slots}{
  \describe{
    \item{\code{enriched.table}:}{Object of class \code{data.frame} containing all the regulatory enrichment analysis. The data frame contains one row for each regulatory factor, while columns are organized as follows:
	\enumerate{
		\item \code{level}, level of analysis on which the enrichment is calculated. The names of the two levels are taken from the object \code{\linkS4class{DEGs}} given as an input to the function \code{GOEnrichment}. 
		\item \code{ID}, ID of the regulatory factor.  
		\item \code{number}, number of genes in the human genome associated to the regulatory factor. 
		\item \code{list}, list of genes in the list of DEGs associated to the regulatory factor. 
		\item \code{pv.fisher}, enrichment p-value calculated with the Fisher test.
		\item \code{pv.fisher.BH}, the Benjamini-Hochberg corrected enrichment p-value calculated with the Fisher test according to the chosen enrichment method.
	}}    
	\item{\code{label.level.enriched}:}{Object of class \code{character} specifying the names of the two levels compared in the experiment.}
}
}
\section{Accessors}{
  \describe{
		\item{enriched.table}{\code{signature(object = "EnrichedSets")}:  displays the slot \code{enriched.table}.}
		\item{label.level.enriched}{\code{signature(object = "EnrichedSets")}: displays the slot label.level.enriched, slot of class \code{character} specifying the names of the two levels compared in the experiment.}
}
}
\section{Methods}{
  \describe{
    \item{Heatmap}{\code{signature(object = "EnrichedSets")}: enables the generation of a heatmap of the top enriched regulatory factors for the first and second level of analysis.}
    \item{label.level.enriched}{\code{signature(object = "EnrichedSets")}: displays the slot \code{label.level.enriched}.}
    \item{Radar}{\code{signature(object = "EnrichedSets")}: enables the generation of a radar plot of the top enriched regulatory factors for the first and second level of analysis.}
    \item{show}{\code{signature(object = "EnrichedSets")}: displays all the two slots of the class.}
	 }
}
\references{
	Alexa A, Rahnenfuhrer J, Lengauer T. Improved scoring of functional groups from gene expression data by decorrelating go graph structure. Bioinformatics 2006, 22(13):1600-7.
	Dassi E et al(2012). AURA: Atlas of UTR Regulatory Activity. Bioinformatics. 28(1):142-4. 
}
\author{
	Erik Dassi, Toma Tebaldi 
}

\seealso{
\code{\linkS4class{GOsets}}
\code{\linkS4class{DEGs}}
}
\examples{
showClass("EnrichedSets")
}
\keyword{classes}
\keyword{GO}
\keyword{enrichment analysis}
