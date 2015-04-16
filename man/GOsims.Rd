\name{GOsims-class}
\docType{class}
\alias{GOsims-class}
\alias{identity.matrix,GOsims-method}
\alias{similarity.matrix,GOsims-method}
\alias{average.similarity.scores,GOsims-method}
\alias{IdentityPlot,GOsims-method}
\alias{SimilarityPlot,GOsims-method}

\title{Class \code{GOsims}}
\description{
A class generated from the function \code{GOComparison()} containing the result of the differential expression analysis at the two expression levels. 
}
\section{Slots}{
 \describe{
 \item{\code{similarity.matrix}:}{Object of class \code{matrix} containing the semantic similarity comparison between the GO terms enriched at the two levels of analysis. The matrix contains one row for each pairwise comparison between GO terms enriched in the two levels, while columns are organized as follows:
	\enumerate{
	\item \code{direction}, direction of the comparison: "first level to second level" means that an enriched GO term from the first level is compared to the most similar GO term enriched in the second level, "second level to first level" means that an enriched GO term from the second level is compared to the most similar GO term enriched in the first level. Level names are taken from the object of class \code{\linkS4class{GOsets}} given as input to the function \code{GOcomparison}.
	\item \code{ontology}, GO ontology of the two terms, either \code{BP},\code{CC} or \code{MF}.
	\item \code{level}, gene expression levels for which the term is enriched, either "first level only", "second level only" or "both levels". Level names are taken from the object of class \code{\linkS4class{GOsets}} given as input to the function \code{GOcomparison}.
	\item \code{start.GO.ID}, Gene Ontology accession ID of the first term of the comparison. 
	\item \code{start.term}, Gene Ontology name of the first term of the comparison.	
	\item \code{end.GO.ID}, Gene Ontology accession ID of the second term of the comparison. 
	\item \code{end.term}, Gene Ontology name of the second term of the comparison.	
	\item \code{similarity.score}, semantic similarity value between the two campared GO terms.}
	}
 \item{\code{identity.matrix}:}{Object of class \code{matrix} containing the identity comparison between the GO terms enriched at the two levels of analysis. The matrix contains one row for each GO term, while columns are organized as follows:
	\enumerate{
	\item \code{ontology}, GO ontology of the term, either \code{BP},\code{CC} or \code{MF}.
	\item \code{level}, gene expression levels for which the term is enriched, either "first level only", "second level only" or "both levels". Level names are taken from the object of class \code{\linkS4class{GOsets}} given as input to the function \code{GOcomparison}.
	\item \code{GO.ID}, Gene Ontology accession ID of the term. 
	\item \code{term}, Gene Ontology name of the term. }
	}
 \item{\code{average.similarity.scores}:}{Object of class \code{vector} containing the general semantic similarity scores between the GO terms enriched at the two levels of analysis. One similarity score, ranging from 0 to 1, is produced for each GO ontology having at least one enriched term.}
 }
}
\section{Accessors}{
 \describe{
	\item{identity.matrix}{\code{signature(object = "GOsims")}: displays the slot \code{identity.matrix}.}
	\item{similarity.matrix}{\code{signature(object = "GOsims")}: displays the slot \code{similarity.matrix}.}
	\item{average.similarity.scores}{\code{signature(object = "GOsims")}: displays the slot \code{average.similarity.scores}.}}
}
\section{Methods}{
 \describe{
	\item{IdentityPlot}{\code{signature(object = "GOsims")}: enables the generation of a barplot where, for each GO ontology, the number of GO terms enriched at different levels are displayed.}
	\item{SimilarityPlot}{\code{signature(object = "GOsims")}: enables the generation of a barplot where, for each GO ontology, the semantic similarity value between GO terms enriched at different levels is displayed.}
}
}
\author{
	Toma Tebaldi, Erik Dassi, Galena Kostoska
}

\seealso{
\code{\linkS4class{TranslatomeDataset}}
\code{\link{computeDEGs}}
\code{\link{GOEnrichment}}
\code{\linkS4class{GOsets}}
}

\examples{
showClass("GOsims")
}
\keyword{classes}
\keyword{GO}
\keyword{enrichment analysis}