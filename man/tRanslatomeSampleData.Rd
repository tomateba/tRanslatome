\name{tRanslatomeSampleData}
\docType{data}
\alias{tRanslatomeSampleData}
\alias{translatome.analysis}
\alias{CCEnrichment}
\alias{limma.DEGs}
\alias{tRanslatomeData}
\alias{expressionMatrix}
\alias{CCComparison}
\alias{regulatory.elements.counts}
\alias{regulatory.elements.regulated}
\title{Sample data set for tRanslatome}
\description{
  matrix with data coming from differentiated and undifferentiated human HepaRG cells 
}

\usage{data(tRanslatomeSampleData)}


\format{tRanslatomeSampleData is a list that has 7 components : 

expressionMatrix is a matrix of 12 columns and 1000 rows containing the microarray signals used as input of the function GetDEGs. Columns are organized in biological triplicates as follows: 
  tot.diff.a, tot.diff.b, tot.diff.c contain the signals coming from the
  transcriptome of the differentiated cell line.

  
  tot.undiff.a, tot.undiff.b, tot.undiff.c contain the signals coming
  from the transcriptome of the undifferentiated cell line.

  
  pol.diff.a, pol.diff.b, pol.diff.c contain the signals coming from the
  translatome of the differentiated cell line.

  
  pol.undiff.a, pol.undiff.b, pol.undiff.c contain the signals coming
  from the translatome of the undifferentiated cell line.
  
  translatome.analysis is an object of class
  \code{\linkS4class{TranslatomeDataset}} generated calling
  \code{newTranslatomeDataset()} on \code{expressionMatrix} (see the
  examples section for the exactcall).
  
 limma.DEGs is an object of class \code{\linkS4class{DEGs}} generated
 calling \code{getDEGs()} on \code{expressionMatrix} (see the examples
 section for the exact call). 
 
 CCEnrichment is an object of class \code{\linkS4class{GOsets}}
 generated calling \code{GOEnrichment()} on \code{limma.DEGs} (see the
 examples section for the exact call).
 
 CCComparison is an object of class \code{\linkS4class{GOsims}}
 generated calling \code{GOComparison()} on \code{CCEnrichment} (see the
 examples section for the exact call).  
 
 regulatory.elements.counts is an object of class \code{data.frame}
 containing the background numbers of regulated and non-regulated genes
 for each post-transcriptional regulatory factor considered by function
 \code{RegulatoryEnrichment}.
 
 regulatory.elements.regulated is an object of class \code{data.frame} containing the list of regulated genes for each post-transcriptional regulatory factor considered by function \code{RegulatoryEnrichment}
 }

\source{
Parent R, Kolippakkam D, Booth G, Beretta L. Mammalian target of rapamycin activation impairs hepatocytic differentiation and targets genes moderating lipid homeostasis and hepatocellular growth. Cancer Res. 2007;67(9):4337-4345}

\examples{
 ##load the tRanslatome sample data
 data(tRanslatomeSampleData) 
}

\seealso{
 \code{\linkS4class{TranslatomeDataset}}
 \code{\link{computeDEGs}}
 \code{\linkS4class{DEGs}}
 \code{\link{RegulatoryEnrichment}}
 \code{\linkS4class{EnrichedSets}}
 \code{\link{GOEnrichment}}
 \code{\linkS4class{GOsets}}
 \code{\link{GOComparison}}
 \code{\linkS4class{GOsims}}
}

\keyword{datasets}
