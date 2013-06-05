\name{DEGs-class}
\docType{class}
\alias{DEGs-class}
\alias{CVplot,DEGs-method}
\alias{DEGs.table,DEGs-method}
\alias{FC.threshold,DEGs-method}
\alias{RegulatoryEnrichment,DEGs-method}
\alias{GOEnrichment,DEGs-method}
\alias{Histogram,DEGs-method}
\alias{label.condition,DEGs-method}
\alias{label.level.DEGs,DEGs-method}
\alias{MAplot,DEGs-method}
\alias{getDEGsMethod,DEGs-method}
\alias{Scatterplot,DEGs-method}
\alias{SDplot,DEGs-method}
\alias{show,DEGs-method}
\alias{significance.threshold,DEGs-method}

\title{Class \code{DEGs}}
\description{
A class generated from the function \code{computeDEGs()} containing the result of the differential expression analysis at the two expression levels. 
}
\section{Slots}{
 \describe{
 \item{\code{DEGs.table}:}{Object of class \code{matrix} containing all the differential expression analysis. The matrix contains one row for each gene, while columns are organized as follows:
	\enumerate{
	\item \code{log2 FC(1st level)}, log2 fold change in the first level,
	\item \code{avg(1st level,control)}, average signal intensities of the control samples in the first level,
	\item \code{sd(1st level,control)}, standard deviation of the signal intensities of the control samples in the first level,
	\item \code{avg(1st level,treated)}, average signal intensities of the treated samples in the first level,
	\item \code{sd(1st level,treated)}, standard deviation of the signal intensities of the treated samples in the first level,
	\item \code{pv.METHOD(1st level)}, the differential expression p-value returned from the chosen method in the first level (the name of this column will change according to the selected method),
	\item \code{pv.METHOD.BH(1st level)}, the Benjamini-Hochberg corrected differential expression p-value returned from the chosen method in the first level (the name of this column will change according to the selected method),
	\item \code{log2 FC(2nd level)}, log2 fold change in the second level,
	\item \code{avg(2nd level,control)}, average signal intensities of the control samples in the second level,
	\item \code{sd(2nd level,control)}, standard deviation of the signal intensities of the control samples in the second level,
	\item \code{avg(2nd level,treated)}, average signal intensities of the treated samples in the second level,
	\item \code{sd(2nd level,treated)}, standard deviation of the signal intensities of the treated samples in the second level,
	\item \code{pv.METHOD(2nd level)}, the differential expression p-value returned from the chosen method in the second level (the name of this column will change according to the selected method),
	\item \code{pv.METHOD.BH(2nd level)}, the Benjamini-Hochberg corrected differential expression p-value returned from the chosen method in the second level (the name of this column will change according to the selected method), 
	\item \code{level1Up}, boolean value set to 1 if the gene is upregulated in the first level, otherwise set to 0, 
	\item \code{level1Down}, boolean value set to 1 if the gene is downregulated in the first level, otherwise set to 0, 
	\item \code{level2Up}, boolean value set to 1 if the gene is upregulated in the second level, otherwise set to 0, 
	\item \code{level2Down}, boolean value set to 1 if the gene is downregulated in the second level, otherwise set to 0, 
	\item \code{DownDown}, boolean value set to 1 if the gene is downregulated in both levels, otherwise set to 0, 
	\item \code{DownUp}, boolean value set to 1 if the gene is downregulated in the first level and up regulated in the second level, otherwise set to 0, 
	\item \code{UpDown}, boolean value set to 1 if the gene is upregulated in the first level and downregulated in the second level, otherwise set to 0, 
	\item \code{UpUp}, boolean value set to 1 if the gene is upregulated in both levels, otherwise set to 0}
	}	
 \item{\code{method}:}{Object of class \code{character} specifying the statistical method used to detect DEGs.}
 \item{\code{significance.threshold}:}{Object of class \code{numeric} specifying the significance threshold used to detect DEGs.}
 \item{\code{FC.threshold}:}{Object of class \code{numeric} specifying the fold change threshold used to detect DEGs.}
 \item{\code{label.level.DEGs}:}{Object of class \code{character} specifying the names of the two levels compared in the experiment.}
 \item{\code{label.condition}:}{Object of class \code{character} specifying the names of the two conditions compared in the experiment.}
 }
}
\section{Accessors}{
 \describe{
 \item{getDEGsMethod}{\code{signature(object = "DEGs")}: displays an object of class \code{character} specifying the statistical method used to detect DEGs.}
 \item{significance.threshold}{\code{signature(object = "DEGs")}: displays an object of class \code{numeric} specifying the significance threshold used to detect DEGs.}
 \item{FC.threshold}{\code{signature(object = "DEGs")}: displays an object of class \code{numeric} specifying the fold change threshold used to detect DEGs.}
 \item{label.level.DEGs}{\code{signature(object = "DEGs")}: displays an object  of class \code{character} specifying the names of the two levels compared in the experiment.}
 \item{label.condition}{\code{signature(object = "DEGs")}: displays an object of class \code{character} specifying the names of the two conditions compared in the experiment.}
 \item{DEGs.table}{\code{signature(object = "DEGs")}: displays an object of class \code{matrix} specifying the table of computed DEGs.}
 }
}
\section{Methods}{
 \describe{
	 \item{CVplot}{\code{signature(object = "DEGs")}: enables the generation of a plot where, for each level, signal coefficients of variation for each gene are displayed against their log2 fold changes. DEGs are color labeled.}
	 \item{GOEnrichment}{\code{signature(object = "DEGs")}: enables the GO enrichment analysis of the differentially expressed genes at each level.}
	 \item{RegulatoryEnrichment}{\code{signature(object = "DEGs")}: enables the enrichment analysis of post-transcriptional regulators (RBPs, miRNAs, ecc) in the DEGs gene list, by means of Fisher test.}
	 \item{Histogram}{\code{signature(object = "DEGs")}: enables the generation of a histogram displaying the number of DEGs belonging to different classes according to their expression behaviour at the two levels.}
	 \item{MAplot}{\code{signature(object = "DEGs")}: enables the generation of a plot where, for each level, average signal intensities for each gene are displayed against their log2 fold changes. DEGs are color labeled.}
	 \item{Scatterplot}{\code{signature(object = "DEGs")}: enables the generation of a plot where fold changes at the first level are displayed for each gene against fold changes at the second level. DEGs are color labeled.}
	 \item{SDplot}{\code{signature(object = "DEGs")}: enables the generation of a plot where, for each level, signal standard deviations for each gene are displayed against their log2 fold changes. DEGs are color labeled.}
	 \item{show}{\code{signature(object = "DEGs")}: displays all the six slots of the class.}
	}
}
\author{
	Toma Tebaldi, Erik Dassi, Galena Kostoska
}

\seealso{
\code{\linkS4class{TranslatomeDataset}}
\code{\link{computeDEGs}}
\code{\link{RegulatoryEnrichment}}
\code{\link{GOEnrichment}}
\code{\linkS4class{EnrichedSets}}
\code{\linkS4class{GOsets}}
}

\examples{
showClass("DEGs")
}

\keyword{classes}
\keyword{DEGs}
\keyword{differential expression}