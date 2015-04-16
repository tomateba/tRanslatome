\name{TranslatomeDataset-class}
\docType{class}
\alias{TranslatomeDataset-class}
\alias{getExprMatrix,TranslatomeDataset-method}
\alias{getConditionA,TranslatomeDataset-method}
\alias{getConditionB,TranslatomeDataset-method}
\alias{getConditionC,TranslatomeDataset-method}
\alias{getConditionD,TranslatomeDataset-method}
\alias{getConditionLabels,TranslatomeDataset-method}
\alias{getDataType,TranslatomeDataset-method}
\alias{getLevelLabels,TranslatomeDataset-method}
\alias{getDEGs,TranslatomeDataset-method}
\alias{computeDEGs,TranslatomeDataset-method}
\alias{show,TranslatomeDataset-method}

\title{Class \code{TranslatomeDataset}}
\description{
A class generated from the function \code{newTranslatomeDataset()} containing the input matrix, the condition vectors and labels and the result of the differential expression analysis at the two expression levels. This class represents an analysis in its entirety, containing all data from input parameters to output results. 
}
\section{Slots}{
 \describe{
	 \item{\code{expr.matrix}:}{Object of class \code{matrix} specifying the normalized signal intensity data, each row representing a gene and each column representing a sample.}
	 \item{\code{cond.a}:}{Object of class \code{character} specifying a vector of column names belonging to expression matrix. These columns contain the signal intensity data coming from the samples of the first expression level of the control condition.}
	 \item{\code{cond.b}:}{Object of class \code{character} specifying a vector of column names belonging to expression matrix. These columns contain the signal intensity data coming from the samples of the first expression level of the treatment condition.}
	 \item{\code{cond.c}:}{Object of class \code{character} specifying a vector of column names belonging to expression matrix. These columns contain the signal intensity data coming from the samples of the second expression level of the control condition.}
	 \item{\code{cond.d}:}{Object of class \code{character} specifying a vector of column names belonging to expression matrix. These columns contain the signal intensity data coming from the samples of the second expression level of the treatment condition.}
	 \item{\code{data.type}:}{Object  of class \code{character} specifying the type of the data contained in \code{exprMatrix}. The default is "array" and the alternative is "ngs".}
	 \item{\code{label.condition}:}{Object of class \code{character} specifying the names given to the two conditions. By default, these values are "control" and "treated", but user can specify others.}
	 \item{\code{label.level}:}{Object of class \code{character} specifying the names given to the two levels. By default levels are named "1st level" and "2nd level", but the user can specify others.}
	 \item{\code{DEGs}:}{Object of class \code{DEGs} specifying the result of the differential expression analysis at the two expression levels.}
 }
}
\section{Accessors}{
 \describe{
	 \item{getExprMatrix}{\code{signature(object = "TranslatomeDataset")}: displays anObject of class \code{matrix} specifying the normalized signal intensity data, each row representing a gene and each column representing a sample.}
	 \item{getConditionA}{\code{signature(object = "TranslatomeDataset")}: displays anObject of class \code{character} specifying a vector of column names belonging to expression matrix. These columns contain the signal intensity data coming from the samples of the first expression level of the control condition.}
	 \item{getConditionB}{\code{signature(object = "TranslatomeDataset")}: displays anObject of class \code{character} specifying a vector of column names belonging to expression matrix. These columns contain the signal intensity data coming from the samples of the first expression level of the treatment condition.}
	 \item{getConditionC}{\code{signature(object = "TranslatomeDataset")}: displays anObject of class \code{character} specifying a vector of column names belonging to expression matrix. These columns contain the signal intensity data coming from the samples of the second expression level of the control condition.}
	 \item{getConditionD}{\code{signature(object = "TranslatomeDataset")}: displays anObject of class \code{character} specifying a vector of column names belonging to expression matrix. These columns contain the signal intensity data coming from the samples of the second expression level of the treatment condition.}
	 \item{getDataType}{\code{signature(object = "TranslatomeDataset")}: displays anObject  of class \code{character} specifying the type of the data contained in \code{exprMatrix}. The default is "array" and the alternative is "ngs".}
	 \item{getConditionLabels}{\code{signature(object = "TranslatomeDataset")}: displays anObject of class \code{character} specifying the names given to the two conditions. By default, these values are "control" and "treated", but user can specify others.}
	 \item{getLevelLabels}{\code{signature(object = "TranslatomeDataset")}: displays anObject of class \code{character} specifying the names given to the two levels. By default levels are named "1st level" and "2nd level", but the user can specify others.}
	 \item{getDEGs}{\code{signature(object = "TranslatomeDataset")}: displays anObject of class \code{DEGs} specifying the result of the differential expression analysis at the two expression levels.}
 }
}
\section{Methods}{
 \describe{
	\item{computeDEGs}{\code{signature(object = "TranslatomeDataset")}: compute the differentially expressed genes at the two levels by means of the chosen method, returning anObject of class \code{DEGs}.}
	\item{show}{\code{signature(object = "TranslatomeDataset")}: displays all the slots of the class.}
 }
}
\author{
	Toma Tebaldi, Erik Dassi, Galena Kostoska
}

\seealso{
\code{\linkS4class{TranslatomeDataset}}
\code{\link{newTranslatomeDataset}}
\code{\link{computeDEGs}}
}

\examples{
showClass("TranslatomeDataset")
}

\keyword{classes}
\keyword{TranslatomeDataset}
\keyword{differential expression}