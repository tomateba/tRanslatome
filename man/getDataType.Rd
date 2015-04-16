\name{getDataType}
\alias{getDataType}

\title{getDataTypeTranslatomeDatasetHelpfile}
\description{
This function displays an object of class \code{character} specifying the type of data represented by \code{expr.matrix}. By default it is set to \code{array}, the other accepted value is \code{ngs}. It takes as input an object of class \code{\linkS4class{TranslatomeDataset}}.
}
\usage{getDataType(object)}
\arguments{
 \item{object}{an object of class \code{\linkS4class{TranslatomeDataset}}.}
 }

\author{
	Toma Tebaldi, Erik Dassi, Galena Kostoska
}

\seealso{
\code{\linkS4class{TranslatomeDataset}}
\code{\link{computeDEGs}}
}

\examples{
data(tRanslatomeSampleData)
getDataType(translatome.analysis)
}

\keyword{TranslatomeDataset}
\keyword{getDataType}