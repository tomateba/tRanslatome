\name{getExprMatrix}
\alias{getExprMatrix}

\title{getExprMatrixTranslatomeDatasetHelpfile}
\description{
This function displays a matrix that contains the normalized signal intensity data, each row representing a gene and each column representing a sample. Row names should correspond to gene names, column names should correspond to sample names. It takes as input an object of class \code{\linkS4class{TranslatomeDataset}}.
}
\usage{getExprMatrix(object)}
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
getExprMatrix(translatome.analysis)
}

\keyword{TranslatomeDataset}
\keyword{getExprMatrix}