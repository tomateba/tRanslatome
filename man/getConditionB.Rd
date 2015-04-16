\name{getConditionB}
\alias{getConditionB}

\title{getConditionBTranslatomeDatasetHelpfile}
\description{
This function displays a character vector of column names (or a numeric vector of columns) belonging to expr.matrix. These columns contain the signal intensity data coming from the samples of the first level of the treatment condition. It takes as input an object of class \code{\linkS4class{TranslatomeDataset}}.
}
\usage{getConditionB(object)}
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
getConditionB(translatome.analysis)
}

\keyword{TranslatomeDataset}
\keyword{getConditionB}