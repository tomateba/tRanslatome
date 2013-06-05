\name{getConditionLabels}
\alias{getConditionLabels}

\title{getConditionLabelsTranslatomeDatasetHelpfile}
\description{
This function displays an object of class \code{character} specifying the names given to the two conditions. By default, the vector is set to \code{c("control","treated")}, but the user can specify other names. It takes as input an object of class \code{\linkS4class{TranslatomeDataset}}.
}
\usage{getConditionLabels(object)}
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
getConditionLabels(translatome.analysis)
}

\keyword{TranslatomeDataset}
\keyword{getConditionLabels}