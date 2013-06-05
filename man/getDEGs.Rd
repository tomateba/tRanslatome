\name{getDEGs}
\alias{getDEGs}

\title{getDEGsTranslatomeDatasetHelpfile}
\description{
This function displays an object of class \code{DEGs} in which each gene is assigned an expression class: up- or down-regulated at the first level, up- or down-regulated at the second level, up-regulated at both levels, down-regulated at both levels, up-regulated at the first level and down-regulated at the second level and vice versa. It takes as input an object of class \code{\linkS4class{TranslatomeDataset}}.
}
\usage{getDEGs(object)}
\arguments{
 \item{object}{an object of class \code{\linkS4class{TranslatomeDataset}}.}
 }

\author{
	Toma Tebaldi, Erik Dassi, Galena Kostoska
}

\seealso{
\code{\linkS4class{TranslatomeDataset}}
\code{\link{computeDEGs}}
\code{\linkS4class{DEGs}}
}

\examples{
data(tRanslatomeSampleData)
getDEGs(translatome.analysis)
}

\keyword{TranslatomeDataset}
\keyword{getDEGs}