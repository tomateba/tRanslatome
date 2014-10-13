\name{label.condition}
\alias{label.condition}

\title{labelConditionDEGsHelpfile}
\description{
This function displays an object of class \code{numeric} specifying the significance threshold used to detect DEGs. It takes as input an object of class \code{\linkS4class{DEGs}}.
}
\usage{label.condition(object)}
\arguments{
 \item{object}{an object of class \code{\linkS4class{DEGs}}.}
 }

\author{
	Toma Tebaldi, Erik Dassi, Galena Kostoska
}

\seealso{
\code{\linkS4class{DEGs}}
\code{\linkS4class{TranslatomeDataset}}
\code{\link{computeDEGs}}
}

\examples{
data(tRanslatomeSampleData)
label.condition(limma.DEGs)
}

\keyword{DEGs}
\keyword{label.condition}