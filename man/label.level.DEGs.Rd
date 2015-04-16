\name{label.level.DEGs}
\alias{label.level.DEGs}

\title{labelLevelDEGsHelpfile}
\description{
This function displays an object of class \code{character} specifying the names of the two levels compared in the experiment. It takes as input an object of class \code{\linkS4class{DEGs}}.
}
\usage{label.level.DEGs(object)}
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
label.level.DEGs(limma.DEGs)
}

\keyword{DEGs}
\keyword{label.level.DEGs}