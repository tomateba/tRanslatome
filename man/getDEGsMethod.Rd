\name{getDEGsMethod}
\alias{getDEGsMethod}

\title{getDEGsMethodDEGsHelpfile}
\description{
This function displays an object of class \code{character} specifying the method that the user employed in the differential expression analysis. It can have one the following values: \code{limma}, \code{t-test}, \code{RP}, \code{SAM}, \code{ANOTA}, \code{DESeq} and \code{none}.By default, this value is set to \code{limma}. It takes as input an object of class \code{\linkS4class{DEGs}}.
}
\usage{getDEGsMethod(object)}
\arguments{
 \item{object}{an object of class \code{\linkS4class{DEGs}}.}
 }

\author{
	Toma Tebaldi, Erik Dassi, Galena Kostoska
}

\seealso{
\code{\linkS4class{DEGs}}
\code{\link{computeDEGs}}
\code{\linkS4class{TranslatomeDataset}}
}

\examples{
data(tRanslatomeSampleData)
getDEGsMethod(limma.DEGs)
}

\keyword{DEGs}
\keyword{getDEGsMethod}