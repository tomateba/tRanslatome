\name{CVplot}
\alias{CVplot}

\title{CVplotDEGsHelpfile}
\description{
The CV plot displays the coefficients of variation of the genes on the x-axis and the log2 fold changes on the y-axis. The upper plot represents a CV plot for the first level of gene expression analysis, whereas the other one refers to the second level. DEGs are color labeled. This function takes as input an object of class \code{\linkS4class{DEGs}}.
}
\usage{CVplot(object, outputformat="on screen",track="")}
\arguments{
 \item{object}{an object of class \code{\linkS4class{DEGs}}.}
 \item{outputformat}{a character string specifying if the plot is saved in jpeg (\code{jpeg}), postscript (\code{postscript}), pdf (\code{pdf}) format, or it is simply displayed on the screen(\code{on screen}). By default this value is \code{on screen}.}
 \item{track}{a character vector of gene names that will be explicitly highlighted in the scatterplot, if they match any gene contained in the object of class \code{\linkS4class{DEGs}}. By default this vector is empty.}
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
CVplot(limma.DEGs, outputformat="on screen",track=c("BLM"))
}

\keyword{DEGs}
\keyword{CV}
\keyword{coefficient of variation}
