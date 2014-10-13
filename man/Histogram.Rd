\name{Histogram}
\alias{Histogram}

\title{HistogramHelpfile}
\description{
This function shows in a histogram the number of genes showing differential expression at both levels, or up/down regulated only at one level.
}
\usage{Histogram(object, plottype ="summary", outputformat="on screen")}
\arguments{
   \item{object}{an object of class \code{\linkS4class{DEGs}}}
   \item{plottype}{a character string specifying whether the histogram should display a summary of DEGs classes (\code{summary}) or detailed classes taking into account the number of genes up or down regulated in the first or second level (\code{detailed}).}
   \item{outputformat}{a character string specifying if the plot is saved in jpeg (\code{jpeg}), postscript (\code{postscript}), pdf (\code{pdf}) format, or it is simply displayed on the screen(\code{on screen}). By default this value is \code{on screen}.}
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
Histogram(limma.DEGs,plottype ="summary", outputformat="on screen")
}

\keyword{DEGs}
\keyword{histogram}
