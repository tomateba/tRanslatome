\name{Scatterplot}
\alias{Scatterplot}

\title{ScatterplotHelpfile}
\description{
This plot shows each gene as dot uniquely determined by its log2 fold change at the first level (represented on the x-axis) and the fold change at the second level (represented on the y-axis).
}
\usage{Scatterplot(object, outputformat="on screen",track="")}
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
Scatterplot(limma.DEGs, outputformat="on screen",track="")
}

\keyword{DEGs}
\keyword{scatterplot}
\keyword{fold change}
