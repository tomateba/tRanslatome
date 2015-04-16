\name{newTranslatomeDataset}
\alias{newTranslatomeDataset}

\title{newTranslatomeDatasetHelpfile}
\description{
This function takes as input a normalized data matrix coming from an high-throughput experiment, along with vectors of column names (or numeric vectors of columns) defining the samples for each of the conditions. It takes as an input a character label specifying the data type(array, ngs) and returns an object of the class \code{\linkS4class{TranslatomeDataset}}, which stores the matrix and the samples definitions and allow to call the \code{computeDEGs} function on it to compute differentially expressed genes in the various conditions.
}
\usage{newTranslatomeDataset(expr.matrix, cond.a, cond.b, cond.c, cond.d,
    data.type="array", label.level=c("1st level","2nd level"),
    label.condition=c("control","treated"))}
\arguments{
 \item{expr.matrix}{a matrix that contains the normalized signal intensity data, each row representing a gene and each column representing a sample. Row names should correspond to gene names, column names should correspond to sample names,}
 \item{cond.a}{a character vector of column names (or a numeric vector of columns) belonging to expr.matrix. These columns contain the signal intensity data coming from the samples of the first level of the control condition (in our example: total RNA, undifferentiated cells),}
 \item{cond.b}{a character vector of column names (or a numeric vector of columns) belonging to expr.matrix. These columns contain the signal intensity data coming from the samples of the first level of the treatment condition(in our example: total RNA, differentiated cells),}
 \item{cond.c}{a character vector of column names (or a numeric vector of columns) belonging to expr.matrix. These columns contain the signal intensity data coming from the samples of the second level of the control condition (in our example: polysomal RNA, undifferentiated cells),}
 \item{cond.d}{a character vector of column names (or a numeric vector of columns) belonging to expr.matrix. These columns contain the signal intensity data coming from the samples of the second level of the treatment group (in our example: polysomal RNA, differentiated cells),}
 \item{data.type}{a character specifying the type of data represented by \code{expr.matrix}. By default it is set to \code{array}, the other accepted value is \code{ngs},}
 \item{label.level}{a character vector specifying the names given to the two levels. By default the vector is set to \code{c("1st level","2nd level")}, but the user can specify other names: in our example the two levels are named "transcriptome" and "translatome",}
 \item{label.condition}{a character vector specifying the names given to the two conditions. By default, the vector is set to \code{c("control","treated")}, but the user can specify other names: in our example the two levels are named "undifferentiated" and "differentiated",}
}
\details{Signals contained in expr.matrix should be previously normalized with standard methods (quantile, percentile shift, ... ) when data is coming from microarrays or in the appropriate cases when it is coming from sequencing experiments.}
\value{
 An object of class \code{\linkS4class{TranslatomeDataset}}
}

\author{
	Toma Tebaldi, Erik Dassi, Galena Kostoska
}

\seealso{
\code{\linkS4class{TranslatomeDataset}}
\code{\link{computeDEGs}}
\code{\linkS4class{DEGs}}
\code{\link{Scatterplot}}
\code{\link{Histogram}}
\code{\link{CVplot}}
\code{\link{MAplot}}
\code{\link{SDplot}}
}

\examples{
data(tRanslatomeSampleData)
translatome.analysis <-
    newTranslatomeDataset(expressionMatrix,
        c("tot.undiff.a","tot.undiff.b","tot.undiff.c"),
        c("tot.diff.a","tot.diff.b","tot.diff.c"),
        c("pol.undiff.a","pol.undiff.b","pol.undiff.c"),
        c("pol.diff.a","pol.diff.b","pol.diff.c"),
        data.type="array", label.level=c("transcriptome","translatome"),
        label.condition=c("undifferentiated","differentiated"))
}

\keyword{TranslatomeDataset}
\keyword{DEGs}
\keyword{differential expression}
\keyword{expression matrix}
