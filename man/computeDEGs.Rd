\name{computeDEGs}
\alias{computeDEGs}

\title{computeDEGsHelpfile}
\description{
This function takes as input an object of the class \code{\linkS4class{TranslatomeDataset}} which contains a normalized data matrix coming from high throughput experiment. 
It takes as an input a character label specifying the method that we want to employ in order to detect DEGs(t-test, ANOTA, DESeq, edgeR, SAM, RP, limma) and returns an object of the class \code{\linkS4class{DEGs}}, in which each gene is assigned an expression class: up- or down-regulated at the first level, up- or down-regulated at the second level, up-regulated at both levels, down-regulated at both levels, up-regulated at the first level and down-regulated at the second level and vice versa.
}
\usage{computeDEGs(object, method="limma", significance.threshold= 0.05,
    FC.threshold= 0, log.transformed = FALSE, mult.cor=TRUE)}
\arguments{
 \item{object}{an object of class \code{\linkS4class{TranslatomeDataset}}}
 \item{method}{a character string that specifies the method that the user wants to employ in the differential expression analysis. It can have one the following values: \code{limma}, \code{t-test}, \code{RP}, \code{SAM}, \code{ANOTA}, \code{DESeq} and \code{none}.By default, this value is set to \code{limma},}
 \item{significance.threshold}{a numeric value specifying the threshold on the statistical significance below which the genes are considered as differentially expressed, the default is set to \code{0.05},}
 \item{FC.threshold}{a numeric value specifying the threshold on the absolute log2 fold change, above which the genes are considered as differentially expressed, the default is set to \code{0},}
  \item{log.transformed}{a boolean variable specifying whether the signals contained in expr.matrix have been previously log2 transformed. By default it is set to \code{FALSE},}
 \item{mult.cor}{a boolean variable specifying whether the significance threshold is applied on the multiple test corrected or on the original p-values obtained from the DEGs detection method. By default it is set to \code{TRUE},}
}
\details{Signals contained in expr.matrix should be previously normalized with standard methods (quantile, percentile shift, ... ) when data is coming from microarrays or in the appropriate cases when it is coming from sequencing experiments.}
\value{
 An object of class \code{\linkS4class{DEGs}}
}

\references{
Smyth GK. (2004) Linear models and empirical Bayes methods for assessing
differential expression in microarray experiments. Stat Appl Genet Mol
Biol., 3:Article3.

Tian L, Greenberg SA, Kong SW, Altschuler J, Kohane IS, Park PJ.(2005)
Discovering statistically significant pathways in expression profiing
studies. Proc Natl Acad Sci USA, 102(38):13544-9.

Breitling R, Armengaud P, Amtmann A, Herzyk P.(2004) Rank products: a
simple, yet powerful, new method to detect differentially regulated
genes in replicated microarray experiments. FEBS Lett., 573(1-3):83-92.

Tusher VG, Tibshirani R, Chu G.(2001) Significance analysis of
microarrays applied to the ionizing radiation response. Proc Natl Acad
Sci USA., 2001, 98(9):5116-21.

Larsson O, Sonenberg N, Nadon R.(2011) anota: Analysis of differential
translation in genome-wide studies. Bioinformatics, 27(10):1440-1.

Anders S, Huber W.(2010) Differential expression analysis for sequence
count data. Genome Biology, 11(10):R106.

Robinson MD, McCarthy DJ, Smyth GK.(2010) edgeR: a Bioconductor package
for differential expression analysis of digital gene expression
data. Bioinformatics, 26(1):139-40.
}

\author{
Toma Tebaldi, Erik Dassi, Galena Kostoska
}

\seealso{
\code{\linkS4class{TranslatomeDataset}}
\code{\linkS4class{DEGs}}
\code{\link{Scatterplot}}
\code{\link{Histogram}}
\code{\link{CVplot}}
\code{\link{MAplot}}
\code{\link{SDplot}}
}

\examples{
data(tRanslatomeSampleData)
computeDEGs(translatome.analysis, method= "limma", ,mult.cor=TRUE) 
}

\keyword{DEGs}
\keyword{differential expression}
\keyword{expression matrix}
