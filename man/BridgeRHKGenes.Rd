\name{BridgeRHKGenes}
\alias{BridgeRHKGenes}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Classic BRIC-seq data analysis pipeline
}
\description{
RNA half-life for each gene is calculated using the RPKM of house-keeping genes.
}
\usage{
BridgeRHKGenes(InputFiles,
                           InforColumn=4, 
                           InforHKGenes=2, 
                           HKGenes=c("GAPDH","PGK1","PPIA","ENO1","ATP5B","ALDOA"), 
                           group, 
                           hour, 
                           RPKMcutoff=0.1,
                           CutoffDataPointNumber = 4,
                           CutoffDataPoint1 = c(1,2),
                           CutoffDataPoint2 = c(8,12),
                           ThresholdHalfLife = c(8,12),
                           CutoffRelExp=0.001,
                           ModelMode="R2_selection")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{InputFiles}{
Vector
}
  \item{InforColumn}{
Integer
}
  \item{InforHKGenes}{
Integer
}
  \item{HKGenes}{
String(gene_name)
}
  \item{group}{
Vector(string)
}
  \item{hour}{
Vector(number)
}
  \item{RPKMcutoff}{
Float
}
  \item{CutoffDataPointNumber}{
Integer
}
  \item{CutoffDataPoint1}{
Vector(Integer or float)
}
  \item{CutoffDataPoint2}{
Vector(Integer or float)
}
  \item{ThresholdHalfLife}{
Vector(Integer or float)
}
  \item{CutoffRelExp}{
Float
}
  \item{ModelMode}{
Bool(True or False)
}
}
\details{
RNA half-life for each gene is calculated using the RPKM of house-keeping genes.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
text files and fig files
}
\references{
https://github.com/Naoto-Imamachi/BRIC-seq_data_analysis/tree/master/BridgeR
}
\author{
Naoto Imamachi
}
\note{
2015-11-05
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (InputFiles,
                           InforColumn=4, 
                           InforHKGenes=2, 
                           HKGenes=c("GAPDH","PGK1","PPIA","ENO1","ATP5B","ALDOA"), 
                           group, 
                           hour, 
                           RPKMcutoff=0.1,
                           CutoffDataPointNumber = 4,
                           CutoffDataPoint1 = c(1,2),
                           CutoffDataPoint2 = c(8,12),
                           ThresholdHalfLife = c(8,12),
                           CutoffRelExp=0.001,
                           ModelMode="R2_selection") 
{
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ BRIC-seq }
\keyword{ Alpha_version }% __ONLY ONE__ keyword per line
