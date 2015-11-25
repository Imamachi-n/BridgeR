\name{BridgeRHalfLifeCalc}
\alias{BridgeRHalfLifeCalc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculation of RNA half-life for each gene
}
\description{
Calculation of RNA half-life for each gene.
}
\usage{
BridgeRHalfLifeCalc(InputFiles = "BridgeR_4_Normalized_expression_dataset.txt", 
                                InforColumn = 4, 
                                group, 
                                hour,
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
  \item{group}{
Vector(string)
}
  \item{hour}{
Vector(number)
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
Calculation of RNA half-life for each gene.
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
2015-11-25
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
function (InputFiles = "BridgeR_4_Normalized_expression_dataset.txt", 
                                InforColumn = 4, 
                                group, 
                                hour,
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
