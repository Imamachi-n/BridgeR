\name{BridgeRHalfLifeCalcR2Select}
\alias{BridgeRHalfLifeCalcR2Select}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculate RNA half-life from BRIC-seq data
}
\description{
Calculate RNA half-life from BRIC-seq data
}
\usage{
BridgeRHalfLifeCalcR2Select(InputFile = "BridgeR_4_Normalized_expression_dataset.txt",
                            group, 
                            hour, 
                            InforColumn = 4, 
                            CutoffDataPointNumber = 4,
                            CutoffDataPoint1 = c(1,2),
                            CutoffDataPoint2 = c(8,12),
                            ThresholdHalfLife = c(8,12),
                            OutputFile = "BridgeR_5C_HalfLife_calculation_R2_selection.txt")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{InputFile}{
File path/name
}
  \item{group}{
Vector(string)
}
  \item{hour}{
Vector(number)
}
  \item{InforColumn}{
Integer
}
  \item{CutoffDataPointNumber}{
Integer
}
  \item{CutoffDataPoint1}{
Vector(Integer)
}
  \item{CutoffDataPoint2}{
Vector(Integer)
}
  \item{ThresholdHalfLife}{
Vector(Integer)
}
  \item{OutputFile}{
File path/name
}
}
\details{
Calculate RNA half-life from BRIC-seq data
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
text file
}
\references{
https://github.com/Naoto-Imamachi/BRIC-seq_data_analysis/tree/master/BridgeR
}
\author{
Naoto Imamachi
}
\note{
2015-11-24
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
function (InputFile = "BridgeR_4_Normalized_expression_dataset.txt",
                            group, 
                            hour, 
                            InforColumn = 4, 
                            CutoffDataPointNumber = 4,
                            CutoffDataPoint1 = c(1,2),
                            CutoffDataPoint2 = c(8,12),
                            ThresholdHalfLife = c(8,12),
                            OutputFile = "BridgeR_5C_HalfLife_calculation_R2_selection.txt")
{
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ BRIC-seq }
\keyword{ Alpha_version }% __ONLY ONE__ keyword per line
