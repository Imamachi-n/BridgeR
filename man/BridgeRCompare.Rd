\name{BridgeRCompare}
\alias{BridgeRCompare}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Comparison of two BRIC-seq datasets}
\description{
Comparison of two BRIC-seq datasets
}
\usage{
BridgeRCompare(InputFile="BridgeR_5C_HalfLife_calculation_R2_selection.txt",
                           InforColumn = 4,
                           group, 
                           hour, 
                           ComparisonFile,
                           CutoffDataPointNumber = 4,
                           ModelMode="R2_selection")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{InputFile}{
Your RPKM value data
}
  \item{InforColumn}{
information column length
}
  \item{group}{
Group
}
  \item{hour}{
Time course
}
  \item{ComparisonFile}{
Group which you want to compare 
}
  \item{CutoffDataPointNumber}{
Integer 
}
  \item{ModelMode}{
Mode 
}
}
\details{
Comparison of two BRIC-seq datasets
}
\value{
Fig files
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
function (InputFile="BridgeR_5C_HalfLife_calculation_R2_selection.txt",
                           InforColumn = 4,
                           group, 
                           hour, 
                           ComparisonFile,
                           CutoffDataPointNumber = 4,
                           ModelMode="R2_selection") 
{
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ BRIC-seq }
\keyword{ Alpha_version }% __ONLY ONE__ keyword per line
