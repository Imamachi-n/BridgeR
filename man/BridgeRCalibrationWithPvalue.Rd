\name{BridgeRCalibrationWithPvalue}
\alias{BridgeRCalibrationWithPvalue}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculate RNA half-life from BRIC-seq data
}
\description{
Calculate RNA half-life from BRIC-seq data
}
\usage{
BridgeRCalibrationWithPvalue(InputFile="BridgeR_6_HalfLife_Pvalue_estimation.txt",
                        group,
                        hour,
                        ComparisonFile,
                        InforColumn=4,
                        CutoffDataPoint=4,
                        OutputFile="BridgeR_7_HalfLife_Calibration.txt")
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
  \item{ComparisonFile}{
Vector(string)
}
  \item{InforColumn}{
Integer
}
  \item{CutoffDataPoint}{
Integer
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
2015-11-20
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
function (InputFile="BridgeR_6_HalfLife_Pvalue_estimation.txt",
                        group,
                        hour,
                        ComparisonFile,
                        InforColumn=4,
                        CutoffDataPoint=4,
                        OutputFile="BridgeR_7_HalfLife_Calibration.txt")
{
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ BRIC-seq }
\keyword{ Alpha_version }% __ONLY ONE__ keyword per line
