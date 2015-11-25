\name{BridgeRHalfLifeCalibration}
\alias{BridgeRHalfLifeCalibration}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calibrate RNA half-life.
}
\description{
Calibrate RNA half-life.
}
\usage{
BridgeRHalfLifeCalibration(InputFile = "BridgeR_5B_HalfLife_calculation_3model.txt", 
                                       group, 
                                       hour, 
                                       ComparisonFile, 
                                       InforColumn = 4, 
                                       OutputFile = "BridgeR_6_Calibrated_HalfLife_calculation_3model.txt")
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
  \item{OutputFile}{
File path/name
}
}
\details{
Calibrate RNA half-life.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
fig file
}
\references{
https://github.com/Naoto-Imamachi/BRIC-seq_data_analysis/tree/master/BridgeR
}
\author{
Naoto Imamachi
}
\note{
2015-11-06
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

#BridgeRHalfLifeComparison(filename=inputfile, group=group, hour=hour, ComparisonFile=compfile, OutputFig=outputfig)

## The function is currently defined as
function(InputFile = "BridgeR_5B_HalfLife_calculation_3model.txt", 
                                       group, 
                                       hour, 
                                       ComparisonFile, 
                                       InforColumn = 4, 
                                       OutputFile = "BridgeR_6_Calibrated_HalfLife_calculation_3model.txt")
{
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ BRIC-seq }
\keyword{ Alpha_version }% __ONLY ONE__ keyword per line
