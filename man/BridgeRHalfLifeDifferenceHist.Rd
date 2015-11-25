\name{BridgeRHalfLifeDifferenceHist}
\alias{BridgeRHalfLifeDifferenceHist}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Draw histgram to compare RNA half-life difference between two conditions.
}
\description{
Draw histgram to compare RNA half-life difference between two conditions.
}
\usage{
BridgeRHalfLifeDifferenceHist(filename = "BridgeR_5C_HalfLife_calculation_R2_selection.txt", 
                                          group, 
                                          hour, 
                                          ComparisonFile, 
                                          InforColumn = 4, 
                                          BinwidthFig = 0.01, 
                                          Calibration = F, 
                                          OutputFig = "BridgeR_7_HalfLife_Difference_LineGraph",
                                          ModelMode = "R2_selection")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{filename}{
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
  \item{BinwidthFig}{
Number(Integer or Float)
}
  \item{Calibration}{
Bool(True or False)
}
  \item{OutputFig}{
File path/name
}
  \item{ModelMode}{
Mode
}
}
\details{
Draw histgram to compare RNA half-life difference between two conditions.
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
#inputfile <- "BridgeR_5_HalfLife_calculation.txt"
#outputfig <- "BridgeR_5_HalfLife_Difference_LineGraph"
#group <- c("Control","knockdown1","knockdown2")
#hour <- c(0,1,2,4,8,12)
#compfile <- c("Control","Knockdown1")

#BridgeRHalfLifeDifferenceHist(filename=inputfile, group=group, hour=hour, ComparisonFile=compfile, OutputFig=outputfig)

## The function is currently defined as
function (filename = "BridgeR_5C_HalfLife_calculation_R2_selection.txt", 
                                          group, 
                                          hour, 
                                          ComparisonFile, 
                                          InforColumn = 4, 
                                          BinwidthFig = 0.01, 
                                          Calibration = F, 
                                          OutputFig = "BridgeR_7_HalfLife_Difference_LineGraph",
                                          ModelMode = "R2_selection")
{
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
