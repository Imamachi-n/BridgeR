\name{BridgeRCalibration}
\alias{BridgeRCalibration}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Calibration of two BRIC-seq datasets}
\description{
Calibration of two BRIC-seq datasets
}
\usage{
BridgeRCalibration(InputFile, InforColumn = 4, group, hour, ComparisonFile)
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
}
\details{
Calibration of two BRIC-seq datasets
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
## The function is currently defined as
function(InputFile="BridgeR_5_HalfLife_calculation.txt", InforColumn=4, group, hour, ComparisonFile){
    library(data.table)
    library(ggplot2)
    BridgeRHalfLifeCalibration(InputFile=InputFile, group=group, hour=hour, ComparisonFile=ComparisonFile, InforColumn=4, OutputFile="BridgeR_6_HalfLife_Calibration")
    BridgeRHalfLifeComparison(filename=InputFile, InforColumn=4, group=group, hour=hour, ComparisonFile=ComparisonFile, LogScale=F, Calibration=T, OutputFig="BridgeR_6_Adjusted_HalfLife_Comparison_ScatteredPlot")
    BridgeRHalfLifeDifferenceHist(filename=InputFile, InforColumn=4, group=group, hour=hour, ComparisonFile=ComparisonFile, BinwidthFig=0.01, Calibration=T, OutputFig="BridgeR_6_Adjusted_HalfLife_Difference_LineGraph")
    BridgeRHalfLifeDifferenceBox(filename=InputFile, InforColumn=4, group=group, hour=hour, ComparisonFile=ComparisonFile, Calibration=T, OutputFig= "BridgeR_6_Adjusted_HalfLife_Comparison_BoxPlot")
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ BRIC-seq }
\keyword{ Alpha_version }% __ONLY ONE__ keyword per line
