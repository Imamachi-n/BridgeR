\name{BridgeRCompare}
\alias{BridgeRCompare}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Comparison of two BRIC-seq datasets}
\description{
Comparison of two BRIC-seq datasets
}
\usage{
BridgeRCompare(InputFile, InforColumn = 4, group, hour, ComparisonFile)
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
#file <- "BridgeR_5_HalfLife_calculation.txt"
#group <- c("Control","Knockdown1","Knockdown2")
#hour <- c(0,1,2,4,8,12)
#compfile <- c("Control","Knockdown1")

#BridgeRCompare(InputFile=file, group=group, hour=hour, ComparisonFile=compfile)

## The function is currently defined as
function (InputFile, InforColumn = 4, group, hour, ComparisonFile) 
{
    library(data.table)
    library(ggplot2)
    BridgeRHalfLifeComparison(filename = InputFile, InforColumn = 4, 
        group = group, hour = hour, ComparisonFile = ComparisonFile, 
        LogScale = F, OutputFig = "BridgeR_5_HalfLife_Comparison_ScatteredPlot")
    BridgeRHalfLifeDistribution(filename = InputFile, InforColumn = 4, 
        group = group, hour = hour, ComparisonFile = ComparisonFile, 
        OutputFig = "BridgeR_5_HalfLife_Distribution_LineGraph")
    BridgeRHalfLifeDifferenceHist(filename = InputFile, InforColumn = 4, 
        group = group, hour = hour, ComparisonFile = ComparisonFile, 
        BinwidthFig = 0.01, OutputFig = "BridgeR_5_HalfLife_Difference_LineGraph")
    BridgeRHalfLifeDifferenceBox(filename = InputFile, InforColumn = 4, 
        group = group, hour = hour, ComparisonFile = ComparisonFile, 
        OutputFig = "BridgeR_5_HalfLife_Comparison_BoxPlot")
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ BRIC-seq }
\keyword{ Alpha_version }% __ONLY ONE__ keyword per line
