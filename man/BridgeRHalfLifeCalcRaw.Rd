\name{BridgeRHalfLifeCalcRaw}
\alias{BridgeRHalfLifeCalcRaw}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculate RNA half-life from BRIC-seq data
}
\description{
Calculate RNA half-life from BRIC-seq data
}
\usage{
BridgeRHalfLifeCalcRaw(filename = "BridgeR_4_Normalized_expression_data.txt", group, hour, InforColumn = 4, CutoffRelExp = 0.1, CutoffDataPoint = 3, OutputFile = "BridgeR_5_HalfLife_calculation.txt")
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
  \item{InforColumn}{
Integer
}
  \item{CutoffRelExp}{
Number(Integer or Float)
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
#inputfile <- "BridgeR_4_Normalized_expression_dataset.txt"
#outputfile <- "BridgeR_5_HalfLife_calculation.txt"
#group <- c("Control","knockdown1","knockdown2")
#hour <- c(0,1,2,4,8,12)

#BridgeRHalfLifeCalculation(filename=inputfile, group=group, hour=hour, OutputFile=outputfile)

## The function is currently defined as
function (filename = "BridgeR_4_Normalized_expression_data.txt", 
    group, hour, InforColumn = 4, CutoffRelExp = 0.1, CutoffDataPoint = 3, 
    OutputFile = "BridgeR_5_HalfLife_calculation.txt") 
{
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ BRIC-seq }
\keyword{ Alpha_version }% __ONLY ONE__ keyword per line
