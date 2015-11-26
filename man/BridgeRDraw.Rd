\name{BridgeRDraw}
\alias{BridgeRDraw}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Draw fitting decay curve for BRIC-seq data
}
\description{
Draw fitting decay curve for BRIC-seq data
}
\usage{
BridgeRDraw(InputFile="BridgeR_5C_HalfLife_calculation_R2_selection.txt",
                        group,
                        hour,
                        ComparisonFile,
                        CutoffRelExp = 0,
                        CutoffDataPoint = 4,
                        InforColumn = 4,
                        GeneInfor = 2,
                        OutputDir = "BridgeRDraw_fig",
                        ModelMode = "R2_selection",
                        DrawMode = "Simple",
                        Color = c("black","red"))
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
  \item{CutoffRelExp}{
Number(Integer or Float)
}
  \item{CutoffDataPoint}{
integer
}
  \item{InforColumn}{
Integer
}
  \item{GeneInfor}{
Integer
}
  \item{OutputDir}{
File path/directory name
}
  \item{ModelMode}{
String
}
  \item{DrawMode}{
String
}
  \item{Color}{
Vector(string)
}
}
\details{
Draw fitting decay curve for BRIC-seq data
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
fig files
}
\references{
https://github.com/Naoto-Imamachi/BRIC-seq_data_analysis/tree/master/BridgeR
}
\author{
Naoto Imamachi
}
\note{
2015-11-26
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
function(InputFile="BridgeR_5C_HalfLife_calculation_R2_selection.txt",
                        group,
                        hour,
                        ComparisonFile,
                        CutoffRelExp = 0,
                        CutoffDataPoint = 4,
                        InforColumn = 4,
                        GeneInfor = 2,
                        OutputDir = "BridgeRDraw_fig",
                        ModelMode = "R2_selection",
                        DrawMode = "Simple",
                        Color = c("black","red"))
{
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ BRIC-seq }
\keyword{ Alpha_version }% __ONLY ONE__ keyword per line
