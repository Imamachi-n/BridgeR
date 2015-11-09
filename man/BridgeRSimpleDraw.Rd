\name{BridgeRSimpleDraw}
\alias{BridgeRSimpleDraw}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Draw fitting decay curve for BRIC-seq data
}
\description{
Draw fitting decay curve for BRIC-seq data
}
\usage{
BridgeRSimpleDraw(filename="BridgeR_4_Normalized_expression_dataset.txt", group, hour, ComparisonFile, CutoffRelExp = 0.1, CutoffDataPoint = 3, InforColumn = 4, OutputDir = "BridgeR_simple_fig")
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
  \item{CutoffRelExp}{
Number(Integer or Float)
}
  \item{CutoffDataPoint}{
integer
}
  \item{InforColumn}{
Integer
}
  \item{OutputDir}{
File path/directory namr
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
2015-11-09
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
function(filename="BridgeR_4_Normalized_expression_dataset.txt", group, hour, ComparisonFile,
         CutoffRelExp = 0.1, CutoffDataPoint = 3, InforColumn = 4, OutputDir = "BridgeR_simple_fig")
{
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ BRIC-seq }
\keyword{ Alpha_version }% __ONLY ONE__ keyword per line
