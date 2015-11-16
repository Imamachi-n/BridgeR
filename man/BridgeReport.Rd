\name{BridgeReport}
\alias{BridgeReport}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Draw fitting decay curve for BRIC-seq data
}
\description{
Draw fitting decay curve for BRIC-seq data
}
\usage{
BridgeReport(filename="BridgeR_4_Normalized_expression_dataset.txt",
             group,
             hour,
             ComparisonFile,
             SearchRow="symbol",
             InforColumn=4)
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
  \item{SearchRow}{
String(header)
}
  \item{InforColumn}{
Integer
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
text and fig files
}
\references{
https://github.com/Naoto-Imamachi/BRIC-seq_data_analysis/tree/master/BridgeR
}
\author{
Naoto Imamachi
}
\note{
2015-11-16
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
function (filename="BridgeR_4_Normalized_expression_dataset.txt",
          group,
          hour,
          ComparisonFile,
          SearchRow="symbol",
          InforColumn=4)
{
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ BRIC-seq }
\keyword{ Alpha_version }% __ONLY ONE__ keyword per line
