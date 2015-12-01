\name{BridgeRHalfSD}
\alias{BridgeRHalfSD}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculate RNA half-life SD
}
\description{
Calculate RNA half-life SD
}
\usage{
BridgeRHalfSD(HalflifeFiles = c("BridgeR_6_HalfLife_Pvalue_estimation_PUM1_study.txt",
                                "BridgeR_6_HalfLife_Pvalue_estimation_PUM2_study.txt",
                                "BridgeR_6_HalfLife_Pvalue_estimation_PUM1_2_study.txt"),
              RPKMFiles = c("siStealth_PUM1_study_compatible_genes_RefSeq_result_mRNA.fpkm_table",
                            "siCTRL_PUM2_study_compatible_genes_RefSeq_result_mRNA.fpkm_table",
                            "siCTRL_PUM1_2_study_RefSeq_compatible_genes_RefSeq_result_mRNA.fpkm_table"),
              hour = c(0,1,2,4,8,12),
              InforColumn=4,
              CompFiles = c(1,1,1),
              CompNames = c("siStealth_PUM1_study","siCTRL_PUM2_study","siCTRL_PUM1_2_study"),
              OutputFile = "HalfLife_RPKM_mean_SD.txt")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{HalflifeFiles}{
File path/name
}
  \item{RPKMFiles}{
File path/name
}
  \item{hour}{
Vector(number)
}
  \item{InforColumn}{
Integer
}
  \item{CompFiles}{
Vector(integer)
}

  \item{CompNames}{
File path/name
}
  \item{OutputFile}{
File path/name
}
}
\details{
Calculate RNA half-life SD
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
2015-12-01
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
function (HalflifeFiles = c("BridgeR_6_HalfLife_Pvalue_estimation_PUM1_study.txt",
                                "BridgeR_6_HalfLife_Pvalue_estimation_PUM2_study.txt",
                                "BridgeR_6_HalfLife_Pvalue_estimation_PUM1_2_study.txt"),
              RPKMFiles = c("siStealth_PUM1_study_compatible_genes_RefSeq_result_mRNA.fpkm_table",
                            "siCTRL_PUM2_study_compatible_genes_RefSeq_result_mRNA.fpkm_table",
                            "siCTRL_PUM1_2_study_RefSeq_compatible_genes_RefSeq_result_mRNA.fpkm_table"),
              hour = c(0,1,2,4,8,12),
              InforColumn=4,
              CompFiles = c(1,1,1),
              CompNames = c("siStealth_PUM1_study","siCTRL_PUM2_study","siCTRL_PUM1_2_study"),
              OutputFile = "HalfLife_RPKM_mean_SD.txt")
{
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ BRIC-seq }
\keyword{ Alpha_version }% __ONLY ONE__ keyword per line
