\name{BridgeRNormalizationForLuc2}
\alias{BridgeRNormalizationForLuc2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimate normalization factor using suctom normalization factor
}
\description{
Estimate normalization factor using suctom normalization factor such as the relative RPKM of house-keeping genes or Luciferase.
}
\usage{
BridgeRNormalizationForLuc2(filename, group, hour, InforColumn = 4, NormFactor = "BridgeR_3_Normalizaion_factor_dataset", OutputFile = "BridgeR_4_Normalized_expression_data.txt")
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
  \item{NormFactor}{
File path/name
}
  \item{OutputFile}{
File path/name
}
}
\details{
Estimate normalization factor using suctom normalization factor such as the relative RPKM of house-keeping genes or Luciferase.
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
#inputfile <- "BridgeR_1_Relative_expression_dataset.txt"
#nf <- "BridgeR_3_Normalizaion_factor_HouseKeepingGenes.txt"
#outputfile <- "BridgeR_4_Normalized_expression_dataset.txt"
#group <- c("Control","knockdown1","knockdown2")
#hour <- c(0,1,2,4,8,12)

#BridgeRNormalizationForLuc2(filename=inputfile, group=group, hour=hour, NormFactor=nf, OutputFile=outputfile)

## The function is currently defined as
function (filename, group, hour, InforColumn = 4, NormFactor = "BridgeR_3_Normalizaion_factor_dataset", 
    OutputFile = "BridgeR_4_Normalized_expression_data.txt") 
{
    group_number <- length(group)
    time_points <- length(hour)
    nf_st <- 2
    nf_ed <- time_points + 1
    nf_Luc2 <- NULL
    normalization_factor <- fread(NormFactor, header = T)[, nf_st:nf_ed, 
        with = F]
    for (a in 1:group_number) {
        if (is.null(nf_Luc2)) {
            nf_Luc2 <- as.vector(as.matrix(normalization_factor[a, 
                ]))
        }
        else {
            nf_Luc2 <- append(nf_Luc2, as.vector(as.matrix(normalization_factor[a, 
                ])))
        }
    }
    input_file <- fread(filename, header = T)
    output_file <- OutputFile
    cat("", file = output_file)
    hour_label <- NULL
    for (a in 1:group_number) {
        if (!is.null(hour_label)) {
            cat("\t", file = output_file, append = T)
        }
        hour_label <- NULL
        for (x in hour) {
            label <- x
            if (x < 10) {
                label <- paste("0", x, sep = "")
            }
            hour_label <- append(hour_label, paste("T", label, 
                "_", a, sep = ""))
        }
        infor_st <- 1 + (a - 1) * (time_points + InforColumn)
        infor_ed <- (InforColumn) * a + (a - 1) * time_points
        infor <- colnames(input_file)[infor_st:infor_ed]
        cat(infor, hour_label, sep = "\t", file = output_file, 
            append = T)
    }
    cat("\n", sep = "", file = output_file, append = T)
    gene_number <- length(input_file[[1]])
    for (x in 1:gene_number) {
        data <- as.vector(as.matrix(input_file[x, ]))
        for (a in 1:group_number) {
            if (a != 1) {
                cat("\t", sep = "", file = output_file, append = T)
            }
            infor_st <- 1 + (a - 1) * (time_points + InforColumn)
            infor_ed <- (InforColumn) * a + (a - 1) * time_points
            exp_st <- infor_ed + 1
            exp_ed <- infor_ed + time_points
            gene_infor <- data[infor_st:infor_ed]
            exp <- as.numeric(data[exp_st:exp_ed])
            nf_Luc2_st <- (time_points) * (a - 1) + 1
            nf_Luc2_ed <- (time_points) * a
            nf_Luc2_exp <- nf_Luc2[nf_Luc2_st:nf_Luc2_ed]
            normalized_exp <- NULL
            normalized_exp <- exp/nf_Luc2_exp
            cat(gene_infor, sep = "\t", file = output_file, append = T)
            cat("\t", sep = "\t", file = output_file, append = T)
            cat(normalized_exp, sep = "\t", file = output_file, 
                append = T)
        }
        cat("\n", sep = "\t", file = output_file, append = T)
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ BRIC-seq }
\keyword{ Alpha_version }% __ONLY ONE__ keyword per line
