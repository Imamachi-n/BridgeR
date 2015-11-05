\name{BridgeRDataSetFromCuffnorm}
\alias{BridgeRDataSetFromCuffnorm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Prepare your BRIC-seq dataset analyzed
}
\description{
Prepare your BRIC-seq dataset analyzed
}
\usage{
BridgeRDataSetFromCuffnorm(CuffnormFiles, group, hour, cutoff = 0.1, InforColumn = 4, OutputFile = "BridgeR_1_Relative_expression_data.txt")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{CuffnormFiles}{
File path/name
}
  \item{group}{
Vector(string)
}
  \item{hour}{
Vector(number)
}
  \item{cutoff}{
Number(integer, float)
}
  \item{InforColumn}{
Integer
}
  \item{OutputFile}{
File path/name
}
}
\details{
Prepare your BRIC-seq dataset analyzed
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
#inputfile <- c("Control_genes_RefSeq_result_mRNA.fpkm_table",
#               "Knockdown_genes_RefSeq_result_mRNA.fpkm_table")
#outputfile <- "BridgeR_1_Relative_expression_dataset.txt"
#group <- c("Control","knockdown1","knockdown2")
#hour <- c(0,1,2,4,8,12)

#BridgeRDataSetFromCuffnorm(CuffnormFiles=inputfile, group=group, hour=hour, OutputFile=outputfile)

## The function is currently defined as
function (CuffnormFiles, group, hour, cutoff = 0.1, InforColumn = 4, 
    OutputFile = "BridgeR_1_Relative_expression_data.txt") 
{
    time_points <- length(hour)
    input_file_numbers <- length(CuffnormFiles)
    input_file <- NULL
    for (filename in CuffnormFiles) {
        if (is.null(input_file)) {
            input_file <- suppressWarnings(fread(filename, header = T))
        }
        else {
            input_file <- cbind(input_file, suppressWarnings(fread(filename, 
                header = T)))
        }
    }
    output_file <- OutputFile
    cat("", file = output_file)
    hour_label <- NULL
    for (a in 1:length(group)) {
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
    sample_size <- length(group)
    for (x in 1:gene_number) {
        data <- as.vector(as.matrix(input_file[x, ]))
        for (a in 1:sample_size) {
            if (a != 1) {
                cat("\t", sep = "", file = output_file, append = T)
            }
            infor_st <- 1 + (a - 1) * (time_points + InforColumn)
            infor_ed <- (InforColumn) * a + (a - 1) * time_points
            gene_infor <- data[infor_st:infor_ed]
            cat(gene_infor, sep = "\t", file = output_file, append = T)
            cat("\t", sep = "", file = output_file, append = T)
            exp_st <- infor_ed + 1
            exp_ed <- infor_ed + time_points
            exp <- data[exp_st:exp_ed]
            exp <- as.numeric(exp)
            start_time <- exp[1]
            if (start_time <= cutoff) {
                cat(rep(0, time_points), sep = "\t", file = output_file, 
                  append = T)
                next
            }
            rel_exp <- NULL
            for (y in 1:time_points) {
                rel_exp <- append(rel_exp, exp[y]/start_time)
            }
            cat(rel_exp, sep = "\t", file = output_file, append = T)
        }
        cat("\n", sep = "", file = output_file, append = T)
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ BRIC-seq }
\keyword{ Alpha_version }% __ONLY ONE__ keyword per line
