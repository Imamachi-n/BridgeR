\name{BridgeRHalfLifeCalculation}
\alias{BridgeRHalfLifeCalculation}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculate RNA half-life from BRIC-seq data
}
\description{
Calculate RNA half-life from BRIC-seq data
}
\usage{
BridgeRHalfLifeCalculation(filename = "BridgeR_4_Normalized_expression_data.txt", group, hour, InforColumn = 4, CutoffRelExp = 0.1, CutoffDataPoint = 3, OutputFile = "BridgeR_5_half-life_calculation.txt")
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
    OutputFile = "BridgeR_5_half-life_calculation.txt") 
{
    time_points <- length(hour)
    group_number <- length(group)
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
            hour_label <- append(hour_label, paste("T", x, "_", 
                a, sep = ""))
        }
        infor_st <- 1 + (a - 1) * (time_points + InforColumn)
        infor_ed <- (InforColumn) * a + (a - 1) * time_points
        infor <- colnames(input_file)[infor_st:infor_ed]
        cat(infor, hour_label, sep = "\t", file = output_file, 
            append = T)
        cat("\t", sep = "", file = output_file, append = T)
        cat("Model", "Decay_rate_coef", "coef_error", "coef_p-value", 
            "R2", "Adjusted_R2", "Residual_standard_error", "half_life", 
            sep = "\t", file = output_file, append = T)
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
            cat(gene_infor, sep = "\t", file = output_file, append = T)
            cat("\t", file = output_file, append = T)
            exp <- as.numeric(data[exp_st:exp_ed])
            cat(exp, sep = "\t", file = output_file, append = T)
            cat("\t", file = output_file, append = T)
            time_point_exp <- data.frame(hour, exp)
            time_point_exp <- time_point_exp[time_point_exp$exp >= 
                CutoffRelExp, ]
            data_point <- length(time_point_exp$exp)
            if (!is.null(time_point_exp)) {
                if (data_point >= CutoffDataPoint) {
                  model <- lm(log(time_point_exp$exp) ~ time_point_exp$hour - 
                    1)
                  model_summary <- summary(model)
                  coef <- -model_summary$coefficients[1]
                  coef_error <- model_summary$coefficients[2]
                  coef_p <- model_summary$coefficients[4]
                  r_squared <- model_summary$r.squared
                  adj_r_squared <- model_summary$adj.r.squared
                  residual_standard_err <- model_summary$sigma
                  half_life <- log(2)/coef
                  if (coef < 0 || half_life >= 24) {
                    half_life <- 24
                  }
                  cat("Exponential_Decay_Model", coef, coef_error, 
                    coef_p, r_squared, adj_r_squared, residual_standard_err, 
                    half_life, sep = "\t", file = output_file, 
                    append = T)
                }
                else {
                  cat("few_data", "NA", "NA", "NA", "NA", "NA", 
                    "NA", "NA", sep = "\t", file = output_file, 
                    append = T)
                }
            }
            else {
                cat("low_expresion", "NA", "NA", "NA", "NA", 
                  "NA", "NA", "NA", sep = "\t", file = output_file, 
                  append = T)
            }
        }
        cat("\n", file = output_file, append = T)
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ BRIC-seq }
\keyword{ Alpha_version }% __ONLY ONE__ keyword per line
