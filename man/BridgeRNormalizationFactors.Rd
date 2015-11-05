\name{BridgeRNormalizationFactors}
\alias{BridgeRNormalizationFactors}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimate normalization factor
}
\description{
Estimate normalization factor
}
\usage{
BridgeRNormalizationFactors(InputFile, group, hour, InforColumn = 4, YMin = -2, YMax = 2, MakeFig = T, figname = "BridgeR_3_Normalizaion_factor", nfname = "BridgeR_3_Normalizaion_factor")
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
  \item{InforColumn}{
Integer
}
  \item{YMin}{
Number(Integer or Float)
}
  \item{YMax}{
Number(Integer or Float)
}
  \item{MakeFig}{
Bool(True or False)
}
  \item{figname}{
File path/name
}
  \item{nfname}{
File path/name
}
}
\details{
Estimate normalization factor
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
#figname <- "BridgeR_3_Normalizaion_factor"
#outputfile <- "BridgeR_3_Normalizaion_factor"
#group <- c("Control","knockdown1","knockdown2")
#hour <- c(0,1,2,4,8,12)

#BridgeRNormalizationFactors(InputFile=inputfile,group=group, hour=hour, figname=figname, nfname=outputfile)

## The function is currently defined as
function (InputFile, group, hour, InforColumn = 4, YMin = -2, 
    YMax = 2, MakeFig = T, figname = "BridgeR_3_Normalizaion_factor", 
    nfname = "BridgeR_3_Normalizaion_factor") 
{
    group_number <- length(group)
    time_points <- length(hour)
    input_file <- fread(InputFile, header = T)
    for (a in 1:group_number) {
        infor_st <- 1 + (a - 1) * (time_points + InforColumn)
        infor_ed <- (InforColumn) * a + (a - 1) * time_points
        exp_st <- infor_ed + 1
        exp_ed <- infor_ed + time_points
        hour_label <- NULL
        for (x in hour) {
            label <- x
            if (x < 10) {
                label <- paste("0", x, sep = "")
            }
            hour_label <- append(hour_label, paste("T", label, 
                "_", a, sep = ""))
        }
        quantile_99_data <- NULL
        quantile_95_data <- NULL
        for (x in exp_st:exp_ed) {
            each_time_exp <- input_file[[x]]
            each_time_quantile_99 <- quantile(each_time_exp, 
                prob = 0.99, na.rm = T)
            each_time_quantile_95 <- quantile(each_time_exp, 
                prob = 0.975, na.rm = T)
            quantile_99_data <- append(quantile_99_data, each_time_quantile_99)
            quantile_95_data <- append(quantile_95_data, each_time_quantile_95)
        }
        output_filename <- paste(nfname, "_", group[a], ".txt", 
            sep = "")
        cat("Percentile", hour_label, sep = "\t", file = output_filename)
        cat("\n", file = output_filename, append = T)
        cat("99\%_percentale", quantile_99_data, sep = "\t", file = output_filename, 
            append = T)
        cat("\n", file = output_filename, append = T)
        cat("97.5\%_percentale", quantile_95_data, sep = "\t", 
            file = output_filename, append = T)
        cat("\n", file = output_filename, append = T)
        if (MakeFig == TRUE) {
            quantile_99_data <- log10(as.vector(quantile_99_data))
            quantile_95_data <- log10(as.vector(quantile_95_data))
            plot_data_quantile_99 <- data.frame(hour, quantile_99_data)
            plot_data_quantile_95 <- data.frame(hour, quantile_95_data)
            figfile <- paste(figname, "_", group[a], ".png", 
                sep = "")
            gene_number <- length(input_file[[1]])
            exp_data <- as.matrix(input_file[, exp_st:exp_ed, 
                with = F])
            exp_data <- t(exp_data)
            exp_data <- factor(exp_data)
            exp_data <- as.numeric(as.character(exp_data))
            exp_data <- log10(exp_data)
            time_data <- as.numeric(rep(hour, gene_number))
            class_data <- NULL
            for (x in 1:gene_number) {
                class_data <- append(class_data, rep(x, time_points))
            }
            plot_data <- data.frame(class_data, time_data, exp_data)
            png(filename = figfile, width = 1200, height = 1200)
            p.scatter <- ggplot()
            p.scatter <- p.scatter + layer(data = plot_data, 
                mapping = aes(x = time_data, y = exp_data, group = class_data), 
                geom = "line", colour = "black", size = 0.02, 
                alpha = 0.05)
            p.scatter <- p.scatter + layer(data = plot_data_quantile_99, 
                mapping = aes(x = hour, y = quantile_99_data), 
                geom = "line", colour = "red", size = 0.5, alpha = 1)
            p.scatter <- p.scatter + layer(data = plot_data_quantile_95, 
                mapping = aes(x = hour, y = quantile_95_data), 
                geom = "line", colour = "blue", size = 0.5, alpha = 1)
            p.scatter <- p.scatter + xlim(0, max(plot_data$time_data)) + 
                ylim(YMin, YMax)
            p.scatter <- p.scatter + ggtitle("All genes distribution")
            p.scatter <- p.scatter + xlab("Time course")
            p.scatter <- p.scatter + ylab("Relative RPKM (Time0 = 1)")
            plot(p.scatter)
            dev.off()
            plot.new()
        }
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ BRIC-seq }
\keyword{ Alpha_version }% __ONLY ONE__ keyword per line
