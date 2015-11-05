\name{BridgeRDatasetChecker}
\alias{BridgeRDatasetChecker}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Check BRIC-seq dataset
}
\description{
Check BRIC-seq dataset
}
\usage{
BridgeRDatasetChecker(InputFile, group, hour, InforColumn = 4, OutputFile = "BridgeR_2_Relative_RPKM_distribution")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{InputFile}{
Vector
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
  \item{OutputFile}{
fig files
}
}
\details{
Check BRIC-seq dataset
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
#outputfile <- "BridgeR_2_Relative_RPKM_distribution"
#group <- c("Control","knockdown1","knockdown2")
#hour <- c(0,1,2,4,8,12)

#BridgeRDatasetChecker(InputFile=inputfile, OutputFile=outputfile, group=group, hour=hour)

## The function is currently defined as
function (InputFile, group, hour, InforColumn = 4, OutputFile = "BridgeR_2_Relative_RPKM_distribution") 
{
    time_points <- length(hour)
    input_file <- fread(InputFile, header = T)
    sample_size <- length(group)
    test_data <- NULL
    if (sample_size == 1) {
        test_data <- input_file[T00_1 == 1, ]
    }
    else if (sample_size == 2) {
        test_data <- input_file[T00_1 == 1 & T00_2 == 1, ]
    }
    else if (sample_size == 3) {
        test_data <- input_file[T00_1 == 1 & T00_2 == 1 & T00_3 == 
            1, ]
    }
    else if (sample_size == 4) {
        test_data <- input_file[T00_1 == 1 & T00_2 == 1 & T00_3 == 
            1 & T00_4 == 1, ]
    }
    merge_fig_data <- NULL
    merge_fig_percentile_data <- NULL
    for (a in 1:sample_size) {
        infor_st <- 1 + (a - 1) * (time_points + InforColumn)
        infor_ed <- (InforColumn) * a + (a - 1) * time_points
        exp_st <- infor_ed + 1
        exp_ed <- infor_ed + time_points
        hour_label <- NULL
        for (x in hour) {
            if (x == 0) {
                next
            }
            label <- x
            if (x < 10) {
                label <- paste("0", x, sep = "")
            }
            hour_label <- append(hour_label, paste(label, "hr_", 
                group[a], sep = ""))
        }
        exp_st <- exp_st + 1
        exp_data <- test_data[, exp_st:exp_ed, with = F]
        exp_percentile_data <- NULL
        time_points_for_fig <- time_points - 1
        for (x in 1:time_points_for_fig) {
            q_data <- test_q(log10(exp_data[[x]]), hour_label[x])
            if (x == 1) {
                exp_percentile_data <- q_data
            }
            else {
                exp_percentile_data <- rbind(exp_percentile_data, 
                  q_data)
            }
        }
        exp_data <- t(exp_data)
        exp_data <- factor(exp_data)
        exp_data <- as.numeric(as.character(exp_data))
        exp_data <- log10(exp_data)
        gene_number <- length(test_data[[1]])
        label_data <- rep(hour_label, gene_number)
        fig_data <- data.frame(exp = exp_data, label = factor(label_data))
        if (a == 1) {
            merge_fig_data <- fig_data
            merge_fig_percentile_data <- exp_percentile_data
        }
        else {
            merge_fig_data <- rbind(merge_fig_data, fig_data)
            merge_fig_percentile_data <- rbind(merge_fig_percentile_data, 
                exp_percentile_data)
        }
        fig_name <- paste(OutputFile, "_Boxplot_", group[a], 
            ".png", sep = "")
        fig_width <- 150 * (time_points - 1)
        png(filename = fig_name, width = fig_width, height = 1200)
        p <- ggplot()
        p <- p + layer(data = fig_data, mapping = aes(x = label, 
            y = exp), geom = "boxplot")
        p <- p + ylim(-2, 2)
        plot(p)
        dev.off()
        plot.new()
        fig_name <- paste(OutputFile, "_Density_", group[a], 
            ".png", sep = "")
        png(filename = fig_name, width = 1300, height = 1000)
        p <- ggplot()
        p <- p + layer(data = fig_data, mapping = aes(x = exp, 
            colour = label), geom = "line", stat = "density", 
            size = 1.2)
        p <- p + xlim(-2, 2) + ylim(0, 7)
        plot(p)
        dev.off()
        plot.new()
        fig_name <- paste(OutputFile, "_Point_", group[a], ".png", 
            sep = "")
        fig_width <- 120 * (time_points - 1)
        png(filename = fig_name, width = fig_width, height = 1200)
        p <- ggplot()
        p <- p + layer(data = exp_percentile_data, mapping = aes(x = name, 
            y = q, colour = factor(factor)), geom = "point", 
            size = 5, shape = 19)
        p <- p + xlab("") + ylab("Relative RPKM (Time0 = 1)")
        p <- p + ylim(-1.5, 1.5)
        plot(p)
        dev.off()
        plot.new()
    }
    fig_name <- NULL
    for (a in 1:sample_size) {
        if (a == 1) {
            fig_name <- paste(OutputFile, "_Boxplot_", group[a], 
                sep = "")
        }
        else {
            fig_name <- paste(fig_name, "_", group[a], sep = "")
        }
    }
    fig_name <- paste(fig_name, ".png", sep = "")
    fig_width <- 150 * (time_points - 1) * sample_size
    png(filename = fig_name, width = fig_width, height = 1200)
    merge_fig_data$label <- factor(merge_fig_data$label, levels = sort(unique(as.character(merge_fig_data$label))))
    p <- ggplot()
    p <- p + layer(data = merge_fig_data, mapping = aes(x = label, 
        y = exp), geom = "boxplot")
    p <- p + ylim(-2, 2)
    plot(p)
    dev.off()
    plot.new()
    fig_name <- NULL
    for (a in 1:sample_size) {
        if (a == 1) {
            fig_name <- paste(OutputFile, "_Point_", group[a], 
                sep = "")
        }
        else {
            fig_name <- paste(fig_name, "_", group[a], sep = "")
        }
    }
    fig_name <- paste(fig_name, ".png", sep = "")
    fig_width <- 110 * (time_points - 1) * sample_size
    png(filename = fig_name, width = fig_width, height = 1200)
    merge_fig_percentile_data$name <- factor(merge_fig_percentile_data$name, 
        levels = sort(unique(as.character(merge_fig_percentile_data$name))))
    p <- ggplot()
    p <- p + layer(data = merge_fig_percentile_data, mapping = aes(x = name, 
        y = q, colour = factor(factor)), geom = "point", size = 5, 
        shape = 19)
    p <- p + xlab("") + ylab("Relative RPKM (Time0 = 1)")
    p <- p + ylim(-1.5, 1.5)
    plot(p)
    dev.off()
    plot.new()
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ BRIC-seq }
\keyword{ Alpha_version }% __ONLY ONE__ keyword per line
