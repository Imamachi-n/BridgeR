\name{BridgeRDrawFittingCurve}
\alias{BridgeRDrawFittingCurve}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Draw fitting decay curve for BRIC-seq data
}
\description{
Draw fitting decay curve for BRIC-seq data
}
\usage{
BridgeRDrawFittingCurve(filename, group, hour, ComparisonFile, CutoffRelExp = 0.1, CutoffDataPoint = 3, InforColumn = 4, OutputDir = "BridgeR_fig", OutputFile = "BridgeR_4_HalfLife_Pvalue.txt")
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
  \item{OutputFile}{
text and fig files
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
#group <- c("Control","knockdown1","knockdown2")
#hour <- c(0,1,2,4,8,12)
#compfile <- c("Control","Knockdown1")

#BridgeRDrawFittingCurve(filename=inputfile, group=group, hour=hour, ComparisonFile=compfile)

## The function is currently defined as
function (filename, group, hour, ComparisonFile, CutoffRelExp = 0.1, 
    CutoffDataPoint = 3, InforColumn = 4, OutputDir = "BridgeR_fig", 
    OutputFile = "BridgeR_4_HalfLife_Pvalue.txt") 
{
    library(data.table)
    library(ggplot2)
    ComparisonFile_name = paste(ComparisonFile, collapse = "_")
    output_dir_name <- paste(OutputDir, ComparisonFile_name, 
        sep = "_")
    dir.create(output_dir_name)
    time_points <- length(hour)
    group_number <- length(group)
    input_file <- fread(filename, header = T)
    comp_file_number <- NULL
    for (a in 1:length(ComparisonFile)) {
        comp_file_number <- append(comp_file_number, which(group == 
            ComparisonFile[a]))
    }
    output_file <- OutputFile
    setwd(output_dir_name)
    cat("", file = output_file)
    hour_label <- NULL
    for (a in comp_file_number) {
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
            "SD_ori", "half_exp_minus", "half_exp_plus", "half_life_SD", 
            sep = "\t", file = output_file, append = T)
    }
    cat("\t", sep = "", file = output_file, append = T)
    cat("p_value(Welch Modified Two-Sample t-Test)", "\n", sep = "\t", 
        file = output_file, append = T)
    gene_number <- length(input_file[[1]])
    for (x in 1:gene_number) {
        data <- as.vector(as.matrix(input_file[x, ]))
        gene_name <- as.character(data[2])
        file_name <- sprintf("\%1$s.png", gene_name)
        paste(output_dir_name, "/", file_name, sep = "")
        png(filename = file_name, width = 640, height = 640)
        p.fitting <- ggplot()
        flg <- 0
        fig_color <- NULL
        N_for_p <- NULL
        SD_for_p <- NULL
        X_for_p <- NULL
        flg_for_p <- 0
        for (a in comp_file_number) {
            if (flg == 0) {
                fig_color <- "black"
            }
            else {
                fig_color <- "red"
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
            time_point_exp_original <- data.frame(hour, exp)
            p.fitting <- p.fitting + layer(data = time_point_exp_original, 
                mapping = aes(x = hour, y = exp), geom = "point", 
                size = 4, shape = 19, colour = fig_color)
            time_point_exp <- time_point_exp_original[time_point_exp_original$exp >= 
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
                  if (coef < 0) {
                    half_life <- Inf
                  }
                  cat("Exponential_Decay_Model", coef, coef_error, 
                    coef_p, r_squared, adj_r_squared, residual_standard_err, 
                    half_life, sep = "\t", file = output_file, 
                    append = T)
                  half_life_exp_lm <- exp(log(0.5))
                  xmin <- min(hour[1])
                  xmax <- max(hour[length(hour)])
                  predicted2 <- data.frame(hour = time_point_exp$hour)
                  predicted2_ribbon <- data.frame(hour = time_point_exp$hour)
                  pred_conf <- predict(model, predicted2, interval = "prediction", 
                    level = 0.95)
                  pred_conf2_SE <- predict(model, predicted2, 
                    se.fit = T)
                  Fit <- pred_conf2_SE$fit
                  df <- pred_conf2_SE$df
                  SE_fit <- pred_conf2_SE$se.fit
                  Residual_fit <- pred_conf2_SE$residual.scale
                  Space_minus <- function(Fit, SE_fit) {
                    Fit - sqrt(SE_fit^2 + Residual_fit^2) * qt(0.95, 
                      df)
                  }
                  Space_plus <- function(Fit, SE_fit) {
                    Fit + sqrt(SE_fit^2 + Residual_fit^2) * qt(0.95, 
                      df)
                  }
                  SE_test <- function(Fit, SE_fit) {
                    sqrt(SE_fit^2 + Residual_fit^2)
                  }
                  test_minus <- Space_minus(Fit, SE_fit)
                  test_plus <- Space_plus(Fit, SE_fit)
                  test_SE <- SE_test(Fit, SE_fit)
                  SE_table <- data.frame(hour = time_point_exp$hour, 
                    exp = SE_fit)
                  SE_model <- lm((SE_table$exp) ~ SE_table$hour - 
                    1)
                  SE_model_coef <- (summary(SE_model))$coefficients[1]
                  SE_half_life <- SE_model_coef * half_life
                  SE_fitting_curve <- sqrt(SE_half_life^2 + Residual_fit^2)
                  SD_fitting_curve <- SE_fitting_curve * sqrt(data_point - 
                    1)
                  half_life_exp_lm_minus <- exp(log(0.5) - SD_fitting_curve)
                  half_life_exp_lm_plus <- exp(log(0.5) + SD_fitting_curve)
                  half_life_plus <- -log(half_life_exp_lm_minus)/coef
                  half_life_minus <- -log(half_life_exp_lm_plus)/coef
                  SD_for_test <- half_life_plus - half_life
                  SD_for_test <- half_life - half_life_minus
                  N_for_p <- append(N_for_p, data_point)
                  SD_for_p <- append(SD_for_p, SD_for_test)
                  X_for_p <- append(X_for_p, half_life)
                  cat("\t", sep = "\t", file = output_file, append = T)
                  cat(SD_fitting_curve, half_life_exp_lm_minus, 
                    half_life_exp_lm_plus, SD_for_test, sep = "\t", 
                    file = output_file, append = T)
                  predicted2$exp <- exp(as.vector(as.matrix(pred_conf[, 
                    1])))
                  predicted2_ribbon$exp_minus <- exp(as.vector(as.matrix(pred_conf[, 
                    2])))
                  predicted2_ribbon$exp_plus <- exp(as.vector(as.matrix(pred_conf[, 
                    3])))
                  p.fitting <- p.fitting + layer(data = predicted2, 
                    mapping = (aes(x = hour, y = exp)), geom = "line", 
                    size = 1.2, colour = fig_color)
                  p.fitting <- p.fitting + layer(data = predicted2_ribbon, 
                    mapping = aes(x = hour, ymin = exp_minus, 
                      ymax = exp_plus), geom = "ribbon", alpha = 0.1, 
                    fill = fig_color)
                  p.fitting <- p.fitting + ggtitle(gene_name)
                  p.fitting <- p.fitting + xlab("Time")
                  p.fitting <- p.fitting + ylab("Relative RPKM (Time0 = 1)")
                  p.fitting <- p.fitting + xlim(0, 12)
                  ybreaks <- seq(0, 10, 0.1)[2:101]
                  p.fitting <- p.fitting + scale_y_log10(breaks = ybreaks, 
                    labels = ybreaks)
                  plot(p.fitting)
                }
                else {
                  cat("few_data", "NA", "NA", "NA", "NA", "NA", 
                    "NA", "NA", "NA", "NA", "NA", "NA", sep = "\t", 
                    file = output_file, append = T)
                  flg_for_p <- 1
                }
            }
            else {
                cat("low_expresion", "NA", "NA", "NA", "NA", 
                  "NA", "NA", "NA", "NA", "NA", "NA", "NA", sep = "\t", 
                  file = output_file, append = T)
                flg_for_p <- 1
            }
            flg = 1
        }
        t_test <- "NA"
        p_value <- "NA"
        if (flg_for_p != 1) {
            t_test <- tsum.test(mean.x = X_for_p[1], s.x = SD_for_p[1], 
                n.x = N_for_p[1], mean.y = X_for_p[2], s.y = SD_for_p[2], 
                n.y = N_for_p[2])
            p_value <- t_test$p.value
        }
        cat("\t", sep = "\t", file = output_file, append = T)
        cat(p_value, "\n", sep = "\t", file = output_file, append = T)
        dev.off()
        plot.new()
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ BRIC-seq }
\keyword{ Alpha_version }% __ONLY ONE__ keyword per line
