#Title: Estimate fitting decay curve
#Auther: Naoto Imamachi
#ver: 1.0.0
#Date: 2015-10-08

###Estimate_normalization_factor_function###
BridgeRHalfLifeCalcR2Select <- function(filename = "BridgeR_4_Normalized_expression_data.txt", group, hour, InforColumn = 4, CutoffRelExp = 0.1, CutoffDataPoint = 3, ExceptTimePoint=c(1), OutputFile = "BridgeR_5C_HalfLife_calculation_R2Selection.txt"){
    ###Prepare_file_infor###
    time_points <- length(hour)
    group_number <- length(group)
    input_file <- fread(filename, header=T)
    output_file <- OutputFile
    
    ###print_header###
    cat("",file=output_file)
    hour_label <- NULL
    for(a in 1:group_number){
        if(!is.null(hour_label)){
            cat("\t", file=output_file, append=T)
        }
        hour_label <- NULL
        for(x in hour){
            hour_label <- append(hour_label, paste("T", x, "_", a, sep=""))
        }
        infor_st <- 1 + (a - 1)*(time_points + InforColumn)
        infor_ed <- (InforColumn)*a + (a - 1)*time_points
        infor <- colnames(input_file)[infor_st:infor_ed]
        cat(infor,hour_label, sep="\t", file=output_file, append=T)
        cat("\t", sep="", file=output_file, append=T)
        cat("Model","Decay_rate_coef","coef_error","coef_p-value","R2","Adjusted_R2","Residual_standard_error","half_life", sep="\t", file=output_file, append=T)
        cat("\t", sep="", file=output_file, append=T)
        cat("Model","Decay_rate_coef","coef_error","coef_p-value","R2","Adjusted_R2","Residual_standard_error","half_life", sep="\t", file=output_file, append=T)
        cat("\t", sep="", file=output_file, append=T)
        cat("Model","Decay_rate_coef","coef_error","coef_p-value","R2","Adjusted_R2","Residual_standard_error","half_life", sep="\t", file=output_file, append=T)
    }
    cat("\n", sep="", file=output_file, append=T)
    
    ###calc_RNA_half_lives###
    half_calc <- function(time_exp_table){
        data_point <- length(time_exp_table$exp)
        if(!is.null(time_exp_table)){
            if(data_point >= CutoffDataPoint){
                if(as.numeric(as.vector(as.matrix(time_exp_table$exp[1]))) > 0){
                    model <- lm(log(time_exp_table$exp) ~ time_exp_table$hour - 1)
                    model_summary <- summary(model)
                    coef <- -model_summary$coefficients[1]
                    coef_error <- model_summary$coefficients[2]
                    coef_p <- model_summary$coefficients[4]
                    r_squared <- model_summary$r.squared
                    adj_r_squared <- model_summary$adj.r.squared
                    residual_standard_err <- model_summary$sigma
                    half_life <- log(2)/coef
                    if(coef < 0 || half_life >= 24){
                        half_life <- 24
                    }
                    #if(coef < 0){
                    #    half_life <- Inf
                    #}
                    cat("Exponential_Decay_Model",coef,coef_error,coef_p,r_squared,adj_r_squared,residual_standard_err,half_life, sep="\t", file=output_file, append=T)
                }else{
                    cat("low_expresion","NA","NA","NA","NA","NA","NA","NA", sep="\t", file=output_file, append=T)
                }
            }else{
                cat("few_data","NA","NA","NA","NA","NA","NA","NA", sep="\t", file=output_file, append=T)
            }
        }else{
            cat("low_expresion","NA","NA","NA","NA","NA","NA","NA", sep="\t", file=output_file, append=T)
        }
    }
    
    gene_number <- length(input_file[[1]])
    for(x in 1:gene_number){
        data <- as.vector(as.matrix(input_file[x,]))
        for(a in 1:group_number){
            if(a != 1){
                cat("\t", sep="", file=output_file, append=T)
            }
            infor_st <- 1 + (a - 1)*(time_points + InforColumn)
            infor_ed <- (InforColumn)*a + (a - 1)*time_points
            exp_st <- infor_ed + 1
            exp_ed <- infor_ed + time_points
            
            gene_infor <- data[infor_st:infor_ed]
            cat(gene_infor, sep="\t", file=output_file, append=T)
            cat("\t", file=output_file, append=T)
            exp <- as.numeric(data[exp_st:exp_ed])
            cat(exp, sep="\t", file=output_file, append=T)
            cat("\t", file=output_file, append=T)
            
            ###Raw_data
            time_point_exp_raw <- data.frame(hour,exp)
            time_point_exp_base <- time_point_exp_raw[time_point_exp_raw$exp > 0, ]
            half_calc(time_point_exp_base)
            cat("\t", sep="", file=output_file, append=T)
            
            ###Except_time_point(Default: CutoffExp <= 0.1)
            time_point_exp_cutoff <- time_point_exp_raw[time_point_exp_raw$exp >= CutoffRelExp, ]
            half_calc(time_point_exp_cutoff)
            cat("\t", sep="", file=output_file, append=T)
            
            ###Except_time_point(you decide)
            time_point_exp_exception <- time_point_exp_raw #Copy
            for(tp in 1:length(ExceptTimePoint)){
                time_point_exp_exception <- time_point_exp_exception[time_point_exp_exception$hour!=ExceptTimePoint[tp],]
            }
            time_point_exp_exception <- time_point_exp_exception[time_point_exp_exception$exp > 0, ]
            half_calc(time_point_exp_exception)
            
        }
        cat("\n", file=output_file, append=T)
    }
}

