#Title: Estimate fitting decay curve
#Auther: Naoto Imamachi
#ver: 1.0.0
#Date: 2015-10-18

###Estimate_normalization_factor_function###
BridgeRHalfLifeCalc3 <- function(filename = "BridgeR_4_Normalized_expression_dataset_20151118.txt", 
                                 group, 
                                 hour, 
                                 InforColumn = 4, 
                                 ThresholdHalfLife = c(8,12),
                                 CutoffDataPoint = 4, 
                                 OutputFile = "BridgeR_5D_HalfLife_calculation_ver4_20151119.txt"){
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
        cat("\t", sep="", file=output_file, append=T)
        cat("Model","R2","half_life", sep="\t", file=output_file, append=T)
    }
    cat("\n", sep="", file=output_file, append=T)
    
    ###calc_RNA_half_lives###
    half_calc <- function(time_exp_table,label,SaveFile=T){
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
                    if(SaveFile == T){
                        cat(label,coef,coef_error,coef_p,r_squared,adj_r_squared,residual_standard_err,half_life, sep="\t", file=output_file, append=T)
                    }else if(SaveFile == F){
                    }
                    return(c(half_life,r_squared))
                }else{
                    if(SaveFile == T){
                        cat("low_expresion","NA","NA","NA","NA","NA","NA","NA", sep="\t", file=output_file, append=T)
                    }else if(SaveFile == F){
                    }
                    return(c("NA","NA"))
                }
            }else{
                if(SaveFile == T){
                    cat("few_data","NA","NA","NA","NA","NA","NA","NA", sep="\t", file=output_file, append=T)
                }else if(SaveFile == F){
                }
                return(c("NA","NA"))
            }
        }else{
            if(SaveFile == T){
                cat("low_expresion","NA","NA","NA","NA","NA","NA","NA", sep="\t", file=output_file, append=T)
            }else if(SaveFile == F){
            }
            return(c("NA","NA"))
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
            test <- half_calc(time_point_exp_base,"Exponential_Decay_Model")
            cat("\t", sep="", file=output_file, append=T)
            
            ###Re-calculation of RNA half-life(Default: ThresholdHalfLife - 12hr)
            R2_list <- NULL
            half_list <- NULL
            label_list <- NULL
            if(test[1] == "NA"){
                cat("Notest","NA","NA","NA","NA","NA","NA","NA", sep="\t", file=output_file, append=T)
                cat("\t", sep="", file=output_file, append=T)
                cat("Notest","NA","NA","NA","NA","NA","NA","NA", sep="\t", file=output_file, append=T)
                cat("\t", sep="", file=output_file, append=T)
                cat("Notest","NA","NA", sep="\t", file=output_file, append=T)
            }else if(test[1] < ThresholdHalfLife[1]){ #Default: <12hr => Delete 8,12hr
                ###Delete 12hr
                time_point_exp_del_1 <- time_point_exp_raw[time_point_exp_raw$hour != 12,]
                time_point_exp_del_1 <- time_point_exp_del_1[time_point_exp_del_1$exp > 0,]
                test1 <- half_calc(time_point_exp_del_1,"Delete_12hr")
                cat("\t", sep="", file=output_file, append=T)
                ###Delete 8,12hr
                time_point_exp_del_2 <- time_point_exp_raw[time_point_exp_raw$hour != 12,]
                time_point_exp_del_2 <- time_point_exp_del_2[time_point_exp_del_2$hour != 8,]
                time_point_exp_del_2 <- time_point_exp_del_2[time_point_exp_del_2$exp > 0,]
                test2 <- half_calc(time_point_exp_del_2,"Delete_8hr_12hr")
                cat("\t", sep="", file=output_file, append=T)
                ###R2_list
                R2_list <- append(R2_list,c(test[2],test1[2],test2[2]))
                half_list <- append(half_list,c(test[1],test1[1],test2[1]))
                label_list <- c("Raw","Delete_12hr","Delete_8hr_12hr")
            }else if(test[1] >= ThresholdHalfLife[1] && test[1] < ThresholdHalfLife[2]){
                ###Delete 12hr
                time_point_exp_del_1 <- time_point_exp_raw[time_point_exp_raw$hour != 12,]
                time_point_exp_del_1 <- time_point_exp_del_1[time_point_exp_del_1$exp > 0,]
                test1 <- half_calc(time_point_exp_del_1,"Delete_12hr")
                cat("\t", sep="", file=output_file, append=T)
                ###Delete 8,12hr
                time_point_exp_del_2 <- time_point_exp_raw[time_point_exp_raw$hour != 12,]
                time_point_exp_del_2 <- time_point_exp_del_2[time_point_exp_del_2$hour != 8,]
                time_point_exp_del_2 <- time_point_exp_del_2[time_point_exp_del_2$exp > 0,]
                test2 <- half_calc(time_point_exp_del_2,"Delete_8hr_12hr")
                cat("\t", sep="", file=output_file, append=T)
                ###Delete 1hr
                time_point_exp_del_3 <- time_point_exp_raw[time_point_exp_raw$hour != 1,]
                time_point_exp_del_3 <- time_point_exp_del_3[time_point_exp_del_3$exp > 0,]
                test3 <- half_calc(time_point_exp_del_3,"Delete_1hr",SaveFile=F)
                ###Delete 1,2hr
                time_point_exp_del_4 <- time_point_exp_raw[time_point_exp_raw$hour != 1,]
                time_point_exp_del_4 <- time_point_exp_del_4[time_point_exp_del_4$hour != 2,]
                time_point_exp_del_4 <- time_point_exp_del_4[time_point_exp_del_4$exp > 0,]
                test4 <- half_calc(time_point_exp_del_4,"Delete_1hr_2hr",SaveFile=F)
                ###R2_list
                R2_list <- append(R2_list,c(test[2],test1[2],test2[2],test3[2],test4[2]))
                half_list <- append(half_list,c(test[1],test1[1],test2[1],test3[1],test4[1]))
                label_list <- c("Raw","Delete_12hr","Delete_8hr_12hr","Delete_1hr","Delete_1hr_2hr")
            }else if(test[1] >= ThresholdHalfLife[2]){ #Default: >=12hr => Delete 1,2hr
                ###Delete 1hr
                time_point_exp_del_1 <- time_point_exp_raw[time_point_exp_raw$hour != 1,]
                time_point_exp_del_1 <- time_point_exp_del_1[time_point_exp_del_1$exp > 0,]
                test1 <- half_calc(time_point_exp_del_1,"Delete_1hr")
                cat("\t", sep="", file=output_file, append=T)
                ###Delete 1,2hr
                time_point_exp_del_2 <- time_point_exp_raw[time_point_exp_raw$hour != 1,]
                time_point_exp_del_2 <- time_point_exp_del_2[time_point_exp_del_2$hour != 2,]
                time_point_exp_del_2 <- time_point_exp_del_2[time_point_exp_del_2$exp > 0,]
                test2 <- half_calc(time_point_exp_del_2,"Delete_1hr_2hr")
                cat("\t", sep="", file=output_file, append=T)
                ###R2_list
                R2_list <- append(R2_list,c(test[2],test1[2],test2[2]))
                half_list <- append(half_list,c(test[1],test1[1],test2[1]))
                label_list <- c("Raw","Delete_1hr","Delete_1hr_2hr")
            }
            ###R2_Selection###
            if(test[1] != "NA"){
                R2_table <- data.frame(label=label_list,R2=R2_list,half=half_list)
                R2_table <- R2_table[R2_table$R2 != "NA",]
                sortlist <- order(R2_table$R2,decreasing=T)
                R2_table <- R2_table[sortlist,]
                result <- as.vector(as.matrix(R2_table[1,]))
                cat(result, sep="\t", file=output_file, append=T)
            }
        }
        cat("\n", file=output_file, append=T)
    }
}

