#Title: Estimate fitting decay curve
#Auther: Naoto Imamachi
#ver: 1.0.0
#Date: 2015-10-24

###Estimate_normalization_factor_function###
BridgeRHalfLifeCalcR2Select <- function(InputFile = "BridgeR_4_Normalized_expression_dataset.txt",
                                        group, 
                                        hour, 
                                        InforColumn = 4, 
                                        CutoffDataPointNumber = 4,
                                        CutoffDataPoint1 = c(1,2),
                                        CutoffDataPoint2 = c(8,12),
                                        ThresholdHalfLife = c(8,12),
                                        OutputFile = "BridgeR_5C_HalfLife_calculation_R2_selection.txt"){
    ###Prepare_file_infor###
    time_points <- length(hour)
    group_number <- length(group)
    input_file <- fread(InputFile, header=T)
    output_file <- paste(OutputFile,".log",sep="")
    log_file <- OutputFile
    
    ###print_header###
    cat("",file=output_file)
    cat("",file=log_file)
    hour_label <- NULL
    CutoffDataPoint1 <- sort(CutoffDataPoint1, decreasing = T)
    CutoffDataPoint2 <- sort(CutoffDataPoint2, decreasing = F)
    CutoffDataPoint1_length <- length(CutoffDataPoint1)
    CutoffDataPoint2_length <- length(CutoffDataPoint2)
    
    for_iteration_number <- CutoffDataPoint1_length + CutoffDataPoint2_length
    for(a in 1:group_number){
        if(!is.null(hour_label)){
            cat("\t", file=output_file, append=T)
            cat("\t", file=log_file, append=T)
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
        cat(infor,hour_label, sep="\t", file=log_file, append=T)
        cat("\t", sep="", file=log_file, append=T)
        
        for_iteration_number2 <- for_iteration_number + 1
        for(number in 1:for_iteration_number2){
            cat("Model","Decay_rate_coef","coef_error","coef_p-value",
                "R2","Adjusted_R2","Residual_standard_error","half_life",
                sep="\t", file=output_file, append=T)
            cat("\t", sep="", file=output_file, append=T)
        }
        cat("Model","R2","half_life", sep="\t", file=output_file, append=T)
        cat("Model","R2","half_life", sep="\t", file=log_file, append=T)
    }
    cat("\n", sep="", file=output_file, append=T)
    cat("\n", sep="", file=log_file, append=T)
    
    ###calc_RNA_half_lives###
    ###Function1#################################
    half_calc <- function(time_exp_table,label){
        data_point <- length(time_exp_table$exp)
        if(!is.null(time_exp_table)){
            if(data_point >= CutoffDataPointNumber){
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
                    cat(label,coef,coef_error,coef_p,r_squared,adj_r_squared,residual_standard_err,half_life, sep="\t", file=output_file, append=T)
                    return(c(half_life,r_squared))
                }else{
                    cat("low_expresion","NA","NA","NA","NA","NA","NA","NA", sep="\t", file=output_file, append=T)
                    return(c("NA","NA"))
                }
            }else{
                cat("few_data","NA","NA","NA","NA","NA","NA","NA", sep="\t", file=output_file, append=T)
                return(c("NA","NA"))
            }
        }else{
            cat("low_expresion","NA","NA","NA","NA","NA","NA","NA", sep="\t", file=output_file, append=T)
            return(c("NA","NA"))
        }
    }
    #############################################
    
    ###Function2#################################
    test_R2 <- function(time_point_exp_raw, cutoff_data_point, halflife_Raw, R2_Raw){
        test_times <- cutoff_data_point
        times_length <- length(test_times)
        times_index <- c(times_length)
        add_index <- times_length
        
        R2_list <- c(R2_Raw)
        half_list <- c(halflife_Raw)
        label_list <- c("Raw")
        for(counter in times_length:1){
            #excepted_time_points
            check_times <- test_times[times_index] #c(24), c(24,12), c(24,12,8)
            time_point_exp_del <- NULL
            time_point_exp_del_label <- paste("Delete_",paste(check_times,collapse="hr_"),"hr",sep="")
            label_list <- append(label_list, time_point_exp_del_label) #
            time_point_exp_del <- time_point_exp_raw
            for(times_list in check_times){
                time_point_exp_del <- time_point_exp_del[time_point_exp_del$hour != as.numeric(times_list),]
            }
            time_point_exp_del <- time_point_exp_del[time_point_exp_del$exp > 0,]
            halflife_R2_result <- half_calc(time_point_exp_del, time_point_exp_del_label)
            cat("\t", sep="", file=output_file, append=T)
            
            R2_list <- append(R2_list, halflife_R2_result[2]) #
            half_list <- append(half_list, halflife_R2_result[1]) #
            
            #Counter
            add_index <- add_index - 1
            times_index <- append(times_index, add_index)
        }
        
        R2_table <- data.frame(label=label_list, R2=R2_list, half=half_list)
        #R2_table <- R2_table[R2_table$R2 != "NA",]
        #sortlist <- order(R2_table$R2, decreasing = T)
        #R2_table <- R2_table[sortlist,]
        #result <- as.vector(as.matrix(R2_table[1,]))
        #cat(result, sep="\t", file=output_file, append=T)
        return(R2_table)
    }
    #############################################
    
    gene_number <- length(input_file[[1]])
    for(x in 1:gene_number){
        data <- as.vector(as.matrix(input_file[x,]))
        for(a in 1:group_number){
            if(a != 1){
                cat("\t", sep="", file=output_file, append=T)
                cat("\t", sep="", file=log_file, append=T)
            }
            infor_st <- 1 + (a - 1)*(time_points + InforColumn)
            infor_ed <- (InforColumn)*a + (a - 1)*time_points
            exp_st <- infor_ed + 1
            exp_ed <- infor_ed + time_points
            
            gene_infor <- data[infor_st:infor_ed]
            cat(gene_infor, sep="\t", file=output_file, append=T)
            cat("\t", file=output_file, append=T)
            cat(gene_infor, sep="\t", file=log_file, append=T)
            cat("\t", file=log_file, append=T)
            
            exp <- as.numeric(data[exp_st:exp_ed])
            cat(exp, sep="\t", file=output_file, append=T)
            cat("\t", file=output_file, append=T)
            cat(exp, sep="\t", file=log_file, append=T)
            cat("\t", file=log_file, append=T)
            
            ###Raw_data
            time_point_exp_raw <- data.frame(hour,exp)
            time_point_exp_base <- time_point_exp_raw[time_point_exp_raw$exp > 0, ]
            test <- half_calc(time_point_exp_base,"Exponential_Decay_Model") #T1/2 - R2
            cat("\t", sep="", file=output_file, append=T)
            
            ###Re-calculation of RNA half-life(Default: ThresholdHalfLife - 12hr)
            R2_list <- NULL
            half_list <- NULL
            label_list <- NULL
            if(test[1] == "NA"){
                for(number in 1:for_iteration_number){
                    cat("Notest","NA","NA","NA","NA","NA","NA","NA", sep="\t", file=output_file, append=T)
                    cat("\t", sep="", file=output_file, append=T)
                }
                cat("Notest","NA","NA", sep="\t", file=output_file, append=T)
                cat("Notest","NA","NA", sep="\t", file=log_file, append=T)
            }else if(test[1] < ThresholdHalfLife[1]){ #Default: <12hr => Delete 8,12hr
                for(number in 1:CutoffDataPoint1_length){
                    cat("Notest","NA","NA","NA","NA","NA","NA","NA", sep="\t", file=output_file, append=T)
                    cat("\t", sep="", file=output_file, append=T)
                }
                R2_table <- test_R2(time_point_exp_raw, CutoffDataPoint2, test[1], test[2]) #test[1], test[2] => T1/2, R2
                
            }else if(test[1] >= ThresholdHalfLife[1] && test[1] < ThresholdHalfLife[2]){
                R2_table1 <- test_R2(time_point_exp_raw, CutoffDataPoint1, test[1], test[2]) #test[1], test[2] => T1/2, R2
                R2_table2 <- test_R2(time_point_exp_raw, CutoffDataPoint2, test[1], test[2]) #test[1], test[2] => T1/2, R2
                R2_table <- rbind(R2_table1, R2_table2)
                
            }else if(test[1] >= ThresholdHalfLife[2]){ #Default: >=12hr => Delete 1,2hr
                for(number in 1:CutoffDataPoint2_length){
                    cat("Notest","NA","NA","NA","NA","NA","NA","NA", sep="\t", file=output_file, append=T)
                    cat("\t", sep="", file=output_file, append=T)
                }
                R2_table <- test_R2(time_point_exp_raw, CutoffDataPoint1, test[1], test[2]) #test[1], test[2] => T1/2, R2
                
            }
            ###R2_Selection###
            if(test[1] != "NA"){
                R2_table <- R2_table[R2_table$R2 != "NA",]
                sortlist <- order(R2_table$R2, decreasing = T)
                R2_table <- R2_table[sortlist,]
                result <- as.vector(as.matrix(R2_table[1,]))
                cat(result, sep="\t", file=output_file, append=T)
                cat(result, sep="\t", file=log_file, append=T)
            }
        }
        cat("\n", file=output_file, append=T)
        cat("\n", file=log_file, append=T)
    }
}

###TEST###
#output_file <- "test.txt"
#CutoffDataPointNumber <- 4
#x <- c(0,1,2,4,8,12)
#y <- c(1,0.9,0.8,0.5,0.3,0.1)
#table <- data.frame(hour=x,exp=y)
#test_data_point <- c(8,12)
#test_R2(table, test_data_point, 4.5, 0.992)

#    R2_list half_list      label_list
#1 0.9919170  4.547501     Delete_12hr
#2 0.9720505  4.378777 Delete_12hr_8hr

#########################
#test <- c(8,12,24,2,12,21,1)
#test_length <- length(test)
#test_index <- c(test_length)
#add_index <- test_length

#for(x in test_length:1){
#    print(test[test_index])
#    add_index <- add_index - 1
#    test_index <- append(test_index,add_index)
#}
