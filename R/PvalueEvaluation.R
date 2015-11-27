#Title: P-value evaluation
#Auther: Naoto Imamachi
#ver: 1.0.0
#Date: 2015-11-20

BridgeRPvalueEvaluation <- function(InputFile="BridgeR_5C_HalfLife_calculation_R2_selection.txt",
                                    group,
                                    hour,
                                    ComparisonFile,
                                    InforColumn=4,
                                    CutoffDataPointNumber = 4,
                                    OutputFile="BridgeR_6_HalfLife_Pvalue_estimation.txt",
                                    Calibration=F){
    ###Prepare_file_infor###
    TimeDataSize <- length(hour)
    group_number <- length(group)
    input_file <- fread(InputFile,header=T)
    output_file <- OutputFile
    comp_file_number <- NULL
    for(index in 1:length(ComparisonFile)){
        comp_file_number <- append(comp_file_number, which(group == ComparisonFile[index]))
    }
    
    ###Print_header###
    cat("", file=output_file) #Reset output file
    halflife_table <- NULL
    correction_value <- NULL
    hour_label <- NULL
    
    for(SampleNum in comp_file_number){
        if(!is.null(hour_label)){
            cat("\t", file=output_file, append=T)
        }
        hour_label <- NULL #Reset hour_label
        for(TimePoints in hour){
            hour_label <- append(hour_label, paste("T",TimePoints,"_",SampleNum, sep=""))
        }
        infor_st <- 1 + (SampleNum - 1)*(TimeDataSize + InforColumn + 3) #3 => Model, R2, half_life
        infor_ed <- (InforColumn)*SampleNum + (SampleNum - 1)*(TimeDataSize + 3) #3 => Model, R2, half_life
        infor_label <- colnames(input_file)[infor_st:infor_ed]
        cat(infor_label,hour_label, sep="\t", file=output_file, append=T)
        cat("\t", sep="", file=output_file, append=T)
        cat("Model","Decay_rate_coef","coef_error","coef_p-value","R2","Adjusted_R2","Residual_standard_error","half_life","half_life_ori",
            "halflife_residual_minus","halflife_residual_minus_R2","half_exp_plus","halflife_residual_plus_R2",
            sep="\t", file=output_file, append=T)
        
        if(Calibration == T){
            halflife_index <- infor_ed + TimeDataSize + 3
            halflife <- input_file[[halflife_index]]
            halflife_table <- cbind(halflife_table, halflife)
        }
    }
    cat("\t", sep="", file=output_file, append=T)
    comp_label <- paste("log2(Relative half-life[",ComparisonFile[2],"/",ComparisonFile[1],"])", sep="")
    cat(comp_label,"p-value(Welch Modified Two-Sample t-test)", sep="\t", file=output_file, append=T)
    cat("\n", sep="", file=output_file, append=T)
    
    if(Calibration == T){
        colnames(halflife_table) <- c("x","y")
        halflife_table <- as.data.frame(halflife_table)
        halflife_table <- halflife_table[halflife_table$x < 24,]
        halflife_table <- halflife_table[halflife_table$y < 24,]
        test_lm <- lm(y ~ x + 0, data=halflife_table)
        correction_value <- round(as.numeric(test_lm$coefficients),digits=3)
        print(paste("Correction value: ",correction_value, sep=""))
    }
    
    ###Calc_p-value###
    pvalue_calc <- function(DataFrame){
        data_point <- length(DataFrame$exp)
        if(!is.null(DataFrame)){
            if(data_point >= CutoffDataPointNumber){
                if(as.numeric(as.vector(as.matrix(DataFrame$exp[1]))) > 0){
                    fitted_model <- lm(log(DataFrame$exp) ~ DataFrame$hour - 1)
                    model_summary <- summary(fitted_model)
                    coef <- -(model_summary$coefficients[1])
                    coef_error <- model_summary$coefficients[2]
                    coef_p <- model_summary$coefficients[4]
                    r_squared <- model_summary$r.squared
                    adj_r_squared <- model_summary$adj.r.squared
                    residual_standard_err <- model_summary$sigma
                    
                    half_life <- log(2)/coef
                    half_life_w <- half_life
                    halflife_value <- log(0.5)
                    
                    if(Calibration == T){
                        if(flg == 1){
                            half_life_w <- half_life/correction_value
                            half_life <- half_life/correction_value
                            halflife_value <- -coef*half_life
                        }
                    }
                    
                    if(coef < 0){
                        half_life <- Inf
                        half_life_w <- 24
                    }else if(half_life > 24){
                        half_life_w <- 24
                    }
                    cat(model,coef,coef_error,coef_p,
                        r_squared,adj_r_squared,residual_standard_err,
                        half_life_w,half_life,
                        sep="\t", file=output_file, append=T)
                    
                    predict_data <- data.frame(hour=DataFrame$hour)
                    pred_conf <- predict(fitted_model, predict_data, se.fit=T)
                    
                    predicted_exp <- pred_conf$fit
                    predicted_exp_SE <- pred_conf$se.fit #Residual for model??
                    predicted_exp_df <- pred_conf$df #Degree of freedom
                    predicted_exp_residual <- pred_conf$residual.scale #Residual for exp data??
                    
                    #residual <- sqrt(predicted_exp_SE^2 + predicted_exp_residual^2) *qt(0.975,df) #95% Prediction interval
                    SE_residual <- sqrt(predicted_exp_SE^2 + predicted_exp_residual^2) #SE
                    SD_residual <- SE_residual * sqrt(predicted_exp_df)
                    predicted_exp_lwr <- predicted_exp - SD_residual
                    predicted_exp_upr <- predicted_exp + SD_residual
                    
                    predicted_exp_lwr_table <- data.frame(hour=DataFrame$hour, exp=predicted_exp_lwr)
                    predicted_exp_lwr_model <- lm(predicted_exp_lwr_table$exp ~ predicted_exp_lwr_table$hour)
                    predicted_exp_lwr_model_summary <- summary(predicted_exp_lwr_model)
                    predicted_exp_lwr_coef <- predicted_exp_lwr_model_summary$coefficients[2]
                    predicted_exp_lwr_intercept <- predicted_exp_lwr_model_summary$coefficients[1]
                    predicted_exp_lwr_R2 <- predicted_exp_lwr_model_summary$r.squared #
                    predicted_exp_lwr_halflife <- (halflife_value - predicted_exp_lwr_intercept)/predicted_exp_lwr_coef
                    
                    predicted_exp_upr_table <- data.frame(hour=DataFrame$hour, exp=predicted_exp_upr)
                    predicted_exp_upr_model <- lm(predicted_exp_upr_table$exp ~ predicted_exp_upr_table$hour)
                    predicted_exp_upr_model_summary <- summary(predicted_exp_upr_model)
                    predicted_exp_upr_coef <- predicted_exp_upr_model_summary$coefficients[2]
                    predicted_exp_upr_intercept <- predicted_exp_upr_model_summary$coefficients[1]
                    predicted_exp_upr_R2 <- predicted_exp_upr_model_summary$r.squared #
                    predicted_exp_upr_halflife <- (halflife_value - predicted_exp_upr_intercept)/predicted_exp_upr_coef
                    
                    halflife_SD_minus <- half_life - predicted_exp_lwr_halflife #
                    halflife_SD_plus <- predicted_exp_upr_halflife - half_life #
                    
                    cat("\t", sep="", file=output_file, append=T)
                    cat(halflife_SD_minus,predicted_exp_lwr_R2,halflife_SD_plus,predicted_exp_upr_R2,
                        sep="\t", file=output_file, append=T)
                    #Half-life, The degree of freedom, SD_minus, SD_plus
                    return(c(half_life_w, half_life, predicted_exp_df, halflife_SD_minus, halflife_SD_plus))
                }
            }
        }
    }
    
    gene_number <- length(input_file[[1]])
    for(index in 1:gene_number){
        data <- as.vector(as.matrix(input_file[index,]))
        
        flg <- 0
        halflife_w <- NULL
        Mean_for_p <- NULL
        N_for_p <- NULL
        SD_minus_for_p <- NULL
        SD_plus_for_p <- NULL
        flg_for_p <- 0
        for(SampleNum in comp_file_number){
            if(flg == 1){
                cat("\t", sep="", file=output_file, append=T)
            }
            infor_st <- 1 + (SampleNum - 1)*(TimeDataSize + InforColumn + 3)
            infor_ed <- (InforColumn)*SampleNum + (SampleNum - 1)*(TimeDataSize + 3)
            exp_st <- infor_ed + 1
            exp_ed <- infor_ed + TimeDataSize
            
            gene_infor <- data[infor_st:infor_ed]
            exp <- as.numeric(data[exp_st:exp_ed])
            
            cat(gene_infor, sep="\t", file=output_file, append=T)
            cat("\t", file=output_file, append=T)
            cat(exp, sep="\t", file=output_file, append=T)
            cat("\t", file=output_file, append=T)
            
            TimePoint_Exp_data_raw <- data.frame(hour,exp)
            
            model_index <- exp_ed + 1 #1: Model, #2:R2, #3:half_life
            model <- data[model_index]
            
            ###Test###
            ttest_infor <- NULL
            if(model == "Notest"){
                cat("NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA", sep="\t", file=output_file, append=T)
                flg_for_p <- 1
            }else if(model == "Raw"){
                TimePoint_Exp_data_base <- TimePoint_Exp_data_raw[TimePoint_Exp_data_raw$exp > 0,]
                ttest_infor <- pvalue_calc(TimePoint_Exp_data_base)
            }else{
                check <- gsub("Delete_","",model)
                check <- gsub("hr","",check)
                check <- as.numeric(strsplit(check,"_")[[1]])
                TimePoint_Exp_data_del <- TimePoint_Exp_data_raw
                for(times_list in check){
                    TimePoint_Exp_data_del <- TimePoint_Exp_data_del[TimePoint_Exp_data_del$hour != times_list,]
                }
                TimePoint_Exp_data_del <- TimePoint_Exp_data_del[TimePoint_Exp_data_del$exp > 0,]
                ttest_infor <- pvalue_calc(TimePoint_Exp_data_del)
            }
            
            halflife_w <- append(halflife_w, ttest_infor[1])
            Mean_for_p <- append(Mean_for_p, ttest_infor[2])
            N_for_p <- append(N_for_p, ttest_infor[3])
            SD_minus_for_p <- append(SD_minus_for_p, ttest_infor[4])
            SD_plus_for_p <- append(SD_plus_for_p, ttest_infor[5])
            
            flg <- 1
        }
        
        halflife_comp <- "NA"
        p_value <- "NA"
        if(flg_for_p != 1){
            if(Mean_for_p[1] <= Mean_for_p[2]){
                cond1_halflife_mean <- Mean_for_p[1]
                cond2_halflife_mean <- Mean_for_p[2]
                cond1_df <- N_for_p[1]
                cond2_df <- N_for_p[2]
                cond1_halflife_SD <- SD_plus_for_p[1] #check!
                cond2_halflife_SD <- SD_minus_for_p[2] #check!
                
                t_test <- tsum.test(mean.x=cond1_halflife_mean, n.x=cond1_df, s.x=cond1_halflife_SD,
                                    mean.y=cond2_halflife_mean, n.y=cond2_df, s.y=cond2_halflife_SD)
                p_value <- t_test$p.value
            }else if(Mean_for_p[1] > Mean_for_p[2]){
                cond1_halflife_mean <- Mean_for_p[1]
                cond2_halflife_mean <- Mean_for_p[2]
                cond1_df <- N_for_p[1]
                cond2_df <- N_for_p[2]
                cond1_halflife_SD <- SD_minus_for_p[1] #check!
                cond2_halflife_SD <- SD_plus_for_p[2] #check!
                
                t_test <- tsum.test(mean.x=cond1_halflife_mean, n.x=cond1_df, s.x=cond1_halflife_SD,
                                    mean.y=cond2_halflife_mean, n.y=cond2_df, s.y=cond2_halflife_SD)
                p_value <- t_test$p.value
            }
            halflife_comp <- log2(halflife_w[2]/halflife_w[1])
        }
        cat("\t", sep="", file=output_file, append=T)
        cat(halflife_comp,p_value, sep="\t", file=output_file, append=T)
        cat("\n", sep="", file=output_file, append=T)
    }
}


###TEST##################
#test <- gsub("Delete_","","Delete_1hr_2hr")
#test <- gsub("hr","",test)
#test <- as.numeric(strsplit(test,"_")[[1]])


#hour <- c(0,1,2,4,8,12)
#exp <- c(1,0.8144329,0.8173235,0.7693528,0.4247929,0.2367443) #
#hour <- c(0,1,2,4)
#exp <- c(1,0.1359368,0.1165081,0.1111112)
#testing_data <- data.frame(hour=hour,exp=exp)
#pred_data <- data.frame(hour=testing_data$hour)
#fitted_model <- lm(log(testing_data$exp) ~ testing_data$hour - 1)
#fitted_model_coef <- -((summary(fitted_model)))$coefficients[1]
#fitted_model_halflife <- log(2)/fitted_model_coef

#pred_conf <- predict(fitted_model, pred_data, se.fit=T)
#predicted_exp <- pred_conf$fit
#predicted_exp_SE <- pred_conf$se.fit #Residual for model??
#predicted_exp_df <- pred_conf$df #Degree of freedom
#predicted_exp_residual <- pred_conf$residual.scale #Residual for exp data??
#SE_residual <- sqrt(predicted_exp_SE^2 + predicted_exp_residual^2) #SE
#SD_residual <- SE_residual * sqrt(predicted_exp_df)
#predicted_exp_lwr <- predicted_exp - SD_residual
#predicted_exp_upr <- predicted_exp + SD_residual

#predicted_exp_lwr_table <- data.frame(hour=hour, exp=predicted_exp_lwr)
#predicted_exp_lwr_model <- lm(predicted_exp_lwr_table$exp ~ predicted_exp_lwr_table$hour)
#predicted_exp_lwr_coef <- ((summary(predicted_exp_lwr_model)))$coefficients[2]
#predicted_exp_lwr_intercept <- ((summary(predicted_exp_lwr_model)))$coefficients[1]
#predicted_exp_lwr_halflife <- (-log(2) - predicted_exp_lwr_intercept)/predicted_exp_lwr_coef

#predicted_exp_upr_table <- data.frame(hour=hour, exp=predicted_exp_upr)
#predicted_exp_upr_model <- lm(predicted_exp_upr_table$exp ~ predicted_exp_upr_table$hour)
#predicted_exp_upr_coef <- ((summary(predicted_exp_upr_model)))$coefficients[2]
#predicted_exp_upr_intercept <- ((summary(predicted_exp_upr_model)))$coefficients[1]
#predicted_exp_upr_halflife <- (-log(2) - predicted_exp_upr_intercept)/predicted_exp_upr_coef



#SD_table <- data.frame(hour=hour, exp=SD_residual)
#SD_model <- lm(SD_table$exp ~ SD_table$hour)
 
#predicted_exp_table <- data.frame(hour=hour, exp=predicted_exp)
#predicted_exp_model <- lm(predicted_exp_table$exp ~ predicted_exp_table$hour - 1)
#predicted_exp_coef <- ((summary(predicted_exp_model)))$coefficients[1]
#predicted_exp_halflife <- log(0.5)/predicted_exp_coef

