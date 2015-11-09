#Title: Estimate fitting decay curve
#Auther: Naoto Imamachi
#ver: 1.0.0
#Date: 2015-10-08

###Estimate_normalization_factor_function###
BridgeRHalfLifeCalcModel3 <- function(filename = "BridgeR_4_Normalized_expression_dataset.txt", group, hour, InforColumn = 4, CutoffRelExp = 0.1, CutoffDataPoint = 3, OutputFile = "BridgeR_5B_HalfLife_calculation_3model.txt"){
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
        cat("status","model1","cor1","a_1","half1_1","half1_2","AIC1","model2","cor2","a_2","b_2","half2",
            "AIC2","model3","cor3","a_3","b_3","c_3","half3","AIC3","selected_model","Selected_R2","Selected_HalfLife",sep="\t",file=output_file, append=T)
    }
    cat("\n", sep="", file=output_file, append=T)
    
    ###calc_RNA_half_lives###
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
            
            #Write_gene_infor
            gene_infor <- data[infor_st:infor_ed]
            cat(gene_infor, sep="\t", file=output_file, append=T)
            cat("\t", file=output_file, append=T)
            
            #Write_time_point_data
            exp <- as.numeric(data[exp_st:exp_ed])
            cat(exp, sep="\t", file=output_file, append=T)
            cat("\t", file=output_file, append=T)
            
            #make_hour-exp_dataframe
            time_point_exp <- data.frame(hour,exp)
            time_point_exp <- time_point_exp[time_point_exp$exp >= CutoffRelExp, ]
            data_point <- length(time_point_exp$exp)
            
            if(!is.null(time_point_exp)){
                if(data_point >= CutoffDataPoint){
                    cat("ok","\t",sep="",file=output_file,append=T)
                    
                    ###Model1: f(t) = exp(-a * t)###
                    optim1 <- function(x){
                        mRNA_exp <- exp(-x * time_point_exp$hour)
                        sum((time_point_exp$exp - mRNA_exp)^2)
                    }
                    
                    out1 <- optim(1,optim1)
                    min1 <- out1$value
                    a_1 <- out1$par[1]
                    half1_1 <- log(2) / a_1
                    
                    model1_pred <- function(x){
                        mRNA_exp <- exp(-x * time_point_exp$hour)
                        (cor(mRNA_exp, time_point_exp$exp, method="pearson"))^2
                    }
                    
                    cor1 <- model1_pred(a_1)
                    
                    model1_half <- function(x){
                        mRNA_half <- exp(-a_1 * x)
                        (mRNA_half - 0.5)^2
                    }
                    
                    out1 <- optim(1,model1_half)
                    half1_2 <- out1$par
                    mRNA_pred <- exp(-a_1 * time_point_exp$hour)
                    s2 <- sum((time_point_exp$exp - mRNA_pred)^2) / data_point
                    AIC1 <- data_point * log(s2) + (2 * 0)
                    
                    cat(min1,cor1,a_1,half1_1,half1_2,AIC1,sep="\t",file=output_file,append=T)
                    cat("\t",file=output_file,append=T)
                    ############################################
                    
                    ###Model2: f(t) = (1 - b)exp(-a * t) + b###
                    optim2 <- function(x){
                        mRNA_exp <- (1.0 - x[2]) * exp(-x[1] * time_point_exp$hour) + x[2]
                        sum((time_point_exp$exp - mRNA_exp)^2)
                    }
                    
                    out2 <- optim(c(1,0),optim2)
                    min2 <- out2$value
                    a_2 <- out2$par[1]
                    b_2 <- out2$par[2]
                    
                    model2_pred <- function(a,b){
                        mRNA_exp <- (1.0 - b) * exp(-a * time_point_exp$hour) + b
                        (cor(mRNA_exp, time_point_exp$exp, method="pearson"))^2
                    }
                    
                    cor2 <- model2_pred(a_2, b_2)
                    
                    model2_half <- function(x){
                        mRNA_half <- (1.0 - b_2) * exp(-a_2 * x) + b_2
                        (mRNA_half - 0.5)^2
                    }
                    
                    out2 <- optim(1,model2_half)
                    half2 <- out2$par
                    
                    if(b_2 >= 0.5){
                        half2 <- Inf
                    }
                    
                    mRNA_pred <- (1.0 - b_2) * exp(-a_2 * time_point_exp$hour) + b_2
                    s2 <- sum((time_point_exp$exp - mRNA_pred)^2) / data_point
                    AIC2 <- data_point * log(s2) + (2 * 1)
                    
                    cat(min2,cor2,a_2,b_2,half2,AIC2,sep="\t",file=output_file,append=T)
                    cat("\t",file=output_file,append=T)
                    ############################################################
                    
                    ###Model3: f(t) = c * exp(-a * t) + (1 - c) * exp(-b * t)###
                    optim3 <- function(x){
                        mRNA_exp <- x[3] * exp(-x[1] * time_point_exp$hour) + (1.0 - x[3]) * exp(-x[2] * time_point_exp$hour)
                        sum((time_point_exp$exp - mRNA_exp)^2)
                    }
                    
                    out3 <- optim(c(1,1,0.1),optim3)
                    min3 <- out3$value
                    a_3 <- out3$par[1]
                    b_3 <- out3$par[2]
                    c_3 <- out3$par[3]
                    
                    model3_pred <- function(a,b,c){
                        mRNA_exp <- c * exp(-a * time_point_exp$hour) + (1.0 - c) * exp(- b * time_point_exp$hour)
                        (cor(mRNA_exp, time_point_exp$exp, method="pearson"))^2
                    }
                    
                    cor3 <- model3_pred(a_3,b_3,c_3)
                    
                    model3_half <- function(x){
                        mRNA_half <- c_3 * exp(-a_3 * x) + (1.0 - c_3) * exp(-b_3 * x)
                        (mRNA_half - 0.5)^2
                    }
                    
                    out3 <- optim(1,model3_half)
                    half3 <- out3$par
                    mRNA_pred <- c_3 * exp(-a_3 * time_point_exp$hour) + (1.0 - c_3) * exp(-b_3 * time_point_exp$hour)
                    s2 <- sum((time_point_exp$exp - mRNA_pred)^2) / data_point
                    AIC3 <- data_point * log(s2) + (2 * 2)
                    
                    cat(min3,cor3,a_3,b_3,c_3,half3,AIC3,sep="\t",file=output_file,append=T)
                    cat("\t",file=output_file,append=T)
                    ############################################################
                    
                    half_table <- data.frame(half=c(as.numeric(half1_2),as.numeric(half2),as.numeric(half3)),AIC=c(as.numeric(AIC1),as.numeric(AIC2),as.numeric(AIC3)))
                    AIC_list <- order(half_table$AIC)
                    #selected_half <- half_table$half[min_AIC_index]
                    
                    AIC_flg <- 0
                    for(min_AIC_index in AIC_list){
                        ###Model1###
                        if(min_AIC_index == 1){
                            selected_half <- half1_2
                            if(selected_half > 24){
                                selected_half <- 24
                            }
                            if(a_1 > 0){
                                cat("model1",cor1,selected_half,sep="\t",file=output_file,append=T)
                                AIC_flg <- 1
                                break
                            }
                        }
                        ###Model2###
                        if(min_AIC_index == 2){
                            selected_half <- half2
                            if(selected_half == "Inf"){
                                selected_half <- 24
                            }else if(selected_half > 24){
                                selected_half <- 24
                            }
                            if(a_2 > 0 && b_2 > 0 && b_2 < 1){
                                cat("model2",cor2,selected_half,sep="\t",file=output_file,append=T)
                                AIC_flg <- 1
                                break
                            }
                        }
                        ###Model3###
                        if(min_AIC_index == 3){
                            selected_half <- half3
                            if(selected_half > 24){
                                selected_half <- 24
                            }
                            if(a_3 > 0 && b_3 > 0 && c_3 > 0 && c_3 < 1){
                                cat("model3",cor3,selected_half,sep="\t",file=output_file,append=T)
                                AIC_flg <- 1
                                break
                            }
                        }
                    }
                    
                    if(AIC_flg == 0){
                        if(a_1 < 0){
                            cat("model1",cor1,24,sep="\t",file=output_file,append=T)
                            next
                        }else{
                            cat("no_good_model","NA",24,sep="\t",file=output_file,append=T)
                            next
                        }
                    }

                }else{
                    cat("few_data","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA",
                        "NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA", sep="\t", file=output_file, append=T)
                }
            }else{
                cat("low_expresion","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA",
                    "NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA", sep="\t", file=output_file, append=T)
            }
        }
        cat("\n", file=output_file, append=T)
    }
}

