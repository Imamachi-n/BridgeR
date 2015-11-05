#Title: Normalization
#Auther: Naoto Imamachi
#ver: 1.0.0
#Date: 2015-10-08

###Estimate_normalization_factor_function###
BridgeRNormalizationForLuc2 <- function(filename, group, hour, InforColumn = 4, NormFactor = "BridgeR_3_Normalizaion_factor_dataset", OutputFile = "BridgeR_4_Normalized_expression_data.txt"){
    ###Load_NormFactor_file###
    group_number <- length(group)
    time_points <- length(hour)
    nf_st <- 2
    nf_ed <- time_points + 1
    nf_Luc2 <- NULL

    normalization_factor <- fread(NormFactor, header=T)[,nf_st:nf_ed,with=F]
    for(a in 1:group_number){
        if(is.null(nf_Luc2)){
            nf_Luc2 <- as.vector(as.matrix(normalization_factor[a,]))
        }else{
            nf_Luc2 <- append(nf_Luc2, as.vector(as.matrix(normalization_factor[a,])))
        }
    }
    
    ###Load_files###
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
            label <- x
            if(x < 10){
                label <- paste("0",x,sep="")
            }
            hour_label <- append(hour_label, paste("T", label, "_", a, sep=""))
        }
        
        infor_st <- 1 + (a - 1)*(time_points + InforColumn)
        infor_ed <- (InforColumn)*a + (a - 1)*time_points
        infor <- colnames(input_file)[infor_st:infor_ed]
        cat(infor,hour_label, sep="\t", file=output_file, append=T)
    }
    cat("\n", sep="", file=output_file, append=T)
    
    ###calc_normalized_Relative_RPKM###
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
            exp <- as.numeric(data[exp_st:exp_ed])
            
            nf_Luc2_st <- (time_points)*(a - 1) + 1
            nf_Luc2_ed <- (time_points)*a
            nf_Luc2_exp <- nf_Luc2[nf_Luc2_st:nf_Luc2_ed]

            normalized_exp <- NULL
            normalized_exp <- exp/nf_Luc2_exp
            
            cat(gene_infor, sep="\t", file=output_file, append=T)
            cat("\t", sep="\t", file=output_file, append=T)
            cat(normalized_exp, sep="\t", file=output_file, append=T)
        }
        cat("\n", sep="\t", file=output_file, append=T)
    }
}
