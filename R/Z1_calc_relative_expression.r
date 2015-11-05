#Title: Calc_relative_expression
#Auther: Naoto Imamachi
#ver: 1.0.0
#Date: 2015-10-07

###Calc_relative_expression_function###
BridgeRDataSetFromCuffnorm <- function(CuffnormFiles, group, hour, cutoff = 0.1, InforColumn = 4, OutputFile = "BridgeR_1_Relative_expression_data.txt"){
    ###Prepare_files###
    time_points <- length(hour)
    input_file_numbers <- length(CuffnormFiles)
    input_file <- NULL
    for(filename in CuffnormFiles){
        if(is.null(input_file)){
            input_file <- suppressWarnings(fread(filename, header=T))
        }else{
            input_file <- cbind(input_file,suppressWarnings(fread(filename, header=T)))
        }
    }
    output_file <- OutputFile
    
    ###print_header###
    cat("",file=output_file)
    hour_label <- NULL
    for(a in 1:length(group)){
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

    ###Read_each_line###
    gene_number <- length(input_file[[1]]) #Total number of genes
    sample_size <- length(group)
    for (x in 1:gene_number){
        data <- as.vector(as.matrix(input_file[x,]))
        for(a in 1:sample_size){
            if(a != 1){
                cat("\t", sep="", file=output_file, append=T)
            }
            ###Infor_data###
            infor_st <- 1 + (a - 1)*(time_points + InforColumn)
            infor_ed <- (InforColumn)*a + (a - 1)*time_points
            gene_infor <- data[infor_st:infor_ed]
            cat(gene_infor, sep="\t", file=output_file, append=T)
            cat("\t", sep="", file=output_file, append=T)
            
            ###Exp_data###
            exp_st <- infor_ed + 1
            exp_ed <- infor_ed + time_points
            exp <- data[exp_st:exp_ed]
            exp <- as.numeric(exp)
            start_time <- exp[1]
            ###if: Time0 <= cutoff_RPKM###
            if(start_time <= cutoff){
                cat(rep(0,time_points), sep="\t", file=output_file, append=T)
                next
            }
            rel_exp <- NULL
            for (y in 1:time_points){
                rel_exp <- append(rel_exp, exp[y]/start_time)
            }
            cat(rel_exp, sep="\t", file=output_file, append=T)
        }
        cat("\n", sep="", file=output_file, append=T)
    }
}
