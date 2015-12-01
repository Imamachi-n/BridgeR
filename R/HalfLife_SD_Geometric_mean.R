#Title: Halflife/RPKM_SD_Geometric_Means
#Auther: Naoto Imamachi
#ver: 1.0.0
#Date: 2015-12-01

BridgeRHalfSD <- function(HalflifeFiles = c("BridgeR_6_HalfLife_Pvalue_estimation_PUM1_study.txt",
                                            "BridgeR_6_HalfLife_Pvalue_estimation_PUM2_study.txt",
                                            "BridgeR_6_HalfLife_Pvalue_estimation_PUM1_2_study.txt"),
                          RPKMFiles = c("siStealth_PUM1_study_compatible_genes_RefSeq_result_mRNA.fpkm_table",
                                        "siCTRL_PUM2_study_compatible_genes_RefSeq_result_mRNA.fpkm_table",
                                        "siCTRL_PUM1_2_study_RefSeq_compatible_genes_RefSeq_result_mRNA.fpkm_table"),
                          hour = c(0,1,2,4,8,12),
                          InforColumn=4,
                          CompFiles = c(1,1,1),
                          CompNames = c("siStealth_PUM1_study","siCTRL_PUM2_study","siCTRL_PUM1_2_study"),
                          OutputFile = "HalfLife_RPKM_mean_SD.txt"
                          ){
    ###Prepare_file_infor###
    TimeDataSize <- length(hour)
    
    Half_index_list <- NULL
    for(SampleIndex in CompFiles){
        infor_st <- 1 + (SampleIndex - 1)*(TimeDataSize + InforColumn + 13) #13 => Infor: Model, R2, half_life, etc...
        infor_ed <- (InforColumn)*SampleIndex + (SampleIndex - 1)*(TimeDataSize + 13) #13 => Infor: Model, R2, half_life, etc...
        Half_index <- infor_ed + TimeDataSize + 8
        Half_index_list <- append(Half_index_list, Half_index)
    }
    
    RPKM_index_list <- NULL
    for(x in 1:length(RPKMFiles)){
        infor_st <- 1
        infor_ed <- InforColumn
        RPKM_index <- infor_ed + 1
        RPKM_index_list <- append(RPKM_index_list, RPKM_index)
    }
    
    ###Infor_data###
    infor_data <- fread(RPKMFiles[1], header=T)[,1:InforColumn,with=F]
    infor_header <- as.vector(as.matrix((fread(RPKMFiles[1], nrows=1, header=F))))[1:4]
    
    ###Read_halflife_data###
    Half_data_table <- NULL
    for(x in 1:length(HalflifeFiles)){
        input_file <- fread(HalflifeFiles[x], header=T)[[Half_index_list[x]]]
        Half_data_table <- cbind(Half_data_table, input_file)
    }
    
    RPKM_data_table <- NULL
    for(x in 1:length(RPKMFiles)){
        input_file <- fread(RPKMFiles[x], header=T)[[RPKM_index_list[x]]]
        RPKM_data_table <- cbind(RPKM_data_table, input_file)
    }
    
    Half_data_table <- as.data.frame(Half_data_table)
    RPKM_data_table <- as.data.frame(RPKM_data_table)
    
    ###Half_variance_distribution###
    Half_variance <- NULL
    for(x in 1:length(Half_data_table[,1])){
        half_data <- as.numeric(as.vector(as.matrix(Half_data_table[x,])))
        half_data <- half_data[!is.na(half_data)]
        if(length(half_data) == 0){
        }else{
            half_data <- half_data - mean(half_data)
        }
        Half_variance <- append(Half_variance, half_data)
    }
    
    Half_variance_table <- as.data.frame(Half_variance)
    ###Fig_data: RPKM(mean) vs Half-life(SD)###
    png(filename="Histgram_HalfLife_variance.png",width = 600, height = 600)
    
    p <- ggplot()
    p <- p + layer(data = Half_variance_table,
                   mapping = aes(x=Half_variance),
                   geom = "bar",
                   geom_params = list(fill = "steelblue", color = "white"),
                   stat = "bin",
                   stat_params = list(binwidth = 0.1))
    #p <- p + xlab("RPKM mean") + ylab("Half-life SD")
    p <- p + xlim(-5,5)
    plot(p)
    
    dev.off() #close_fig
    plot.new()
    
    
    ####Calc_mean & SD function###
    test_func <- function(table){
        list(means=mean(table,na.rm = T), sds=sd(table,na.rm = T))
    }
    
    halflife_data <- as.data.frame(t(Half_data_table))
    halflife_sd <- as.data.frame(t(sapply(halflife_data, test_func)))
    
    RPKM_data <- as.data.frame(t(RPKM_data_table))
    RPKM_sd <- as.data.frame(t(sapply(RPKM_data, test_func)))
    
    fig_data <- cbind(halflife_sd,RPKM_sd)
    colnames(fig_data) <- c("Half_mean","Half_SD","RPKM_mean","RPKM_SD")
    
    ###fig_color###
    color_vec <- NULL
    true_data <- 0
    for(x in 1:length(fig_data[,1])){
        if(is.na(fig_data[x,2])){
            color_vec <- append(color_vec, "black")
            next
        }
        if(as.numeric(fig_data[x,2]) > 1.5){
            color_vec <- append(color_vec, "red")
        }else{
            color_vec <- append(color_vec, "black")
            true_data <- true_data + 1
        }
    }
    
    false_data <- length(color_vec[color_vec == "red"])
    
    print(paste("True data: ", as.character(true_data), sep=""))
    print(paste("False data: ", as.character(false_data), sep=""))
    
    result_data <- cbind(infor_data,Half_data_table,RPKM_data_table,fig_data)
    half_names <- paste("HalfLife_",CompNames,sep="")
    RPKM_names <- paste("RPKM_",CompNames,sep="")
    colnames(result_data) <- c(infor_header,half_names,RPKM_names,c("HalfLife_mean","HalfLife_SD","RPKM_mean","RPKM_SD"))
    
    result_table_pre <- data.frame(lapply(result_data, as.character))
    result_table <- write.table(result_table_pre, file=OutputFile, sep="\t", row.names = F)
    
    ###Fig_data: RPKM(mean) vs Half-life(SD)###
    png(filename="RPKM_mean_vs_HalfLife_SD.png",width = 600, height = 600)
    
    p <- ggplot()
    p <- p + layer(data = fig_data,
                   mapping = aes(x=as.numeric(RPKM_mean), y=as.numeric(Half_SD)),
                   geom = "point",
                   alpha=0.4)
    p <- p + xlab("RPKM mean") + ylab("Half-life SD")
    p <- p + xlim(0,100) + ylim(0,10)
    #p <- p + scale_colour_manual(values=c("black","red")) + theme(legend.position="none")
    plot(p)
    
    dev.off() #close_fig
    plot.new()
    
    ###Fig_data: Half-life(SD) vs Half-life(mean)###
    png(filename="HalfLife_mean_vs_HalfLife_SD.png",width = 600, height = 600)
    
    p <- ggplot()
    p <- p + layer(data = fig_data,
                   mapping = aes(x=as.numeric(Half_mean), y=as.numeric(Half_SD)), #, colour=factor(color_vec)),
                   geom = "point",
                   alpha=0.4)
    p <- p + xlab("Half-life mean") + ylab("Half-life SD")
    p <- p + xlim(0,24) + ylim(0,15)
    #p <- p + scale_colour_manual(values=c("black","red")) + theme(legend.position="none")
    plot(p)
    
    dev.off()
    plot.new()
}

###TEST###
#setwd("C:/Users/Naoto/Documents/github/test_R_libraries/BridgeR_test/2015-11-30_PUM1_study_CTRL_comparison/test")
#BridgeRHalfSD()
