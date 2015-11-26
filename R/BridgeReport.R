#Title: BridgeReport: Data visualization for RNA decay curve with shiny library
#Auther: Naoto Imamachi
#ver: 1.0.0
#Date: 2015-11-18

BridgeReport <- function(filename1 = "siStealth_compatible_genes_RefSeq_result_mRNA.fpkm_table",
                         filename2 = "siPUM1_compatible_genes_RefSeq_result_mRNA.fpkm_table",
                         filename3 = "BridgeR_5C_HalfLife_calculation_R2_selection.txt",
                         group,
                         hour,
                         ComparisonFile,
                         SearchRow = "symbol",
                         InforColumn = 4,
                         Color = c("black","red"),
                         CutoffDataPoint1 = c(1,2),
                         CutoffDataPoint2 = c(8,12))
{
    time_points <- length(hour) 
    group_number <- length(group)
    rpkm_file1 <- fread(filename1, header=T)
    rpkm_file2 <- fread(filename2, header=T)
    input_file <- fread(filename3, header=T)
    
    table_header <- NULL
    for(x in 1:time_points){
        table_header <- append(table_header,paste(hour[x],"hr",sep=""))
    }
    
    ###Prepare_file_infor###
    comp_file_number <- NULL
    for(a in 1:length(ComparisonFile)){
        comp_file_number <- append(comp_file_number, which(group == ComparisonFile[a]))
    }
    setkeyv(input_file,SearchRow)
    setkeyv(rpkm_file1,SearchRow)
    setkeyv(rpkm_file2,SearchRow)
    
    ###Prepare_except_time_point###
    CutoffDataPoint1 <- sort(CutoffDataPoint1, decreasing = T)
    CutoffDataPoint2 <- sort(CutoffDataPoint2, decreasing = F)

    delete_times_list <- function(test_times){
        times_length <- length(test_times)
        times_index <- c(times_length)
        add_index <- times_length
        label_list <- NULL
        
        for(counter in times_length:1){
            if(length(times_index) != 1){
                check_times <- test_times[times_index]
                del_label <- paste("Delete ",paste(check_times,collapse = "-"),"hr",sep="")
                label_list <- append(label_list, del_label)
            }

            #Counter
            add_index <- add_index - 1
            times_index <- append(times_index, add_index)
        }
        return(label_list)
    }
    
    label_list0 <- c("Raw")
    label_list1 <- delete_times_list(CutoffDataPoint1)
    label_list2 <- delete_times_list(CutoffDataPoint2)
    label_list3 <- NULL
    
    label_name <- paste("Delete ",hour[-1],"hr",sep="")
    label_list3 <- append(label_list3, label_name)

    label_list_all <- c(label_list0, label_list1, label_list2, label_list3)
    select_list <- NULL
    for(x in 1:length(label_list_all)){
        select_list <- append(select_list,list(x))
    }
    names(select_list) <- label_list_all

    ###Function1###
    select_exp_data <- function(time_point_exp, input_select_id){
        time_point_exp_for_test <- time_point_exp
        select_name <- names(select_list[as.numeric(input_select_id)]) #Raw, Delete 1h-2h...
        if(select_name == "Raw"){
        }else{
            select_name <- gsub("Delete ", "", select_name)
            select_name <- gsub("hr", "", select_name)
            select_name <- as.numeric(as.vector(strsplit(select_name, "-")[[1]]))
            for(i in select_name){
                time_point_exp_for_test <- time_point_exp_for_test[time_point_exp_for_test$hour != i,]
            }
        }
        return(time_point_exp_for_test)
    }
    
    rpkm_exp_st <- InforColumn + 1
    rpkm_exp_ed <- InforColumn + time_points
    
    ui <- fluidPage(
        titlePanel("BridgeReport ver 0.2.0"),
        
        sidebarLayout(
            position="left",
            sidebarPanel(
                helpText("Select a gene sysbol to examine. Information will be collected from your BRIC-seq dataset."),
                textInput("text",
                          label = "Input gene symbol",
                          value = "Enter text..."),
                sliderInput("range_x", 
                            label = "X-axis(Time course):",
                            min = 0, max = 12, value = c(0, 12)),
                sliderInput("range_y", 
                            label = "Y-axis(Relative RNA remaining):",
                            min = 0.001, max = 10, value = c(0.1,1.5)),
                selectInput("select1", label = paste(ComparisonFile[1]," (Select time points)",sep=""), 
                            choices = select_list, selected = 1),
                selectInput("select2", label = paste(ComparisonFile[2]," (Select time points)",sep=""), 
                            choices = select_list, selected = 1)
            ),
            
            mainPanel(
                #h1("test"),
                #textOutput("text1"),
                
                fluidRow(
                    column(width = 5,
                           plotOutput("plot1",
                                      click = "plot1_click",
                                      brush = brushOpts(id = "plot1_brush")
                           )
                    )
                ),
                
                tableOutput("mytable1")
                
            )
        )
    )
    
    server <- function(input, output) {
        r_squared_list <- NULL
        half_life_list <- NULL
        model_list <- NULL
        
        output$plot1 <- renderPlot({
            data <- as.vector(as.matrix(input_file[input$text]))
            gene_name <- as.character(input$text)
            
            ###Prepare_ggplot2###
            p <- ggplot()
            
            flg <- 0
            fig_color <- NULL
            for(a in comp_file_number){
                if(flg == 0){
                    fig_color <- Color[1]
                }else if(flg == 1){
                    fig_color <- Color[2]
                }
                infor_st <- 1 + (a - 1)*(time_points + InforColumn + 3)
                infor_ed <- (InforColumn)*a + (a - 1)*(time_points + 3)
                exp_st <- infor_ed + 1
                exp_ed <- infor_ed + time_points
                model_index <- infor_ed + time_points + 1
                
                exp <- as.numeric(data[exp_st:exp_ed])
                model_name <- as.character(data[model_index])
                model_list <- append(model_list, model_name)
                time_point_exp_original <- data.frame(hour,exp)
                
                # For storing which rows have been excluded
                #vals <- reactiveValues(
                #    keeprows = rep(TRUE, nrow(time_point_exp_original))
                #)
                
                time_point_exp <- time_point_exp_original[time_point_exp_original$exp > 0, ]
                
                if(flg == 0){
                    time_point_exp_for_test <- select_exp_data(time_point_exp, input$select1)
                }else if(flg == 1){
                    time_point_exp_for_test <- select_exp_data(time_point_exp, input$select2)
                }

                p <- p + layer(data=time_point_exp_for_test, 
                               mapping=aes(x=hour, y=exp), 
                               geom="point",
                               size=4,
                               shape=19,
                               colour=fig_color)
                
                data_point <- length(time_point_exp_for_test$exp)
                if(!is.null(time_point_exp_for_test)){
                    if(data_point >= 3){
                        if(as.numeric(as.vector(as.matrix(time_point_exp_for_test$exp[1]))) > 0){
                            model <- lm(log(time_point_exp_for_test$exp) ~ time_point_exp_for_test$hour - 1)
                            model_summary <- summary(model)
                            coef <- -model_summary$coefficients[1]
                            half_life <- log(2)/coef
                            if(coef < 0 || half_life >= 24){
                                half_life <- 24
                            }
                            r_squared <- model_summary$r.squared
                            half_life <- round(half_life,digits=3)
                            r_squared <- round(r_squared,digits=3)
                            half_life_list <- append(half_life_list,half_life)
                            r_squared_list <- append(r_squared_list,r_squared)
                            
                            fig_data <- data.frame(hour=time_point_exp_for_test$hour)
                            predicted <- as.numeric(as.vector(as.matrix(predict(model, fig_data))))
                            fig_data$exp <- exp(as.vector(as.matrix(predicted)))
                            p <- p + layer(data=fig_data,
                                           mapping=(aes(x=hour, y=exp)),
                                           geom="line",
                                           size=1.2,
                                           colour=fig_color)
                        }else{
                            half_life_list <- append(half_life_list,"NA")
                            r_squared_list <- append(r_squared_list,"NA")
                        }
                    }else{
                        half_life_list <- append(half_life_list,"NA")
                        r_squared_list <- append(r_squared_list,"NA")
                    }
                }else{
                    half_life_list <- append(half_life_list,"NA")
                    r_squared_list <- append(r_squared_list,"NA")
                }
                flg <- 1
            }
            p <- p + ggtitle(gene_name)
            p <- p + xlab("Time")
            p <- p + ylab("Relative RPKM (Time0 = 1)")
            p <- p + xlim(input$range_x[1],input$range_x[2])
            ybreaks1 <- seq(0,10,0.1)[2:101]
            ybreaks2 <- seq(0,0.1,0.01)[2:10]
            ybreaks <- c(ybreaks1,ybreaks2)
            p <- p + scale_y_log10(breaks=ybreaks, labels=ybreaks, limits=c(input$range_y[1],input$range_y[2]))
            
            if(is.null(r_squared_list) || is.null(half_life_list)){
                
                table_data <- data.frame(R2=c("NA","NA"), HalfLife=c("NA","NA"))
            }else{
                
                rpkm_data1 <- as.character(round(as.numeric(as.vector(as.matrix(rpkm_file1[input$text]))[rpkm_exp_st:rpkm_exp_ed]),digits=3))
                rpkm_data2 <- as.character(round(as.numeric(as.vector(as.matrix(rpkm_file2[input$text]))[rpkm_exp_st:rpkm_exp_ed]),digits=3))
                table_data <- data.frame(R2=as.character(r_squared_list), HalfLife=as.character(half_life_list), model=model_list)
                table_data <- cbind(table_data,rbind(rpkm_data1,rpkm_data2))
                colnames(table_data) <- c("R2","Half-life","Model",table_header)
                rownames(table_data) <- as.character(ComparisonFile)
            }
            
            output$mytable1 = renderTable({
                table_data
            })
            
            p
        })
    }
    
    app <- shinyApp(ui, server) # shinyApp() is also OK!
    runApp(app)
}
