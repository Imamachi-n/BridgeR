#Title: BridgeReport: Data visualization for RNA decay curve with shiny library
#Auther: Naoto Imamachi
#ver: 1.0.0
#Date: 2015-11-18

BridgeReport <- function(filename1="siStealth_compatible_genes_RefSeq_result_mRNA.fpkm_table",
                         filename2="siPUM1_compatible_genes_RefSeq_result_mRNA.fpkm_table",
                         filename3="BridgeR_4_Normalized_expression_dataset.txt",
                         group,
                         hour,
                         ComparisonFile,
                         SearchRow="symbol",
                         InforColumn=4)
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
    
    rpkm_exp_st <- InforColumn + 1
    rpkm_exp_ed <- InforColumn + time_points
    
    ui <- fluidPage(
        titlePanel("BridgeReport ver 0.1.0"),
        
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
                selectInput("select", label = "Select decay curve calculation method", 
                            choices = list("Raw" = 1,
                                           "Except 1hr" = 2,
                                           "Except 1-2hr" = 3,
                                           "Except 2hr" = 4,
                                           "Except 4hr" = 5,
                                           "Except 8hr" = 6,
                                           "Except 12hr" = 7,
                                           "Except 8-12hr" = 8,
                                           "Except 1,12hr" = 9,
                                           "Cutoff>=0.1" = 10,
                                           "Cutoff>=0.05" = 11,
                                           "Cutoff>=0.01" = 12), selected = 1)
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
        
        output$plot1 <- renderPlot({
            data <- as.vector(as.matrix(input_file[input$text]))
            gene_name <- as.character(data[2])
            
            ###Prepare_ggplot2###
            p <- ggplot()
            
            flg <- 0
            fig_color <- NULL
            for(a in comp_file_number){
                if(flg == 0){
                    fig_color <- "black"
                }else{
                    fig_color <- "red"
                }
                infor_st <- 1 + (a - 1)*(time_points + InforColumn)
                infor_ed <- (InforColumn)*a + (a - 1)*(time_points)
                exp_st <- infor_ed + 1
                exp_ed <- infor_ed + time_points
                
                exp <- as.numeric(data[exp_st:exp_ed])
                time_point_exp_original <- data.frame(hour,exp)
                
                # For storing which rows have been excluded
                #vals <- reactiveValues(
                #    keeprows = rep(TRUE, nrow(time_point_exp_original))
                #)
                
                time_point_exp <- time_point_exp_original[time_point_exp_original$exp > 0, ]
                if(input$select == 1){
                    
                }else if(input$select == 2){
                    time_point_exp <- time_point_exp[time_point_exp$hour!=1,]
                }else if(input$select == 3){
                    time_point_exp <- time_point_exp[time_point_exp$hour!=1,]
                    time_point_exp <- time_point_exp[time_point_exp$hour!=2,]
                }else if(input$select == 4){
                    time_point_exp <- time_point_exp[time_point_exp$hour!=2,]
                }else if(input$select == 5){
                    time_point_exp <- time_point_exp[time_point_exp$hour!=4,]
                }else if(input$select == 6){
                    time_point_exp <- time_point_exp[time_point_exp$hour!=8,]
                }else if(input$select == 7){
                    time_point_exp <- time_point_exp[time_point_exp$hour!=12,]
                }else if(input$select == 8){
                    time_point_exp <- time_point_exp[time_point_exp$hour!=12,]
                    time_point_exp <- time_point_exp[time_point_exp$hour!=8,]
                }else if(input$select == 9){
                    time_point_exp <- time_point_exp[time_point_exp$hour!=1,]
                    time_point_exp <- time_point_exp[time_point_exp$hour!=12,]
                }else if(input$select == 10){
                    time_point_exp <- time_point_exp[time_point_exp$exp >= 0.1, ]
                }else if(input$select == 11){
                    time_point_exp <- time_point_exp[time_point_exp$exp >= 0.05, ]
                }else if(input$select == 12){
                    time_point_exp <- time_point_exp[time_point_exp$exp >= 0.01, ]
                }
                
                p <- p + layer(data=time_point_exp, 
                               mapping=aes(x=hour, y=exp), 
                               geom="point",
                               size=4,
                               shape=19,
                               colour=fig_color)
                
                data_point <- length(time_point_exp$exp)
                if(!is.null(time_point_exp)){
                    if(data_point >= 3){
                        if(as.numeric(as.vector(as.matrix(time_point_exp$exp[1]))) > 0){
                            model <- lm(log(time_point_exp$exp) ~ time_point_exp$hour - 1)
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
                            
                            fig_data <- data.frame(hour=time_point_exp$hour)
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
                
                rpkm_data1 <- round(as.numeric(as.vector(as.matrix(rpkm_file1[input$text]))[rpkm_exp_st:rpkm_exp_ed]),digits=3)
                rpkm_data2 <- round(as.numeric(as.vector(as.matrix(rpkm_file2[input$text]))[rpkm_exp_st:rpkm_exp_ed]),digits=3)
                table_data <- data.frame(R2=r_squared_list, HalfLife=half_life_list)
                table_data <- cbind(table_data,rbind(rpkm_data1,rpkm_data2))
                colnames(table_data) <- c("R2","Half-life",table_header)
                rownames(table_data) <- ComparisonFile
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
