# Installing packages, load packages

if(!require("dplyr")) install.packages("dplyr")
if(!require("stringr")) install.packages("stringr")
if(!require("ggplot2")) install.packages("ggplot2")
if(!require("ggpubr")) install.packages("ggpubr")
if(!require("shiny")) install.packages("shiny")
if(!require("shinyFiles")) install.packages("shinyFiles")
if(!require("egg")) install.packages("egg")
if(!require("shinythemes")) install.packages("shinythemes")
library(shinythemes)
library(egg)
library(shiny)
library(shinyFiles)
library(dplyr)
library(stringr)
library(ggplot2)
library(ggpubr)

# Starts a web page, choose a theme
ui <- fluidPage(theme = shinytheme("journal"),
                # Navigation bar name
                navbarPage("SNPulse",
                           # First tab
                           tabPanel("Home",
                                    sidebarPanel(
                                      fileInput(inputId = "bim",
                                                label = "Upload the .bim file", 
                                                multiple = FALSE,
                                                accept = ".bim"),
                                      fileInput(inputId = "ped"
                                                ,label = "Upload the .ped file", 
                                                multiple = FALSE,
                                                accept = ".ped"),
                                      fileInput(inputId = "fam",
                                                label = "Upload the .fam file", 
                                                multiple = FALSE,
                                                accept = ".fam"),
                                      fileInput(inputId = "txt",
                                                label = "Upload the .txt file (contains the ID and mean date)", 
                                                multiple = FALSE,
                                                accept = c("t")),
                                      numericInput(inputId = "step", 
                                                   label = "Choose time step", 
                                                   value = 2500),
                                      textInput(inputId = "snp", label = "Name of SNP", value = ""),
                                      actionButton("submit", "Submit"),
                                    ), #sidebarPanel
                                    
                                    mainPanel(
                                      # Show the output (plot)
                                      plotOutput(outputId = "output", width = "80%", height = 660),
                                      # Add a download button
                                      downloadButton("download")
                                    ) # mainPanel
                           ), # tabPanel #Home
                           
                           tabPanel("Readme", 
                                    titlePanel("Readme"), 
                                    div(includeMarkdown("readme.md"), 
                                        align="justify")
                           ) #tabPanel #About
                ) # navabarPage
) #fluidPage

server <- function(input, output, session) {
  options(shiny.maxRequestSize=2000*1024^2)
  
  # Function 0
  make_ped <- function(snp,ped,bim){
    # Read bim and ped file
    bim <- read.table(bim$datapath, header=FALSE)
    ped <- read.delim(ped$datapath, header=FALSE) 
    
    # Extract snp name and snp list
    snp_name = snp
    snp_list = bim$V2
    
    # Use snp name and snp list to get the col_index_snp
    col_index_snp=match(snp_name,snp_list)+6
    
    # Extract the first 6 cols and the allele count of the snp
    snp_ped = select(ped, all_of(c(1,2,col_index_snp)))
    
    return(snp_ped)
  }
  
  # Function 1
  make_primary_data <- function(txt_file, fam_file) {
    # Read txt file, this file has ID and mean date
    date_data <- read.table(txt_file$datapath, header = FALSE) 
    # Make sure that is a numeric data
    date_data$V2 <- as.numeric(date_data$V2) 
    
    # Read fam file
    fam_data <- read.table(fam_file$datapath, header = FALSE) 
    filter_data <- unique(date_data[date_data$V1 %in% fam_data$V2,])
    return(filter_data)
  }
  
  # Function 2
  create_time_vec <- function(step,data){
    # Make step always positive
    step = abs(step)
    # Create end_time
    round_digit = nchar(max(data$V2))-2
    end_time = round(max(data$V2),-round_digit)
    # Start time
    start_time = 0
    # Make time vector
    time_vec = seq(start_time, end_time+step, step)
    return(time_vec)
  }
  
  # Function 3
  make_secondary_data <- function(time_vec,ped,primary){
    #Making return data
    ped_data_list=list()
    
    # Add data to the return data
    for (i in 1:(length(time_vec)-1)) { 
      # Create time data, add to a time_data list
      time_temp = primary[primary$V2 >= time_vec[i] & primary$V2 < time_vec[i+1], ]
      
      # Separate ped data, 
      ped_temp = ped[ped$V2 %in% time_temp$V1,]
      ped_data_list[[i]]=ped_temp
    }
    return(ped_data_list)
  }
  
  # Function 4
  allele_data <- function(snp, bim, ped, time_vec){
    bim <- read.table(bim$datapath, header=FALSE)
    
    snp_name = snp
    snp_list=bim$V2
    
    # Find column index
    col_index_allele=match(snp_name,snp_list)
    
    # Retrieve the minor and major
    minor_temp = as.character(bim$V5[col_index_allele])
    major_temp = as.character(bim$V6[col_index_allele])
    
    #Create vector to save the data
    minor_count_vec = NULL
    major_count_vec = NULL
    maf_vec = NULL
    
    for (i in 1:length(ped)){
      # Retrieve ped data
      sep_ped_temp=ped[[i]][,3]
      
      # Look for minor and major allele in bim file
      minor_temp = as.character(bim$V5[col_index_allele])
      major_temp = as.character(bim$V6[col_index_allele])
      
      # Paste the allele data to a string
      allele_string = paste(sep_ped_temp, collapse=' ')
      
      # Count minor, major and calculate maf
      # Minor count and append the the current count to a vector
      minor_count = str_count(allele_string, minor_temp)
      minor_count_vec = c(minor_count_vec,minor_count)
      
      # Minor count and append the the current count to a vector
      major_count = str_count(allele_string, major_temp)
      major_count_vec = c(major_count_vec, major_count)
      
      # Calculate maf and append the maf to a vector
      if (minor_count+major_count==0){
        maf = 0
      } else {
        maf = minor_count/(minor_count+major_count)
      }
      maf_vec = c(maf_vec, maf) 
    }
    
    #Label minor and major allele
    base_lab <- c("A","C","G","T")
    # Label minor
    result_minor <- tryCatch({
      as.numeric(minor_temp)
    }, error = function(e) {
      return(NULL)
    })
    if (is.null(result_minor)){
      minor_lab=base_temp
    } else (minor_lab=base_lab[result_minor])
    
    #Label major
    result_major <- tryCatch({
      as.numeric(major_temp)
    }, error = function(e) {
      return(NULL)
    })
    if (is.null(result_major)){
      major_lab = base_temp
    } else (major_lab = base_lab[result_major])
    
    # Create new time_vec, remove the year 0
    time_vec_temp <- time_vec[-1]
    
    # Return a data frame
    final_data_frame <- data.frame(minor_count_vec, major_count_vec,
                                   time_vec_temp, maf_vec)
    colnames(final_data_frame) <- c(minor_lab,major_lab,"Year","MAF")
    
    #Return the allele data
    return(c(final_data_frame,snp_name))
  }
  
  # Function 5
  plotting_maf <- function(data,time_vec){
    # Extract the minor, major, snp for label later 
    snp_name = data[5]
    time_vec = time_vec
    
    minor_lab=colnames(as.data.frame(data[1:3][1]))
    major_lab=colnames(as.data.frame(data[1:3][2]))
    label = paste(snp_name," - ",minor_lab,"/",major_lab,sep="")
    
    # Make data for bar plot
    plot_data_bar <- as.data.frame(data[1:3])
    colnames(plot_data_bar) <- c("Minor","Major","Year")
    
    # Transverse the plot_data_bar for stacked bar plotting
    plot_data_bar_long <- reshape2::melt(plot_data_bar, 
                                         id.vars= "Year", 
                                         variable.name = "Allele", 
                                         value.name = "Count")
    
    # Make data for line plot
    plot_data_line <- as.data.frame(data[3:4])
    
    # Find time_vec resolution, to make the plot beautiful
    # Part of the trickiness here is that padding works differently 
    # for points/lines vs. bars, which have width around the data point. 
    # So if we want our points to be aligned with our bars, we need to tweak that.
    # To estimate what this should be, I note that, for example the time step 
    # resolution is 2500, the bars will be 90% of that width (2,250), so each side
    # of the bar extends 1,125 beyond the range that the points need. 
    # So I add 1,125 to the point padding, to this line expand = expansion(add = res)
    
    step = max(time_vec)/(length(time_vec)-1)
    res <- rep((0.9*step)/2+500,2)
    
    # Plot line+point plot
    plot1 <- ggplot() +
      geom_point(data = plot_data_line, aes(x = Year, y = MAF)) +
      geom_line(data = plot_data_line, aes(x = Year, y = MAF)) +
      theme_classic() +
      scale_x_reverse(expand = expansion(add = res)) +
      scale_y_continuous(expand = expansion(mult = c(0,0.01)), limits = c(0,1.1)) +
      ggtitle(label) +
      ylab(paste("MAF (",minor_lab,")",sep="")) +
      theme(text = element_text(size = 22 )) +
      theme(plot.title = element_text(hjust = 0.5,size=30)) +
      theme(axis.title.x = element_blank(), 
            axis.text.x = element_blank()) + #remove the y axis label
      theme(axis.ticks.length.x = unit(0, "cm")) + # set the length of the divider of value in the x axis to 0
      theme(plot.margin = unit(c(0,0,0,0), "cm")) # expand the plot so when ggarange there is no gap in between two plots
    
    # Plot stacked bar plot
    plot2 <- ggplot() +
      geom_bar(data = plot_data_bar_long,stat = "identity", aes(x=Year,y=Count,fill=Allele)) + 
      scale_fill_manual(values = c("Minor"="#f6b26b","Major"="#2986cc")) +
      scale_x_reverse(expand = expansion(add = c(500, 500))) + # expand the bar and also reverse the x axis value
      scale_y_reverse(expand = expansion(mult = c(0.05, 0))) + # expand the y axis upward and downward so the x-axis kinda disappear when merge two plot
      theme_classic() + theme(legend.position = "right") +
      theme(text = element_text(size = 22)) +
      theme(plot.margin = unit(c(0,0,0,0.09), "cm")) # Expand the plot so when merge it seamless
    
    # Combine plots
    plot3 <- egg::ggarrange(plot1, plot2, heights = c(0.25, 0.4), ncol=1)
    
    # Return result
    return(plot3)
  }
  
  newped <- reactive({
    # Use function 0
    req(input$bim, input$ped, input$snp)
    
    # Error handling
    newped <- tryCatch({
      make_ped(snp = input$snp, bim = input$bim, ped = input$ped)
    }, error = function(e) {
      return("error1")
    })
    newped
  })
  
  primary <- reactive({
    # Use function 1
    req(input$txt, input$fam)
    
    # Error handling
    primary <- tryCatch({
      make_primary_data(txt_file = input$txt, fam_file = input$fam)
    }, error = function(e) {
      return("error2")
    })
    primary
  })
  
  time_vector <- reactive({
    # Use function 2
    req(input$step)
    
    # Error handling
    if (is.null(primary())){
      time_vector = "error3"
    } else if (primary()=="error2") {
      time_vector = "error2"
    } else {
      time_vector <- create_time_vec(step = input$step, data = primary())
    }
    time_vector
  })
  
  ped_list <- reactive({
    # Use function 3
    # Error handling
    ped_list <- tryCatch({
      make_secondary_data(time_vec = time_vector(), 
                          ped = newped(), primary = primary())
    }, error = function(e) {
      return(NULL)
    })
    ped_list
  })
  
  plot_data <- reactive({
    # Use function 4
    req(input$snp, input$bim)
    # Error handling
    plot_data <- tryCatch({
      allele_data(snp = input$snp, bim = input$bim, ped = ped_list(),
                  time_vec = time_vector())
    }, error = function(e) {
      return(NULL)
    })
    plot_data
  })
  
  plot <- reactive({
    # Use function 5
    # Error handling
    plot <- tryCatch({
      plotting_maf(plot_data(), time_vector())
    }, error = function(e) {
      return(NULL)
    })
    plot
    
    error <-reactive({
      error=NULL
      if (newped()=="error1"){
        error=c("1. Input SNP could not be found in .bim file",
                "2. Wrong input or not right format of .bim or .ped file")
      } else if(primary()=="error2"){
        error=c("Wrong input or not right format of.txt or .fam file")
      } else if(time_vector()=="error3"){
        error=c("No common ID was found in between .txt and .fam file")
      }
      error
    })
  })
 
  #OUTPUT 
  output$output <- renderPlot({
    if (input$submit > 0) {
      if (!is.null(plot())){
        print(plot())
      } else {
        plot(c(1,2),c(1,2), xlab = "This is an error display plot", ylab = "Same as x-axis")
        text(x=1.5,
             y=1.9,
             labels="MAF plot can not be made because error happened",
             col="red",
             cex=1)
        text(x=1.5,
             y=1.5,
             labels=error)
      }
    }
  }) #output$output
  
  output$download <- downloadHandler(
    filename = function() {
      paste0(input$snp,".png")
    },
    content = function(file){
      png(file,
          height = 2048,
          width = 2048,
          res = 300)
      print(plot())
      dev.off()
    }) #output$download
  
} #server

shinyApp(ui, server)