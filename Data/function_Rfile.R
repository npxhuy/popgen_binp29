# Create a main data
# Take txt file in the working directory (that has the ID of individuals and its mean date)
# Take the .fam file in the working directory
# Filter out the similar ID between those two and return a filter data of ID and mean date
make_primary_data <- function(){
  # Read txt file, this file has ID and mean date
  date_data=read.table(list.files(pattern = "\\.txt$"),header=FALSE)              # OPTIONAL TO TAKE USER INPUT, write if else for the correct format of txt file
  # Make sure that is a numeric data
  date_data$V2 <- as.numeric(date_data$V2) 
  
  # Read fam file
  fam_data=read.table(list.files(pattern = "\\.fam$"),header=FALSE)               # OPTIONAL TO TAKE USER INPUT
  filter_data = unique(date_data[date_data$V1 %in% fam_data$V2,])
  return(filter_data)
}

# Create a time step vector function, take the step as an number and data as the
# main data that have the mean date of individuals and its ID
create_time_vec <- function(step,data){
  round_digit=nchar(max(data$V2))-2
  end_time=round(max(data$V2),-round_digit)
  start_time=0
  time_vec=seq(start_time,end_time+step,step)
  return(time_vec)
}

#Example run
time_vec <- create_time_vec(step = 2000,
                            data = make_primary_data())

# Make main data for calculate MAF
# Take time_vec from the create_time_vec function
# ped is the ped file in the working directory
# primary is from the make_primary_data function
make_secondary_data <- function(time_vec,ped,primary){
  #time_data_list=list()
  ped_data_list=list()
  for (i in 1:(length(time_vec)-1)) { 
    # Create time data, add to a time_data list
    time_temp=primary[primary$V2 >= time_vec[i] & primary$V2 < time_vec[i+1], ]
    # time_data_list[[i]]=time_temp                                             # I dont think this data is important tho, maybe can add it later to the function as return
    
    # Separate ped data, 
    ped_temp=ped[ped$V2 %in% time_temp$V1,]
    ped_data_list[[i]]=ped_temp
  }
  return(ped_data_list)
}

#Example run
ped_data_list <- make_secondary_data(time_vec = time_vec,
                                     ped = ped_data,
                                     primary = make_primary_data())


# The bim is the .bim from the working directory
# The ped1 is the original ped file
# The ped2 is the generated file from make_secondary_data
# The time_vec is generated from create_time_vec
allele_data <- function(bim, ped, time_vec){
  snp_name=readline("Input: ")
  snps=bim$V2
  # Ask for user input using readline
  col_index_snp=match(snp_name,snps)+6
  col_index_allele=match(snp_name,snps)
  
  # Retrieve the minor and major
  minor_temp=as.character(bim$V5[col_index_allele])
  major_temp=as.character(bim$V6[col_index_allele])
  
  #Create vector to save the data
  minor_count_vec=NULL
  major_count_vec=NULL
  maf_vec=NULL
  
  for (i in 1:length(ped)){
    # Retrieve ped data
    sep_ped_temp=ped[[i]][,col_index_snp]
    
    # Look for minor and major allele in bim file
    minor_temp=as.character(bim$V5[col_index_allele])
    major_temp=as.character(bim$V6[col_index_allele])
    
    # Paste the allele data to a string
    allele_string=paste(sep_ped_temp, collapse=' ')
    
    # Count minor, major and calculate maf
    # Minor count and append the the current count to a vector
    minor_count=str_count(allele_string, minor_temp)
    minor_count_vec = c(minor_count_vec,minor_count)
    
    # Minor count and append the the current count to a vector
    major_count = str_count(allele_string, major_temp)
    major_count_vec = c(major_count_vec, major_count)
    
    # Calculate maf and append the maf to a vector
    if (minor_count+major_count==0){
      maf=0
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
    major_lab=base_temp
  } else (major_lab=base_lab[result_major])
  
  # Create new time_vec, remove the year 0
  time_vec_temp <- time_vec[-1]
  
  # Return a data frame
  final_data_frame <- data.frame(minor_count_vec,major_count_vec,time_vec_temp,maf_vec)
  colnames(final_data_frame) <- c(minor_lab,major_lab,"Year","MAF")
  
  #Return the allele data
  return(final_data_frame)
}

plotting_maf <- function(data){
  
  #Extract the minor lab for label later
  minor_lab=colnames(data[1])
  major_lab=colnames(data[2])
  
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
  plot1 <- ggplot()+
    geom_point(data=plot_data_line, aes(x = Year, y = MAF))+
    geom_line(data=plot_data_line, aes(x = Year, y = MAF))+
    theme_classic()+
    scale_x_reverse(expand = expansion(add = res)) +
    scale_y_continuous(expand = expansion(mult = c(0,0.01)), limits = c(0,1.1)) +
    theme(axis.title.x = element_blank(), axis.text.x = element_blank()) + #remove the y axis label
    theme(axis.ticks.length.x = unit(0, "cm")) + #set the length of the divider of value in the x axis to 0
    theme(plot.margin = unit(c(0,0,0,0), "cm")) #expand the plot so when ggarange there is no gap in between two plots
  
  # Plot stacked bar plot
  plot2 <- ggplot()+
    geom_bar(data=plot_data_bar_long,stat = "identity", aes(x=Year,y=Count,fill=Allele)) + 
    scale_fill_manual(values = c("Minor"="green","Major"="blue")) +
    scale_x_reverse(expand = expansion(add = c(500, 500))) + #expand the bar and also reverse the x axis value
    scale_y_reverse(expand = expansion(mult = c(0.05, 0))) + #expand the y axis upward and downward so the x-axis kinda disappear when merge two plot
    theme_classic() + theme(legend.position = "none") + #Remove legend
    theme(plot.margin = unit(c(0,0,0,0), "cm")) #Expand the plot so when merge it seamless
  
  # Combine plots
  plot3 <- ggarrange(plot1, plot2, heights = c(0.25, 0.4))
  
  # Return result
  return(plot3)
}


allele_data_1 <- allele_data(bim = bim_data,ped = ped_data_list, time_vec = time_vec )





