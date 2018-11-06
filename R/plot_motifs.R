# NOTE: Visualization is best in R Markdown. HIGHLY RECOMMENDED TO RUN THIS function ON MARKDOWN.
# NOTE: Currently only for One-Dimensional data.


# Analysis on Time Series is a very 'visual' domain, and it helps a lot if we could visualize the motifs found in the time series.
# And this function does exactly that, giving a series of images of all the motifs found to compare and gain a better understanding of the findings.

plot_motifs <- function(motif_object, actual_data){
  
  
  # @Details: 
  # Plots all the motifs found, i.e. outputs from find_motif() function. Each pair of image represent a motif at two different location according to the 
  # the nearest neighbour index (motif_idx) from the output object of find_motifs() function. 
  # 
  # @Parameters:
  # 1. motif_object : The output object from the function find_motif().
  # 2. actual_data : The vector of which you want to visualize the motifs of.
  
  
  #Checking whether the inputs are of the correct object type: 
  if(typeof(motif_object) != "list"){
    stop("Error: First argument must be an object of class `list`. Input the output object from find_motif function")
  }
  
  #Creating essential parameters: 
  number_of_motifs = length(motif_object$motif$motif_idx)
  window_size = motif_object$w
  
  #Running a for loop for plotting motifs: 
  for (i in 1:number_of_motifs){
    #Taking the motif number for giving the title to the plot
    motif_number = i
    
    for (j in 1:2){
      start <- motif_object$motif$motif_idx[[i]][j]
      end <- motif_object$motif$motif_idx[[i]][j]+window_size-1
      plot(actual_data[start:end], type='o', main = paste("Motif", motif_number, sep = "-"), xaxt="n", xlab = paste("Index starting from", start, sep="-"),
           ylab = "Data")
      #To remove the numbering on the x-axis, as it is misleading/confuses with the Index at which motif are present: 
      axis(1, at = c(window_size), las=2)
    }
  }
}