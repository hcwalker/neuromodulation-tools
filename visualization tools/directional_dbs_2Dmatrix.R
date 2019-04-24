## visualize a single variable across directional DBS contacts 
## in a scaleable, square 2D color matrix

install.packages("dplyr")
install.packages("plot3D")
library(dplyr)
library(plot3D)
setwd("C:/R processes")

{
directional_rows <- 2 # number of directional rows
directional_segments <- 4 # number of directional contact segments per row
participants <- 10 # number of participants
mydata <- c(1:(directional_rows*directional_segments*participants)) # insert a vector with the data of interest sorted by participant and contact number
spacer <- rep(NA, directional_segments) # spacer with number of NAs equal in length to directional segments
mydata_spaced <- mydata 
for (i in 1:(participants)) { # append spaces between each participant
  mydata_spaced <- append(mydata_spaced, spacer, 
  after=((directional_rows*directional_segments*i)+((i-1)*length(spacer))))
}
mydata_spaced
x <- 1:(directional_rows*participants+(participants-1)) # dimensions for initial matrix
y <- 1:directional_segments
mydata_matrix_long <- matrix (nrow=max(y), ncol=max(x)+1, data=mydata_spaced)
mydata_df_long <- mydata_matrix_long %>% t() %>% as.data.frame() # transpose to long format
mydata_df_long
empty_rows <- which(is.na(mydata_df_long[,1])) # identify empty rows
dividers <- ceiling((length(empty_rows))^0.5) # define indices to slice the data
dimensions <- (ceiling(participants^0.5))^2 # dimensions of a square output matrix
padding <- dimensions - participants # check to see if padding is needed
if (padding > 0) { # if padding is needed
  for (i in 1:padding) { # iterate through number of pads needed
    for (i in 1:(directional_rows + 1)) { # insert NAs for each pad
      mydata_df_long[nrow(mydata_df_long)+1,] <- c(rep(NA,directional_segments))
    }
  }
}
empty_rows <- which(is.na(mydata_df_long[,1])) # identify empty rows 
dividers <- seq(empty_rows[dividers],max(empty_rows),by=empty_rows[dividers]) # identify empty rows for slicing the dataframe
mydata_df_long["space"] <- "NA" # add padding between columns of participants
mydata_df_wide <- slice(mydata_df_long, (1:dividers[1])) # first iteration wide format dataframe
for (i in 1:(length(dividers)-1)) { # subsequent iterations make wide format dataframe
  temporary <- slice(mydata_df_long, (dividers[i]+1):dividers[i+1])
  mydata_df_wide <- cbind(mydata_df_wide, temporary)
  }
mydata_df_wide <- mydata_df_wide[-nrow(mydata_df_wide),] # remove empty padding row
mydata_df_wide <- mydata_df_wide[,-ncol(mydata_df_wide)] # remove empty padding column
y <- 1:nrow(mydata_df_wide) # identify dimensions of new matrix
x <- 1:ncol(mydata_df_wide)
mydata_df_wide <- mydata_df_wide %>% t() %>% as.numeric() # transpose and change to numeric
mydata_matrix <- matrix(nrow=max(x), ncol=max(y), data=mydata_df_wide)
empty_rows <- which(is.na(mydata_matrix[,1])) # identify empty rows
empty_rows_diff <- empty_rows %>% diff() %>% append(FALSE, after=0) # identify consecutive rows with padding NAs, add an element for symmetry
empty_rows <- data.frame(empty_rows,empty_rows_diff)
empty_rows <- empty_rows %>% 
  dplyr::filter(empty_rows_diff==1) %>% 
  dplyr::select(empty_rows) %>% 
  unlist(use.names = FALSE)
if (length(empty_rows)>0) { # remove padding rows
  mydata_matrix <- mydata_matrix[-empty_rows,]
  mydata_matrix <- mydata_matrix[-nrow(mydata_matrix),]
}
par(mfrow = c(1, 1))
par(fg="white") # plot the data
image2D(mydata_matrix, border="white", 
        main="therapeutic window (mA)",
        axes=FALSE, xlab=NA, ylab=NA)
par(fg="black")
}
