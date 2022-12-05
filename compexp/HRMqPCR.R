##Functions to analyse HRM and qPCR data


loadHrmDataFile <- function(filename) {
    #' Load the data from an HRM table (csv format) exported from the BioRad
    #' software.
    #' It assumes columns with "X" in their name are to be dropped.
    data <- read.csv(filename, header = T)
    ## Drop columns with "X" in their name
    data <- data[, !grepl("X", names(data))]
    return(data)
}

#I modified Matthieu's function, find maximum of negative first derivative (makes a nicer plot)
findMeltingTempOneSample <- function(temperature, rawRfu, range = c(77, 83),
                                    step = 1e-4) {
    #' Find the melting temperature for a sample
    #'
    #' @param temperature Vector with the temperatures of the melting curve
    #' @param rawRfu Vector with the rfu data of the melting curve
    #' @param range Range in which the minimum of the first-order derivative is
    #'   searched for
    #' @param step Step used for the interpolation grid
    f <- splinefun(temperature, rawRfu)
    x <- seq(range[1], range[2], by = step)
    d <- -f(x, deriv = 1) #Negative first derivative
    out <- x[which(d == max(d))] #Find the maximum
    if (length(out) > 1) return(NA)
    return(out)
}


findMeltingTemp <- function(hrmData, samples = NULL, range = c(77, 83), step = 1e-4) {
    #' Determine the melting temperature for all samples in a HRM data table
    #'
    #' @param hrmData Data frame with raw RFU data, loaded with
    #'   loadHrmFile() for example
    #' @param samples If NULL, calculate temperature for all samples, otherwise
    #'   only for specified samples
    #' @param range Range in which the minimum of the first-order derivative is
    #'   searched for
    #' @param step Step used for the interpolation grid
    sampleNames <- names(hrmData)[names(hrmData) != "Temperature"]
    tempValues <- rep(NA, length = length(sampleNames))
    for (i in 1:length(sampleNames)) {
        temp <- hrmData[, "Temperature"]
        sampleData <- hrmData[, sampleNames[i]]
        tempValues[i] <- findMeltingTempOneSample(temp, sampleData, range, step)
    }
    names(tempValues) <- sampleNames
    if (is.null(samples)) {
        return(tempValues)
    } else {
        return(tempValues[samples])
    }
}


#Function to normalize RFU data
RFUnorm <- function(RFUdat) {
    min.rfu <- min(RFUdat) #Get the minimum
    minsub <- RFUdat - min.rfu #Substract minimum
    RFU.norm <- minsub/(max(minsub)) #Divide subtracted by maximum
    return(RFU.norm) #RFU from 1 to 0
}

#This is not needed anymore
comp.1.deriv <- function(ycol, xcol) {
res.vec <- rep(0, length(ycol))
for(i in 1:(length(ycol)-1)) {
    #Compute negative first order derivative
    y1 <- ycol[i]
    y2 <- ycol[i+1]
    x1 <- xcol[i]
    x2 <- xcol[i+1]
    k <- (y2 - y1)/(x2 - x1)
    res.vec[i] <- -k
}
return(res.vec)
}


##This function converts obeserved y to proportion
##Then use x = (y - a)/B to estimate unknown proportions
#Function to calculate proportion from standard curve, input posterior dist. returns distribution
convert2prop <- function(post, y) {
    prop <- (y - post$a) / post$B
    prop[prop > 1] <- 1 #Proportion is between 1 and 0
    prop[prop < 0] <- 0
    return(prop)
}
