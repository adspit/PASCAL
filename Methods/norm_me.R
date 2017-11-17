norm_me = function (input, method = "median", smoothOutliers = TRUE, ...) 
{
  values <- c()
  if (method == "none") {
    cat("Skipping normalization ... \n")
  }
  else {
    if (method == "median") {
      cat("Applying median normalization ... \n")
      matrixValues <- matrix(0, ncol = ncol(copynumber(input)), nrow = nrow(copynumber(input)))
      for (i in 1:ncol(input)) {
        aux = copynumber(input)[copynumber(input)[, i] > -0.1 & copynumber(input)[, i] <  0.1,i]
        if (length(aux)==0)
        {matrixValues[,i] = 0}
        else
        {matrixValues[,i] = median(aux)}
      }
    }
     else if (method == "mode") {
       cat("Applying mode normalization ... \n")
       for (i in 1:ncol(input)) {
         density <- density(copynumber(input[, i]))
         value <- density$x[which(density$y == max(density$y))]
         values <- c(values, value)
       }
       matrixValues <- matrix(rep(values, nrow(input)), ncol = ncol(input), 
                              byrow = TRUE)
     }
    # matrixValues <- matrix(rep(values, nrow(input)), ncol = ncol(input), 
    #                        byrow = TRUE)
    copynumber(input) <- copynumber(input) - matrixValues
  }
  if (smoothOutliers) {
    cat("Smoothing outliers ... \n")
    CNA.object <- DNAcopy::smooth.CNA(DNAcopy::CNA(copynumber(input), 
                                                   chromosomes(input), bpstart(input), data.type = "logratio"), 
                                      ...)
    for (i in 1:ncol(input)) {
      copynumber(input)[, i] <- CNA.object[[i + 2]]
    }
  }
  input
}