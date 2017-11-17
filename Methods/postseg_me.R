
postseg_me = 

function (segmentData, inter = c(-0.1, 0.1)) 
{
  seg <- segmented(segmentData)
  matrixValues <- matrix(0, ncol = ncol(seg), nrow = nrow(seg))
  for (i in 1:ncol(seg)) {
    aux = seg[seg[, i] > -0.1 & seg[, i] <  0.1,i]
    
    if(length(aux)==0) {matrixValues[,i] = 0}
    else
    { matrixValues[,i] = median(aux)}
  }

  seg <- seg - matrixValues
  countlevall <- apply(seg, 2, function(x) {
    as.data.frame(table(x))
  })
  intcount <- function(int, sv) {
    sv1 <- as.numeric(as.vector(sv[, 1]))
    wh <- which(sv1 <= int[2] & sv1 >= int[1])
    return(sum(sv[wh, 2]))
  }
  postsegnorm <- function(segvec, int = inter, intnr = 3) {
    intlength <- (int[2] - int[1])/2
    gri <- intlength/intnr
    intst <- int[1] + (0:intnr) * gri
    intend <- intst + intlength
    ints <- cbind(intst, intend)
    intct <- apply(ints, 1, intcount, sv = segvec)
    whmax <- which.max(intct)
    return(ints[whmax, ])
  }
  postsegnorm_rec <- function(segvec, int, intnr = 3) {
    newint <- postsegnorm(segvec, int, intnr)
    newint <- postsegnorm(segvec, newint, intnr)
    newint <- postsegnorm(segvec, newint, intnr)
    newint <- postsegnorm(segvec, newint, intnr)
    newint <- postsegnorm(segvec, newint, intnr)
    return(newint[1] + (newint[2] - newint[1])/2)
  }
  listres <- lapply(countlevall, postsegnorm_rec, int = inter)
  vecres <- c()
  for (i in 1:length(listres)) {
    vecres <- c(vecres, listres[[i]])
  }
  segmented(segmentData) <- t(t(seg) - vecres)
  copynumber(segmentData) <- t(t(copynumber(segmentData) - 
                                   matrixValues) - vecres)
  return(segmentData)
}
