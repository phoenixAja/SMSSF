fh <- read.csv("hdac8_consensus.txt", header=FALSE, stringsAsFactors = FALSE, sep=" ")


Vina <- na.omit(as.data.frame(cbind(as.numeric(fh$V1), as.numeric(fh$V3))))
AD4 <- na.omit(as.data.frame(cbind(as.numeric(fh$V1), as.numeric(fh$V4))))
OpenEye_FredGauss_npr <- na.omit(as.data.frame(cbind(as.numeric(fh$V1), as.numeric(fh$V5)))) 
OpenEye_HybridGauss_npr <- na.omit(as.data.frame(cbind(as.numeric(fh$V1), as.numeric(fh$V6)))) 
OpenEye_FredGauss_MOE <- na.omit(as.data.frame(cbind(as.numeric(fh$V1), as.numeric(fh$V7))))                
OpenEye_HybriddGauss_MOE <- na.omit(as.data.frame(cbind(as.numeric(fh$V1), as.numeric(fh$V8))))  

prgrms <- list(Vina, AD4, OpenEye_FredGauss_npr, OpenEye_HybridGauss_npr, OpenEye_FredGauss_MOE, OpenEye_HybriddGauss_MOE)

for (i in prgrms){
  names <- c("vina", "AD4", "OE_FG_npr", "OE_HG_npr", "OE_FG_MOE", "OE_HG_MOE")
  for (j in seq_along(names)){
    nfh <- change_file(i)
    write.table(nfh, file=paste0(names[j],"Processed.dat"), row.names = FALSE, col.names = FALSE, sep=" ")
  }
}                                                                     
                                               
replace_negs <- function(energy){
  pos_free_energy <- energy*-1
  return(pos_free_energy)
}

change_file <- function(prgrm_fh){
  neg_E <- as.numeric(as.character(prgrm_fh[,2]))
  pos_E <- sapply(neg_E, replace_negs)
  prgrm_fh[,2] <- NULL
  new_frame <- cbind(prgrm_fh, pos_E)
}
