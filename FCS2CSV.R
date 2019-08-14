require(flowCore)

bin <- commandArgs(TRUE)

bincomtudo <- read.FCS(paste(bin ,".fcs", sep=""))
#cotovelodophelps <- data.frame(exprs(bin))
out <- paste(bin,".txt", sep="")
write.csv(data.frame(exprs(bincomtudo)),file = out)