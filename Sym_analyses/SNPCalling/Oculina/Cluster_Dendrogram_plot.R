
setwd("/Users/hannahaichelman/Documents/BU/Host_Buffering/MPCC_2018/Sym_analyses/SNPCalling/Oculina")
#clustering / PCoA based on identity by state (IBS) 
#(for low and/or uneven coverage)


#First we start by reading in the bams files.
bams=read.table("bams_brown")[,1] # list of bam files
goods=c(1:length(bams))

ma = as.matrix(read.table("Oarb_browns.ibsMat"))
ma=ma[goods,goods]
dimnames(ma)=list(bams[goods],bams[goods])
hc=hclust(as.dist(ma),"ave")
plot(hc,cex=0.7) 
abline(h=0.17, col='darkgray', lty = 2, lwd=2) # add a vertical line to indicate cut-off for clones
# This cut-off aligns well with what James Fifer has used in the past (0.15), 
# this means that any samples that cluster together below that threshold are clones
