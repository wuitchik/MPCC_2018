# Set working directory
setwd("/Users/hannahaichelman/Documents/BU/Host_Buffering/MPCC_2018/Sym_analyses/SNPCalling/Oculina")
#clustering / PCoA based on identity by state (IBS) 
#(for low and/or uneven coverage)


#First we start by reading in the bams files.
bams=read.table("bams_brown")[,1] # list of bam files for symbiotic/brown oculina individuals only
# bams=read.table("bams")[,1] # list of bam files for all oculina individuals
# bams=read.table("bams_noA")[,1] # list of bam files for all oculina individuals except genotype A that looks weird in dendrogram

goods=c(1:length(bams))

ma = as.matrix(read.table("Oarb_browns.ibsMat")) # identity by state matrix for symbiotic/brown oculina individuals only
# ma = as.matrix(read.table("Oarb_newref.ibsMat")) # all oculina individuals
# ma = as.matrix(read.table("Oarb_newref_noA.ibsMat")) # all oculina individuals except genotype A

ma=ma[goods,goods]
dimnames(ma)=list(bams[goods],bams[goods])
hc=hclust(as.dist(ma),"ave")
plot(hc,cex=0.7) 
abline(h=0.17, col='darkgray', lty = 2, lwd=2) # add a vertical line to indicate cut-off for clones
# This cut-off aligns well with what James Fifer has used in the past (0.15), 
# this means that any samples that cluster together below that threshold are clones

# Sample names are read as follows: 
# OI3_H_B = Oculina, genet I, fragment 3, heat treatment, brown/symbiotic
# This analysis indicates that for brown/symbiotic individuals, we only have one unexpected clonal group (M and L)
# When we look at all individuals, we see that genotype A arranges strangely across the dendrogram, we remove from analysis for this reason
# For white/aposymbiotic individuals, we see that there are two more clonal groups: N+O+P and H+K.
