#First we start by reading in the bams files.
bams=read.table("Astrangia/bam_list2")[,1] # list of bam files
goods=c(1:length(bams))

ma = as.matrix(read.table("Astrangia/coral_ang.ibsMat"))
ma=ma[goods,goods]
dimnames(ma)=list(bams[goods],bams[goods])
hc=hclust(as.dist(ma),"ave")
plot(hc,cex=0.7) 
