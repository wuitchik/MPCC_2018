library(adegenet)


#### Host Data ####
countDataHost <- read.table("/Users/hannahaichelman/Dropbox/BU/Host_Buffering/MPCC_2018/Sym_analyses/CountsFiles/Oculina_counts_newref_host.txt")

names(countDataHost) = c("OA4_C_W", "OA5_F_W", "OA6_H_W", 
                         "OC4_F_B",	"OC5_H_B",	"OC9_C_B",	
                         "OD4_C_B",	"OD5_F_B",	"OD6_H_B",
                         "OE10_F_W", "OE11_C_W", "OE3_H_W",
                         "OF7_C_B",	"OF8_F_B",	"OF9_H_B",
                         "OH11_F_W", "OH15_H_W", "OH1_C_W",
                         "OI1_C_B",	"OI2_F_B",	"OI3_H_B",	
                         "OJ13_C_B",	"OJ14_F_B",	"OJ15_H_B",	
                         "OK1_C_W", "OK2_F_W", "OK3_H_W",
                         "OL6_C_B",	"OL7_F_B",	"OL8_H_B",	
                         "OM1_C_B",	"OM2_F_B",	"OM3_H_B",
                         "ON4_C_W", "ON5_F_W", "ON6_H_W",
                         "OO1_C_W", "OO2_F_W",
                         "OP1_C_W", "OP2_F_W", "OP3_H_W",
                         "OQ11_H_W", "OQ1_C_W", "OQ4_F_W",
                         "OR7_C_B",	"OR8_F_B",	"OR9_H_B")

treatmentHost = as.factor(sapply(strsplit(colnames(countDataHost), split = "_"), "[[", 2)) %>%
  revalue(c("C" = "control", "F" = "cold", "H" = "heat"))

genotypeHost  = as.factor(sapply(strsplit(colnames(countDataHost), split = ""), "[[", 2))

symState  = as.factor(sapply(strsplit(colnames(countDataHost), split = "_"), "[[", 3)) %>%
  revalue(c("B" = "sym", "W" = "apo"))

expDesign_Host = data.frame(colnames(countDataHost), treatmentHost, genotypeHost, symState)
expDesign_Host$type = "host"
expDesign_Host$treat_type = paste(expDesign_Host$treatmentHost,expDesign_Host$symState, sep = "_")
names(expDesign_Host) = c("sample", "treatment", "genotype", "sym_state", "type", "treat_type")

colSums(countDataHost)
# remove clones and outlier genotype A from experimental design and counts
expDesign_Host_clones_removed = expDesign_Host %>%
  filter(genotype!="P") %>%
  filter(genotype!="O") %>%
  filter(genotype!="L") %>%
  filter(genotype!="K") %>%
  filter(genotype!="A") #strange on cluster dendrogram

countDataHost_clones_removed = countDataHost %>%
  select(-c(OP1_C_W, OP2_F_W, OP3_H_W, OO1_C_W, OO2_F_W, OL6_C_B, OL7_F_B, OL8_H_B, OK1_C_W, OK2_F_W, OK3_H_W, OA4_C_W, OA5_F_W, OA6_H_W))

# Run DESeq
dds<-DESeqDataSetFromMatrix(countData=countDataHost_clones_removed, colData=expDesign_Host_clones_removed, design=~treat_type) #can only test for the main effects of treatment
dds = DESeq(dds)

# rlog transformation
rlogged = rlog(dds, blind = TRUE)
vst = vst(dds, blind=TRUE)


dat_rlog = as.data.frame(assay(rlogged))
dat_vst = as.data.frame(assay(vst))

# Now run DAPC - look for optimal number of clusters
clus = find.clusters(t(dat_rlog), max.n.clust=30) # 20 pc's and 2 clusters
names(clus)
clus$grp
# freezing corals all in one group, heat and control in the other. lines up well with what we see in the pca 

dapc1 = dapc(t(dat_rlog), clus$grp)

scatter(dapc1, bg = "white", scree.da = FALSE, legend = TRUE, solid = .4)


#### Symbiont in Host Data ####

countDataSym <- read.table("/Users/hannahaichelman/Dropbox/BU/Host_Buffering/MPCC_2018/Sym_analyses/CountsFiles/Oculina_Counts_newref_sym.txt")

names(countDataSym) = c("OA4", "OA5", "OA6", 
                        "OC4_F_B",	"OC5_H_B",	"OC9_C_B",	
                        "OD4_C_B",	"OD5_F_B",	"OD6_H_B",
                        "OE10", "OE11", "OE3",
                        "OF7_C_B",	"OF8_F_B",	"OF9_H_B",
                        "OH11", "OH15", "OH1",
                        "OI1_C_B",	"OI2_F_B",	"OI3_H_B",	
                        "OJ13_C_B",	"OJ14_F_B",	"OJ15_H_B",	
                        "OK1", "OK2", "OK3",
                        "OL6_C_B",	"OL7_F_B",	"OL8_H_B",	
                        "OM1_C_B",	"OM2_F_B",	"OM3_H_B",
                        "ON4", "ON5", "ON6",
                        "OO1", "OO2",
                        "OP1", "OP2", "OP3",
                        "OQ11", "OQ1", "OQ4",
                        "OR7_C_B",	"OR8_F_B",	"OR9_H_B")

countDataSym_brown = countDataSym %>%
  select(-c(OA4, OA5, OA6, OE10, OE11, OE3, OH11, OH15, OH1, OK1, OK2, OK3, ON4, ON5, ON6, OO1, OO2, OP1, OP2, OP3, OQ11, OQ1, OQ4))

# We found a pair of clones in the symbiotic Oculina dataset, M's and L's, so removing same clone pair from the host and other analyses (L's)
countDataSym_brown_noclone = countDataSym_brown %>%
  select(-c(OL6_C_B,  OL7_F_B,  OL8_H_B))

treatmentSym = as.factor(sapply(strsplit(colnames(countDataSym_brown_noclone), split = "_"), "[[", 2)) %>%
  revalue(c("C" = "control", "F" = "cold", "H" = "heat"))

genotypeSym  = as.factor(sapply(strsplit(colnames(countDataSym_brown_noclone), split = ""), "[[", 2))

expDesign_Sym = data.frame(colnames(countDataSym_brown_noclone), treatmentSym, genotypeSym)
expDesign_Sym$type = "sym_inhost"
expDesign_Sym$treat_type = paste(expDesign_Sym$treatmentSym,expDesign_Sym$type, sep = "_")
names(expDesign_Sym) = c("sample", "treatment", "genotype", "type", "treat_type")

dds_Sym<-DESeqDataSetFromMatrix(countData=countDataSym_brown_noclone, colData=expDesign_Sym, design=~treatment) #can only test for the main effects of treatment
dds_Sym = DESeq(dds_Sym)

# rlog transformation
rlogged_syminhost = rlog(dds, blind = TRUE)
vst_syminhost = varianceStabilizingTransformation(dds, blind=TRUE)


dat_rlog_syminhost = as.data.frame(assay(rlogged_syminhost))
dat_vst_syminhost = as.data.frame(assay(vst_syminhost))

# Now run DAPC - look for optimal number of clusters
clus_syminhost = find.clusters(t(dat_rlog_syminhost), max.n.clust=20) # 20 pc's and 2 clusters
names(clus_syminhost)
clus_syminhost$grp
# freezing corals all in one group, heat and control in the other. lines up well with what we see in the pca 

dapc1 = dapc(t(dat_rlog_syminhost), clus_syminhost$grp)

scatter(dapc1, bg = "white", scree.da = FALSE, legend = TRUE, solid = .4)


