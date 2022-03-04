# load libraries
library(tidyverse)

# CHOOSE SET OF INPUT FILES TO LOOP RUN -----------------------------------
data_path <- "Astrangia"   # path to the data
results_files <- dir(data_path, pattern = "*.csv") # get file names
names <- sub("\\.csv.*", "", results_files)

# Load all result files from folder
for(i in names){
  filepath = file.path(data_path,paste(i,".csv", sep="")) # sets up .csv files
  assign(i, read.delim(filepath, sep = ",") %>% # loads up files
    mutate(mutated_p = -log(pvalue), # modifies to a signed p-value
           mutated_p_updown = ifelse(log2FoldChange < 0,
                                     mutated_p*-1,
                                     mutated_p*1)) %>%
    select(X, mutated_p_updown) %>%
    rename(GeneID = X) %>%
    na.omit() 
    )
}

# Write these new modified files
for( i in 1:length(names)) {
  write.table(get(names[i]), 
             paste0(names[i],
             "_modified_pvalues.csv"),
             sep = ",",
             row.names = FALSE, 
            col.names = FALSE,
            quote = FALSE)
}

# Now we can loop our GO analyses

# SET BASIC VARS ---------------------------------------------------------
goDatabase = "go.obo"
goAnnotations = "astrangia_iso2go.tab"
divisions = c('CC', 'MF', 'BP')
source("gomwu.functions.R")

####################  OVERALL SETS #################### 
inputFiles = paste0(names, "_modified_pvalues.csv")

# LOOP THROUGH SELECTED INPUT FILES AND GO DIVISIONS ---------------------------------------
for (goDivision in divisions){
  print('--------------')
  print('--------------')
  print("------WORKING ON A NEW GO CATEGORY-----")
  print(goDivision)
  print('--------------')
  print('--------------')
  for (input in inputFiles){
    print('--------------')

    print("input being worked on")
    print(input)
    #print(paste0(input, '...', sep='.'))
    gomwuStats(input, goDatabase, goAnnotations, goDivision,
               perlPath ="perl", # replace with full path to perl executable if it is not in your system's PATH already
               largest =0.1,  # a GO category will not be considered if it contains more than this fraction of the total number of genes
               smallest = 10,   # a GO category should contain at least this many genes to be considered
               clusterCutHeight = 0.25, # threshold for merging similar (gene-sharing) terms. See README for details.
               # Alternative=ALTERNATIVE # by default the MWU test is two-tailed; specify "g" or "l" of you want to test for "greater" or "less" instead.
               # Module=TRUE,Alternative="g" # un-remark this if you are analyzing a SIGNED WGCNA module (values: 0 for not in module genes, kME for in-module genes). In the call to gomwuPlot below, specify absValue=0.001 (count number of "good genes" that fall into the module)
               #	Module=TRUE # un-remark this if you are analyzing an UNSIGNED WGCNA module 
    )
  }
}




# record how many significant and save thsoe ---------------------------------------------
library(tidyverse)
divRec = c()
inRec = c()
sigRec = c()
for (goDivision in divisions){
  for (input in inputFiles){
    resName = paste(paste('MWU', goDivision, sep = "_"), input, sep = "_")
    go.res = read.table(resName, header = T)
    totSig = sum(go.res$p.adj < 0.1)
    divRec = append(divRec, goDivision)
    inRec = append(inRec, input)
    sigRec = append(sigRec, totSig)
    sig=go.res[go.res$p.adj < 0.1,] %>% 
      arrange(pval)
    sigOut=paste( c('./resultsFiles/', goDivision, input, '_sigGos.tsv'), collapse='')
    if (nrow(sig)>0){
      sig %>% 
        write_tsv(path=sigOut)
    }
  }
}
res = tibble('goDivision'=divRec,
             'input'=inRec,
             'nSig'=sigRec)
res %>% 
  write_tsv(path='./resultsFiles/gomwu_results_summary.tsv')



# OUTPUT THE HIGH-STRESS GO RESULTS -------------------------------------------
input='corStress_For_MWU.csv'
goDivisions=c('BP', 'MF', 'CC')
go.res = data.frame()
for (goDivision in goDivisions){
  print(goDivision)
  resName = paste(paste('MWU', goDivision, sep = "_"), input, sep = "_")
  divDf = read.table(resName, header = T)
  divDf$goDivision=goDivision
  go.res=rbind(go.res, divDf)
}
go.res %>% 
  arrange(goDivision, pval) %>% 
  write_tsv('../results_tables/all_high_stress_go_mwu_results.tsv')


#all high stress
input='corStress_For_MWU.csv';goDivision='BP';LVL=1e-9
input='corStress_For_MWU.csv';goDivision='MF';LVL=1e-3
input='corStress_For_MWU.csv';goDivision='CC';LVL=1e-3

results=gomwuPlot(input,goAnnotations,goDivision,
                  # absValue=0.0001,
                  absValue=-log(0.05, 10),  # genes with the measure value exceeding this will be counted as "good genes". Specify absValue=0.001 if you are doing Fisher's exact test for standard GO enrichment or analyzing a WGCNA module (all non-zero genes = "good genes").
                  level1=LVL, # FDR threshold for plotting. Specify level1=1 to plot all GO categories containing genes exceeding the absValue.
                  level2=LVL, # FDR cutoff to print in regular (not italic) font.
                  level3=LVL, # FDR cutoff to print in large bold font.
                  txtsize=1.2,    # decrease to fit more on one page, or increase (after rescaling the plot so the tree fits the text) for better "word cloud" effect
                  treeHeight=0.5, # height of the hierarchical clustering tree
                  #	colors=c("dodgerblue2","firebrick1","skyblue2","lightcoral") # these are default colors, un-remar and change if needed
)




# plot results for each module with any significant enrichment -------------------
sig = res %>% 
  filter(nSig>=2)

for (i in 1:nrow(sig)){
  row=sig[i,]
  goDivision=row['goDivision']
  input = row['input']
  figFileName = paste('./resultsFiles/', sep='', paste(paste(goDivision, input, sep='_'), 'tree.pdf', sep='_'))
  pdf(figFileName)
  gomwuPlot(input,goAnnotations,goDivision,
            absValue=0.00001,  # genes with the measure value exceeding this will be counted as "good genes". Specify absValue=0.001 if you are doing Fisher's exact test for standard GO enrichment or analyzing a WGCNA module (all non-zero genes = "good genes").
            level1=0.05, # FDR threshold for plotting. Specify level1=1 to plot all GO categories containing genes exceeding the absValue.
            level2=0.05, # FDR cutoff to print in regular (not italic) font.
            level3=0.001, # FDR cutoff to print in large bold font.
            txtsize=1.2,    # decrease to fit more on one page, or increase (after rescaling the plot so the tree fits the text) for better "word cloud" effect
            treeHeight=0.5, # height of the hierarchical clustering tree
            #	colors=c("dodgerblue2","firebrick1","skyblue","lightcoral") # these are default colors, un-remar and change if needed
  )
  dev.off()
}


#Code below is included in separate script look_at_genes_within_GO_terms.R

#------------ LOOK AT GENES WITHIN GO TERMS ------------#
library(cowplot)
library(tidyverse)
#upload the results for each genes output by the functions above
#re-pick input if you need to
input = 'brown_moduleInput.csv'
goDivision = 'CC'

#upload GO enrichment results
resName = paste(paste('MWU', goDivision, sep = "_"), input, sep = "_")
go.res=read.table(resName, header = T)
go.res = go.res %>% 
  arrange(pval)
sig = go.res[go.res$p.adj < 0.1,]
head(sig, n=20)


#upload the gene-GO associations
geneResName=paste(goDivision, input, sep='_')
gene.res=read.table(geneResName, header=T)
head(gene.res)


#search the GO results for a particular term
searchString = 'ribo'
sig[grep(searchString, sig$name),] 


#select the GO term from your results file (sigGo above)
go="GO:0006979" #oxidative stress
go="GO:0000302" #response to reactive oxygen species
go='GO:0004930' #G-protein coupled receptor activity
go='GO:0007409;GO:0048812' #neuron projection morphogenesis
go='GO:0043009;GO:0009792;GO:0009790'

#subset for that GO term
go.genes = gene.res[gene.res$term == go, 'seq']
length(go.genes)

#get gene names
geneSet = as.character(go.genes)

#GATHER GENE NAMES
adat = read.table('../../metadata/Amillepora_euk.emapper.annotations.tsv',
                  sep='\t',
                  header = TRUE)
rownames(adat) = adat$query_name
annots = adat[geneSet,] %>% 
  dplyr::select(query_name,eggNOG.annot) %>% 
  dplyr::rename(gene=query_name)
nrow(annots)
length(geneSet)


#merge with deseq results
ll=load('../../correlated_only/deseqResults/stress_deseqResults.Rdata')
ll
mdat = res %>% 
  data.frame() %>% 
  mutate(gene = rownames(res)) %>% 
  right_join(annots, by = 'gene')
head(mdat)
geneSet = mdat$gene
geneNames=mdat$eggNOG.annot

#check the log2 values match expectation from heatmap
mdat %>% 
  ggplot(aes(y=log2FoldChange)) +
  geom_boxplot()


#------------ BUILD HEATMAP FOR A GIVEN GO TERM ------------#
library(pheatmap)


#you'll need the variance stabilized counts
ll=load('../../largeIgnored/strongStressOnly_project_controlled.Rdata')
ll
rld.df=data.frame(t(datExpr))
head(rld.df)


#subset the variance stabilized counts to include only the GO of interest
g=rld.df[geneSet,]
g[1:10,1:10]
cdat2 = coldata[colnames(g),] %>% 
  mutate(stress=if_else(treat=='control',
                        'control',
                        'stress'))
cdat2=cdat2[with(cdat2,order(cdat2$stress, cdat2$my_title)),]
orderedRuns = cdat2$Run
orderedTreats = cdat2$stress
orderedG = g[,orderedRuns]

pheatmap(test,
         cluster_cols=FALSE,
         cluster_rows=TRUE,
         labels_col=paste(orderedTreats, orderedRuns, sep='.')[1:10])






