# GO analysis (Fishers)

# load libraries
library(tidyverse)


# CHOOSE SET OF INPUT FILES TO loop through from Venn folder -----------------------------------
data_path <- "~/Documents/GitHub/MPCC_2018/Host_analyses/GO_venn"
results_files <- dir("~/Documents/GitHub/MPCC_2018/Host_analyses/Venn_diagrams/", pattern = "*.csv") # get file names
names <- sub("\\.csv.*", "", results_files)
iso2go = read.delim("astrangia_iso2go.tab", sep = "\t")


# Load all result files from folder and mutate to have a 1 for presence, and 0 for absence
for(i in names){
  filepath = file.path("~/Documents/GitHub/MPCC_2018/Host_analyses/Venn_diagrams/",paste(i,".csv", sep="")) # sets up .csv files
  assign(i, read.delim(filepath, sep = ",") %>% # loads up files
           select(Protein_ID) %>%
           mutate(Fishers = 1) %>%
           
          # mutate(mutated_p = -log(pvalue), # modifies to a signed p-value
          #        mutated_p_updown = ifelse(log2FoldChange < 0,
          #                                  mutated_p*-1,
          #                                  mutated_p*1)) %>%
          # select(X, mutated_p_updown) %>%
          # rename(GeneID = X) %>%
          # na.omit() 
  )
}
