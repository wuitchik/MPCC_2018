## Pipeline to call snips from O.arbuscula Tagseq data for the purpose of assessing clonality of samples

# Working directory

/projectnb/coral/MPCC_2018/oculina_nov_2021/SNP_analysis

# Modules

module load python3
module load bowtie2
module load htslib/1.9
module load samtools/1.9
module load angsd
module load picard

# creating bowtie2 index for your genome
bowtie2-build coral_refgenome.fasta coral_refgenome.fasta

# trimming
for file in *.fastq
    do echo  "../../tagseq_files/tagseq_clipper.pl $file | fastx_clipper -a AAAAAAAA -l 20 -Q33 | fastx_clipper -a AGATCGGAAG -l 20 -Q33 | fastq_quality_filter -Q33 -q 20 -p 90 > ${file/.fastq/}.trim" >> trimming
done 

../../scc6_qsub_launcher.py -N trim -P coral -M wuitchik@bu.edu -j y -h_rt 24:00:00 -jobsfile trimming
qsub trim_array.qsub

# mapping 

for file in *.trim
    do echo "bowtie2 -x coral_refgenome.fasta -U $file --local -p 4 -S ${file/.trim/}.sam">> maps
done

../../scc6_qsub_launcher.py -N mapping -P coral -M wuitchik@bu.edu -j y -h_rt 24:00:00 -jobsfile maps
qsub mapping_array.qsub

#converting .sam to .bam with proper sorting

for file in *.sam
do echo "samtools sort -O bam -o ${file/.sam/}.bam $file" >> sortConvert
done 

/projectnb/coral/MPCC_2018/scc6_qsub_launcher.py -N bamming -P coral -M wuitchik@bu.edu -j y -h_rt 24:00:00 -jobsfile sortConvert
qsub bamming_array.qsub

# picard index files

picard CreateSequenceDictionary R=coral_refgenome.fasta O=coral.dict

# samtool index files 

samtools faidx coral_refgenome.fasta 

samtools view P3.bam | head

# running angsd

FILTERS="-uniqueOnly 1 -remove_bads 1 -minMapQ 20 -baq 1 -ref coral_refgenome.fasta -maxDepth 470" 
TODO="-doQsDist 1 -doDepth 1 -doCounts 1 -dumpCounts 2"

for file in *.bam
do echo "$file" >> bam_list
done 

angsd -b bam_list -GL 1 $FILTERS $TODO -P 8 -out dd 

#Make new file called plotQC.R and copy script from Misha's website (https://github.com/z0on/2bRAD_denovo/blob/master/plotQC.R) into it

module load R
Rscript plotQC.R dd

A4.bam       0.03824740
K1.bam       0.05653386
R9.bam       0.05936636
N6.bam       0.06758977
J15.bam      0.06973381
M3.bam       0.07146659
C5.bam       0.07156565
K3.bam       0.07218629
A5.bam       0.07387817
L8.bam       0.07921364
D6.bam       0.07941149
F7.bam       0.07980265
D4.bam       0.08143480
R7.bam       0.08534722
M2.bam       0.08631755
E11.bam      0.08670974
I3.bam       0.08808522
Q1.bam       0.08935118
P3.bam       0.09012561
L6.bam       0.09387959
H1.bam       0.09403230
Q11.bam      0.09427447
H11.bam      0.09485621
C4.bam       0.09550095
E3.bam       0.09579730
K2.bam       0.09773928
C9.bam       0.10110722
F9.bam       0.10296546
O1.bam       0.10299229
O2.bam       0.10384705
P1.bam       0.10485341
A6.bam       0.10544847
N4.bam       0.11063041
F8.bam       0.11079286
M1.bam       0.11243517
R8.bam       0.11320965
P2.bam       0.11479815
L7.bam       0.11485547
E10.bam      0.11583710
J13.bam      0.11635389
I2.bam       0.12188022
H15.bam      0.12375534
I1.bam       0.12483439
D5.bam       0.13159315
J14.bam      0.14367725
Q4.bam       0.15549067
N5.bam       0.17563195

# Same filters as Hannah
FILTERS="-uniqueOnly 1 -remove_bads 1 -minMapQ 20 -minQ 25 -dosnpstat 1 -doHWE 1 -sb_pval 1e-5 -hetbias_pval 1e-5 -skipTriallelic 1 -minInd 36 -snp_pval 1e-5 -minMaf 0.05"
TODO="-doMajorMinor 1 -ref coral_refgenome.fasta -doMaf 1 -doCounts 1 -makeMatrix 1 -doIBS 1 -doCov 1 -doGeno 8 -doBcf 1 -doPost 1 -doGlf 2"

angsd -b bam_list -GL 1 $FILTERS $TODO -P 8 -out coral_ang

# ADMIXTURE

module load admixture
NSITES=`zcat coral_ang.beagle.gz | wc -l`
echo $NSITES

## 33939

for K in 2 3 4; 
do 
NGSadmix -likes coral_ang.beagle.gz -K $K -P 10 -o admix;
done