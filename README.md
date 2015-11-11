# chipPipeline

README

The below commands cover how to run the code using a unix system

Java dependencies:

bnkit_1.8 -> java v1.8

bnkit_1.7 -> java v1.7

Python dependencies:

pysam

pybedtools


R dependencies:
R version 3.2.2
install.packages("ggplot2")
install.packages("gridExtra")
install.packages("tableplot")
install.packages("RColorBrewer")
source("http://bioconductor.org/biocLite.R")
biocLite("Rsamtools")
biocLite('GenomicAlignments')
biocLite("TxDb.Hsapiens.UCSC.hg19.knownGene")
biocLite("ChIPseeker")

-------------------------------------------------------------------------
Example data files:
wgEncodeH1hescSP1.narrowpeak - downloaded from FactorBook. A ChIP-peak file
wgEncodeHaibTfbsH1hescSp1Pcr1xAlnRep1.bam - raw read file used to call ChIP-peaks (must be indexed)
Genome - hg19

-------------------------------------------------------------------------
Example command:
python extractChipPeaks.py wgEncodeH1hescSP1.narrowPeak 500 hg19.fa \
wgEncodeHaibTfbsH1hescSp1Pcr1xAlnRep1.bam \
20 t1 bnkit_1.7.jar 2 12

-------------------------------------------------------------------------
Example output:
t1_clustering.pdf - summary result figures
wgEncodeH1hescSP1_central_500_t1.out - bed file of peaks with locations as 500bp windows around summit
wgEncodeH1hescSP1_central_500_t1.fa - fasta file of peaks with locations as 500bp windows around summit
wgEncodeH1hescSP1_central_500_t1.f.fa - formatted fasta file
wgEncodeH1hescSP1_seg20_500_t1.out - segmented ChIP-peaks with read counts
wgEncodeH1hescSP1_seg20_500_t1.out_comp.out - complexity results for different cluster sizes - used to identify MDL
wgEncodeH1hescSP1_seg20_500_t1.out_alpha_5 - the alpha values describing the resulting Dirichlet clusters
wgEncodeH1hescSP1_seg20_500_t1.out_bin_0_5 - indexes of peaks belonging to this cluster
wgEncodeH1hescSP1_seg20_500_t1.out_bin_c5_0 - bed results for peaks belonging to this cluster

-------------------------------------------------------------------------

To generate the files necessary to perform Dirichlet clustering:

python extractChipPeaks.py chip_result.narrowPeak 500 \
genome.fa chip_raw_reads.bam 20 t1 

-------------------------------------------------------------------------

To generate the files necessary to perform Dirichlet clustering and
to perform the clustering across a range of cluster sizes. Only
the optimal cluster size is reported:

python extractChipPeaks.py chip_result.narrowPeak 500 \
genome.fa chip_raw_reads.bam 20 t1 bnkit.jar 2 10

-------------------------------------------------------------------------

To perform Dirichlet clustering and only record the optimal cluster
size results:

java -jar bnkit.jar chip_result_seg.out 2 10 false

-------------------------------------------------------------------------

To perform Dirichlet clustering and record ALL cluster results:

java -jar bnkit.jar chip_result_seg.out 2 10 true

-------------------------------------------------------------------------

After clustering, the result files for the optimal cluster size
are processed so they can be passed into R to generate figures.

The python script will print a set of variables that must be replaced
in the Rscript (generateFigures.R) to generate the result plots
