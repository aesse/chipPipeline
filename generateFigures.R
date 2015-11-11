# setwd("~/Research15/Clustering/MDD")
library(ggplot2)
library(gridExtra)
library(tableplot)
library(RColorBrewer)
require(TxDb.Hsapiens.UCSC.hg19.knownGene)
require(ChIPseeker)

norm <- function(x, y) {
  x/y
}
getPercent <- function(x,y) {
  print(x)
  x <- as.numeric(x)
  y <- as.numeric(y)
  (x/y)*100
}

binary <- function(x) {
  if (x > 0.05) {
    x = 1
  } else {
    x = 0
  }
}

colRamp <- c("#1F78B4","#A6CEE3", "#B2DF8A", "#33A02C",
             "#FB9A99", "#E31A1C", "#FDBF6F", "#FF7F00",
             "#CAB2D6", "#6A3D9A", "#FFFF99", "#B15928",
             "gray50", "navyblue", "chocolate", "hotpink" )

#####################################################################

#######REPLACE BELOW ARGUMENTS###########

cluster = 5
folder = ''
name = 't1'
fileBase = 'wgEncodeH1hescSP1_seg20_500_t1.out'
narrow = 'wgEncodeH1hescSP1.narrowPeak'

#####################################################################

segments = 25
window = 500

segs = rep(0, segments)
for (j in 1:segments) {
  segs[j] <- paste("s",j,sep='')
}
head = rep(0, cluster)
for (j in 0:cluster) {
  head[j+1] <- paste("c",j,sep='')
}



####### Optimal Cluster #########
compList <- read.table(paste(fileBase, "_comp.out", sep=""), header = FALSE, sep='\t')
colnames(compList) <- c("Cluster", "DescriptionLength")
min <- which.min(compList$DescriptionLength) + 1

optimalClus <- ggplot(compList, aes(x=Cluster, y=DescriptionLength)) + 
  geom_point(size=1.2) + geom_line(size=3) +
  geom_vline(xintercept=min, colour='red') +
  scale_x_discrete(limits=seq(2,length(compList$DescriptionLength)+1,1)) +
  theme_bw() + theme(text = element_text(size=45),axis.title.y=element_text(vjust=3.5),axis.title.x=element_text(vjust=.05),plot.title=element_text(vjust=3.5)) +
  theme(plot.margin = unit(c(1,1,1,1), "cm")) +
  ggtitle(paste(name, "- search cluster performance", sep=" ")) +
  xlab("Cluster size") + ylab("Description Length")

######### Plot Clusters #################
data =  read.table(paste(fileBase, "_alpha_", cluster, ".out", sep=""), header = FALSE, sep=',', row.names=1)
totals = apply(data, 1, sum)

orig = as.data.frame(data)
datat = as.data.frame(t(data))
normal = mapply(norm, datat,totals) 
normal = as.data.frame(t(normal))

members = rep(0, cluster)
clusCounts = rep(0, cluster)
for (i in 0:(cluster-1)) {
  dataCounts =  read.table(paste(fileBase, "_bin_",i,"_", cluster, ".out", sep=''), header = FALSE, sep='\t')
  members[i+1] = paste("c",i,"-",dim(dataCounts)[1],sep = "")
  clusCounts[i+1] = dim(dataCounts)[1]
}

for (i in 1:(cluster)) {
  if (i == 1) {
    m <- matrix(c(normal[i,], orig[i,], rep(head[i], segments), segs, rep(totals[i], segments)), ncol=5, nrow=length(segs))
  } else {
    n <- matrix(c(normal[i,], orig[i,], rep(head[i], segments), segs, rep(totals[i], segments)), ncol=5, nrow=length(segs))
    m <- rbind(m, n)
  }
}

dftm <- data.frame(as.numeric(m[,1]), as.numeric(m[,2]), unlist(m[,3]), unlist(m[,4]), as.numeric(m[,5]))
colnames(dftm) <- c("Peak", "oAlpha", "Cluster", "Segments", "Alpha")

plotClusters <- ggplot(dftm, aes(x=Segments, y=Peak, color=Cluster, group = Cluster)) + 
  scale_x_discrete(limits=segs) + geom_line(aes(size=Alpha)) +
  scale_size(range = c(1,7), limits = c(50, 2000000), breaks = c(10000, 250000, 500000, 1000000, 1500000, 1750000, 2000000), name = "AlphaSum") +
  scale_colour_manual(values = colRamp,labels = members) +
  theme_bw() + theme(text = element_text(size=45),axis.title.y=element_text(vjust=3.5),axis.title.x=element_text(vjust=.05),plot.title=element_text(vjust=3.5)) +
  theme(plot.margin = unit(c(1,1,1,1), "cm"), axis.text.x=element_text(angle = 45, hjust=1)) +
  theme(legend.key.size = unit(2, "cm")) +
  guides(colour = guide_legend(override.aes = list(size=3))) +
  ggtitle(paste(name, "-", cluster, "clusters", sep=" ")) +
  ylab("Peak Height (Normalised)") 

#########################################
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

labelFiles <- function(fileBaseClus) {
  files <- c(rep("0", cluster))
  ids <- c(rep("0", cluster))
  for (i in 0:(cluster-1)) {
    # for (i in c(0,1,2,3,5,6,7,9)) {
    print(i)
    fName = paste(fileBaseClus,i,".out", sep='')
    id = paste("c",cluster,"_", i, sep='')
    ids[i+1] = id
    files[i+1] = fName
  }
  
  fs <- as.list(files)
  names(fs) <- ids
  return(fs)
}

tidyAnnotations <- function(peakAnnoList) {
  for (p in 0:(cluster-1)) {
    print(p)
    entryName <- paste('c', cluster, '_', p, sep="")
    entry <- as.matrix(as.data.frame(with(peakAnnoList, get(entryName))))
    clusterAssign = rep(p, dim(entry)[1])
    entry <- cbind(entry, clusterAssign)
    if (p == 0) {
      annotations <- entry
    } else {
      annotations <- rbind(annotations, entry)
    }
  }
  annotations[,9][grepl("Intron",annotations[,9])] <- "Intron"
  annotations[,9][grepl("Exon",annotations[,9])] <- "Exon"
  annotations[,9][grepl("Promoter",annotations[,9])] <- "Promoter"
  annotations[,9][grepl("Downstream",annotations[,9])] <- "Downstream"
  annotations <- data.frame(annotations[,1], as.numeric(annotations[,2]), as.numeric(annotations[,3]), annotations[,5], annotations[,9], as.numeric(annotations[,17]), annotations[,18])
  colnames(annotations) <- c('chrm', 'start', 'end', 'strand', 'annotation', 'distanceToTSS', 'cluster')
  return(annotations)
}

plotPercentages <- function(annotations, plotLabel) {
  tbl <- table(annotations$cluster, annotations$annotation)
  chi <- chisq.test(tbl)
  rTotals <- apply(tbl, 1, sum)
  cTotals <- apply(tbl, 2, sum)
  pvalMat <- matrix(rep(0, length(tbl)), nrow = dim(tbl)[1], ncol=dim(tbl)[2])
  perTbl <- mapply(getPercent, cTotals, sum(cTotals))
  expected <- matrix(c(rep('Exp', length(colnames(tbl))), colnames(tbl),perTbl, rep(1, length(colnames(tbl))),rep(1, length(colnames(tbl))),rep(1, length(colnames(tbl)))), ncol=6)
  percentMat <- matrix(rep(0, length(tbl[1,])*length(tbl[,1])*6), ncol=6)
  k = 1
  for (i in 1:length(rTotals)) { #rows/clusters
    for (j in 1:length(tbl[i,])) { #columns/locations
      percent <- (tbl[i,j] / rTotals[i])*100
      
      rowTotal = as.double(rTotals[i])
      colTotal = as.double(cTotals[j])
      sumCol = as.double(sum(cTotals))
      expVal = (colTotal*rowTotal)/sumCol
      
      a <- tbl[i,j]
      b <- rTotals[i] - a
      c <- cTotals[j] - a
      d <- sum(cTotals) - a - b - c
      conTbl <- matrix(c(a,c,b,d), ncol = 2)
      chi <- fisher.test(conTbl, alternative="t")
      pvalMat[i,j] <- chi$p.value
      
      if(a > expVal) {
        change = "Over"
      } else {
        change = "Under"
      } 
      
      if (a == 0) {
        absCount = "0 Count"
      } else {
        absCount = "Real"
      }
      
      percentMat[k,] <- c(paste('c',i-1,sep=""), colnames(tbl)[j], format(round(percent, 2), nsmall = 2),format(chi$p.value, scientific=T, nsmall=2, digits = 3), change, absCount)
      k = k+1
    }
  }
  percentMat <- rbind(percentMat, expected)
  percent <- data.frame(percentMat[,1], percentMat[,2], as.numeric(percentMat[,3]), as.numeric(percentMat[,4]), percentMat[,5], percentMat[,6])
  colnames(percent) <- c('Cluster', 'Location', 'Percent', "pval", "change", "Counts")
  
  xmid <- (length(rTotals)+1)/2 + 0.5
  genLoc <- ggplot(percent, aes(x=Cluster, y=Percent)) +
    geom_bar(aes(fill=Location), stat="identity") + 
    theme_bw() + theme(text = element_text(size=45),axis.title.y=element_text(vjust=3.5),axis.title.x=element_text(vjust=.05),plot.title=element_text(vjust=3.5)) +
    theme(plot.margin = unit(c(1,1,1,1), "cm")) +
    theme(legend.key.size = unit(2, "cm")) + 
    scale_fill_manual(values = colRamp) +
    ggtitle(plotLabel) +
    annotate("text", x = xmid, y = 105, size=12, label = paste("Chi square p-value: ", format(chi$p.value, digits=3, scientific=T)))
  
  percent = percent[percent$Cluster != "Exp",]
  
  pvalTab <- ggplot(percent[percent$pval < 0.05,], aes(x = Cluster, y = Location, label = format(pval, nsmall = 1), colour=change)) +
    geom_text(size = 12) + 
    theme_bw() + theme(text = element_text(size=45),axis.title.y=element_text(vjust=3.5),axis.title.x=element_text(vjust=.05),plot.title=element_text(vjust=3.5)) +
    theme(plot.margin = unit(c(1,1,1,1), "cm")) + 
    #     theme(legend.position="none") + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    theme(axis.ticks = element_blank()) + theme(legend.key.size = unit(0.5, "cm")) +
    geom_vline(xintercept=seq(1.5, cluster-0.5, 1)) +
    scale_colour_manual(values = c("green","red", "blue")) + 
    geom_text(data = percent[percent$pval >= 0.05,], aes(x = Cluster, y = Location), colour = "black", size = 12) + 
    #     geom_rect(data = percent[percent$count == "0 Count",], aes(x = Cluster, y=Location, colour = Counts)) +
    ggtitle("Individual Fisher tests of each classification")
  
  result <- list("genLoc" = genLoc, "pvalTab" = pvalTab, "percent" = percent)
  return(result)
}

fileBaseClus = paste(fileBase, "_bin_c",cluster,"_", sep="")
cfs <- labelFiles(fileBaseClus)
cPeakAnnoList <- lapply(cfs, annotatePeak, TxDb=txdb, tssRegion=c(-3000, 1000))
cAnnotations <- tidyAnnotations(cPeakAnnoList)
# cPercentRes <- plotPercentages(cAnnotations, paste(name, " - c", cluster, " Genomic Locations", sep=""))
cPercentRes <- plotPercentages(cAnnotations, "Genomic Locations of Binding Sites")
cGenLoc <- cPercentRes$genLoc
cpvalTab <- cPercentRes$pvalTab
cPercent <- cPercentRes$percent

####################################################

pdf(paste(name, "_clustering.pdf", sep=""), width=30, height=45)
grid.arrange(optimalClus, plotClusters, cGenLoc, cpvalTab, ncol=1)
dev.off()
