################################################################################
#
# basicRnaSeqAnalysis.R calls differentially expressed genes and makes figures
#
# Author:       Fred P. Davis, NIAMS/NIH. fred.davis@nih.gov
#
# This script assumes that kallisto has already been run on the RNA-seq samples
# To run this program, start R and type these commands at the prompt:
#
# R> source("../src/R/basicRnaSeqAnalysis.R")
# R> dat <- paperRun(returnData=TRUE)
# R> makePaperFigures(dat)
#
#
################################################################################


startup <- function() {

   library(grDevices)
   library(Matrix)
   library(readr)
   library(readxl)
   library(magrittr)
   library(dplyr)

   library(calibrate)
   library(pheatmap)

   library(grid)

   library(sleuth)
   library(digest)
   library(ggplot2)

   library(png)

   library(tidyr)
   library(RColorBrewer)
   library(statmod)

}




main <- function(dat, makeFigs = TRUE, returnData = TRUE) {

   startup()

# 0. Load Data

   if (missing(dat))                    dat <- list()

   if (! "specs" %in% names(dat))       dat$specs <- setSpecs()

   if (! "edat" %in% names(dat)) {
      print("Loading expression")
      dat$edat <- loadExpr(specs       = dat$specs,
                           printTables = TRUE)
      if (returnData)    return(dat)
   }

   if (! "deg" %in% names(dat)) {
      print("Loading DEG")
      dat <- defineDEG(dat)
      if (returnData)    return(dat)
   }

   if (makeFigs) {
      makeFigures(dat)
   }

   if (returnData) {
      return(dat)
   } else {
      return(1)
   }

}


makeFigures <- function(dat, figList) {

   if (missing(figList)) {figList <- "all"}

# Table of DEG genes

# Scatterplot of all pairwise comparisons, highlighitng/#'ing DEG's and naming top-10s

# Heatmap of DEG genes
   if (any(c("all","2D") %in% figList)) {
       plotHmap.bulkReporter(dat$dat.bulk,
                             figName    = "2D",
                             dataset    = "reporter",
                             nlMode     = "fracMax",
                             repAvg        = TRUE,
                             clusterRows   = TRUE,
                             decontamEx    = TRUE,
                             show_rownames = FALSE)
   }


# GEO supp expression table
   if (any(c("all","GEO_table") %in% figList)) {
      makeGeoExpressionTable(dat$dat.bulk)
   }

}



loadExpr <- function(specs, printTables = TRUE){

   txExprMat<-NULL
   for (i in 1:nrow(specs$rnaSamples)) {
      print(paste("Loading expression: ",
                  specs$rnaSamples$sampleName[i],
                  date()))
   
#target_id       length  eff_length      est_counts      tpm
#FBtr0082757     1844    1595    9       0.718348
   
      kallistoDir <- specs$kallistoBaseDir
      if (specs$rnaSamples$stranded[i] == "yes") {
         kallistoDir <- specs$kallistoBaseDir.stranded }

      curAbund<-read.table(
         paste(kallistoDir,"/",
               specs$rnaSamples$runID[i],"/",specs$rnaSamples$sampleName[i],
               "/abundance.tsv",sep=""),
         header=TRUE,sep="\t",
         colClasses=c("character","NULL","NULL","numeric","numeric"))

      col1<-paste("est_counts.",specs$rnaSamples$sampleName[i],sep="")
      col2<-paste("tpm.",specs$rnaSamples$sampleName[i],sep="")
   
      colnames(curAbund)<-c( "transcript_id", col1,col2)

      if (is.null(txExprMat)) {
         txExprMat<-curAbund
      } else {
         txExprMat[,col1]<-curAbund[,col1]
         txExprMat[,col2]<-curAbund[,col2]
      }
   }
   txExprMat<-merge(specs$transcriptInfo, txExprMat,
                    all.y=TRUE, by="transcript_id")
   
   allCols<-colnames(txExprMat)
   tpmCols<-allCols[grep("tpm.",allCols)]
   estCountCols<-allCols[grep("est_counts.",allCols)]
   
   # Reorder transcript matrix
   newTxExprMat<-txExprMat[, c("transcript_id", "transcript_name",
                               "gene_id", "gene_name",
                               tpmCols, estCountCols)]
   
   # Format to max 3 decimal digits
   if (printTables) {
      x<-newTxExprMat
      for (i in c(tpmCols, estCountCols)) {
         x[, i]<-sprintf("%.3f", x[, i])}
   
      gz1<-gzfile(paste0(specs$outDir, "/transcriptExpressionMatrix.txt.gz"),
                  "w")
      write.table(x, file=gz1, col.names=TRUE, row.names=FALSE,
                  quote=FALSE,sep="\t")
      close(gz1)
      x<-NULL
   }
   
   
   # Remove ERCC, rRNA entries; renormalize TPM to 1M total
   cleanTxExprMat<-newTxExprMat[
      grep("ERCC|rRNA", newTxExprMat$gene_name, invert = TRUE),]
   
   print("Renormalizing expression table after removing ERCC, rRNA")
   for (i in (1:length(tpmCols))) {
      cleanTxExprMat[,tpmCols[i]]<-(1E6 * cleanTxExprMat[,tpmCols[i]] / 
                                    sum(cleanTxExprMat[,tpmCols[i]]))
   }
   
   
   # Format to max 3 decimal digits
   if (printTables) {
      x<-cleanTxExprMat
      for (i in c(tpmCols, estCountCols)) {
         x[,i]<-sprintf("%.3f", x[,i])}
   
      gz1<-gzfile(paste0(specs$outDir,
                         "/transcriptExpressionMatrix.no_RRNA_ERCC.txt.gz"),
                  "w")
      write.table(x, file=gz1, col.names=TRUE, row.names=FALSE,
                  quote=FALSE,sep="\t")
      close(gz1)
      x<-NULL
   }
   
   
   print("Aggregating to gene expression matrix")
   geneExprMat<- txExprMat[, c("gene_id", "gene_name",
                               tpmCols, estCountCols)] %>% 
                 group_by(gene_id, gene_name) %>%
                 summarise_all(funs(sum))

   # Format to max 3 decimal digits
   if (printTables) {
      x <- geneExprMat

      for (i in c(tpmCols, estCountCols)) { x[[i]] <- sprintf("%.3f", x[[i]]) }
   
      gz1<-gzfile(paste0(specs$outDir, "/geneExpressionMatrix.txt.gz"), "w")
      write.table(x, file=gz1, col.names=TRUE, row.names=FALSE,
                  quote=FALSE, sep="\t")
      close(gz1)
      x<-NULL
   }
   
   
   
   print("Aggregating clean (no rRNA) gene expression matrix")
   cleanGeneExprMat<- cleanTxExprMat[,c("gene_id","gene_name",
                                        tpmCols,estCountCols)]%>% 
                 group_by(gene_id, gene_name) %>%
                 summarise_all(funs(sum))
   
   # Format to max 3 decimal digits
   if (printTables) {
      x<-cleanGeneExprMat
      for (i in c(tpmCols, estCountCols)) {
         x[[i]]<-sprintf("%.3f", x[[i]])}
   
      gz1<-gzfile(paste0(specs$outDir,
                         "/geneExpressionMatrix.no_RRNA_ERCC.txt.gz"),"w")
      write.table(x, file=gz1, col.names=TRUE,
                  row.names=FALSE, quote=FALSE,sep="\t")
      close(gz1)
      x<-NULL
   }
   
   #head ../../results/RNAseq/kallisto/opticlobe_2015/Lawf1_d1-rep1-2015/abundance.txt
   #target_id       length  eff_length      est_counts      tpm
   #FBtr0082757     1844    1595    9       0.718348

   return(list(
      geneExpr          = geneExprMat,
      txExpr            = txExprMat,
      geneExpr.clean    = cleanGeneExprMat,
      txExpr.clean      = cleanTxExprMat
   ))

}



defineDEG <- function(dat,
                      plotDEG=TRUE ){

# consistent DEG pairwise SLEUTH comparison

   celltypePairs <- dat$specs$deg.celltypePairs

   genes.deg <- list()
   t2g <- data.frame(
      target_id = dat$specs$transcriptInfo$transcript_id,
      ens_gene  = dat$specs$transcriptInfo$gene_id,
      ext_gene  = dat$specs$transcriptInfo$gene_name,
      stringsAsFactors = FALSE
   )

   allSamples <- dat$specs$rnaSamples

   for (cellPair in celltypePairs) {
   
      print(paste0("comparing ",cellPair$cellType1," to ",cellPair$cellType2))
      cell12 <- paste0(cellPair$cellType1,"_vs_",cellPair$cellType2)
      cell1 <- cellPair$cellType1
      cell2 <- cellPair$cellType2


## HERENOW: more flexible sample specs

      samples1   <- allSamples$sampleName[allSamples$cellType == cellPair$cellType1 &
                                          allSamples$dataType2 == cellPair$dataType2]

      samples2   <- allSamples$sampleName[allSamples$cellType == cellPair$cellType2 &
                                          allSamples$dataType2 == cellPair$dataType2]


      print(paste0("samples 1: ", paste0(samples1, collapse = ", ")))
      print(paste0("samples 2: ", paste0(samples2, collapse = ", ")))
   
      curSamples <- c(samples1, samples2)
         
      curSampleInfo <- allSamples[match(curSamples, allSamples$sampleName), ]
      print(curSampleInfo)


# Weirdo sleuth access business; coefficient exists in sleuth_wt() for the 
# condition with alphanumercally lesser name

      cond1 <- paste0(1,".",cellPair$cellType1)
      cond2 <- paste0(2,".",cellPair$cellType2)
         
      curSampleInfo$condition <- c(rep(cond1, length(samples1)),
                                   rep(cond2, length(samples2)))
         
      curSampleInfo$path <- paste0(dat$specs$kallistoBaseDir, "/",
                                      curSampleInfo$runID, "/",
                                      curSampleInfo$sampleID)

# sleuth info https://rawgit.com/pachterlab/sleuth/master/inst/doc/intro.html
   
      comparisonString <- paste0( paste0(sort(samples1), collapse = "_"),
                                  "_vs_",
                                  paste0(sort(samples2), collapse = "_") )
   
      comparisonMD5 <- digest(comparisonString, algo = "md5")
      curSleuthPath <- paste0(dat$specs$sleuthBaseDir, "/", comparisonMD5)
   
      comparisonFC <- data.frame(
         gene_name      = dat$edat$geneExpr$gene_name,
         tpm.celltype2  = apply(dat$edat$geneExpr[, paste0("tpm.", samples2)],
                                1, mean),
         tpm.celltype1  = apply(dat$edat$geneExpr[, paste0("tpm.", samples1)],
                                1, mean),
         stringsAsFactors = FALSE
      )
   
      exprGenes.up <- comparisonFC[
         comparisonFC$tpm.celltype1 >= dat$specs$thresh$exprGenes.minTPM &
         comparisonFC$tpm.celltype1 >= comparisonFC$tpm.celltype2 *
             dat$specs$thresh$deGenes.FC ,
         "gene_name"]
   
      exprGenes.down <- comparisonFC[
         comparisonFC$tpm.celltype2 >= dat$specs$thresh$exprGenes.minTPM &
         comparisonFC$tpm.celltype2 >= comparisonFC$tpm.celltype1 * 
            dat$specs$thresh$deGenes.FC ,
         "gene_name"]
   
      comparisonFC$log2fc.celltype1_vs_celltype2 <- round(
         (1 + comparisonFC$tpm.celltype1) / (1 + comparisonFC$tpm.celltype2),
         digits = 2 )
   
      if (file.exists(paste0(curSleuthPath, "/done.txt"))){
   
         sleuthData <- readSleuthOutput(resultsDir = curSleuthPath)
         genes.deg[[cell12]]$upGenes <- sleuthData$upGenes
         genes.deg[[cell12]]$downGenes <- sleuthData$downGenes
         genes.deg[[cell12]]$sleuthResults <- sleuthData$sleuthResults
         genes.deg[[cell12]]$designMat <- sleuthData$designMat
         if (!is.na(sleuthData$so)) genes.deg[[cell12]]$so <- sleuthData$so
   
      } else {
   
         designMat <- data.frame(
            sample    = curSampleInfo$sampleName,
            condition = curSampleInfo$condition,
            path      = curSampleInfo$path,
            stringsAsFactors = FALSE
         )
         genes.deg$designMat <- designMat
   
         curDesignMat <- designMat[designMat$condition %in% c(cond1, cond2),]
         genes.deg[[cell12]]$designMat <- curDesignMat
         print(curDesignMat)
      
         print(paste("sleuth_prep()", date()))
         genes.deg[[cell12]]$so <- sleuth_prep(curDesignMat, ~ condition,
                                              target_mapping = t2g)
      
         print(paste("sleuth_fit()", date()))
         genes.deg[[cell12]]$so <- sleuth_fit(genes.deg[[cell12]]$so)
      
         print(paste("sleuth_wt()", date()))
         genes.deg[[cell12]]$so <- sleuth_wt(genes.deg[[cell12]]$so,
                                            paste0('condition',cond2))
      
         print(paste("sleuth_results()", date()))
         genes.deg[[cell12]]$sleuthResults <- sleuth_results(
            genes.deg[[cell12]]$so, paste0('condition',cond2))
         print(paste("DONE", date()))
   
         genes.deg[[cell12]]$upGenes <- unique(na.omit(
           genes.deg[[cell12]]$sleuthResults[
               genes.deg[[cell12]]$sleuthResults$qval <=
               dat$specs$thresh$sleuth.qval &
               genes.deg[[cell12]]$sleuthResults$b < 0, ])$ext_gene)
   
         genes.deg[[cell12]]$downGenes <- unique(na.omit(
           genes.deg[[cell12]]$sleuthResults[
               genes.deg[[cell12]]$sleuthResults$qval <=
               dat$specs$thresh$sleuth.qval &
               genes.deg[[cell12]]$sleuthResults$b > 0, ])$ext_gene)
   
         writeSleuthOutput(
            sleuthResults = genes.deg[[cell12]]$sleuthResults,
            upGenes       = genes.deg[[cell12]]$upGenes,
            downGenes     = genes.deg[[cell12]]$downGenes,
            designMat     = genes.deg[[cell12]]$designMat,
            outDir        = curSleuthPath,
            comparisonName= comparisonString)
     }
   
      genes.deg[[cell12]]$fc <- comparisonFC[ comparisonFC$gene_name %in%
         unique(c(genes.deg[[cell12]]$upGenes,
                  genes.deg[[cell12]]$downGenes)), ]
   
      genes.deg[[cell12]]$upGenes.minExpr <- sort(intersect(
         genes.deg[[cell12]]$upGenes,
         exprGenes.up))
   
      genes.deg[[cell12]]$downGenes.minExpr <- sort(intersect(
         genes.deg[[cell12]]$downGenes,
         exprGenes.down))

      if (plotDEG) {
         outFn <- paste0(dat$paths$outDir,"/bulkDEGscatter.",cell12,".pdf")
         pdf(outFn)
         plot(1 + comparisonFC$tpm.celltype2,
              1 + comparisonFC$tpm.celltype1,
              pch=20,
              cex=0.5,
              log="xy",
              xlab = paste0(cell2," (TPM)"),
              ylab = paste0(cell1," (TPM)"),
              main=cell12)
         degGenes.up <- comparisonFC[comparisonFC$gene_name %in%
                                     genes.deg[[cell12]]$upGenes.minExpr,]
         degGenes.down  <- comparisonFC[comparisonFC$gene_name %in%
                                        genes.deg[[cell12]]$downGenes.minExpr,]
                                 
         points(1 + degGenes.up$tpm.celltype2,
                1 + degGenes.up$tpm.celltype1,
                pch=20,cex=0.5,
                col="#f1a340")

         points(1 + degGenes.down$tpm.celltype2,
                1 + degGenes.down$tpm.celltype1,
                pch=20,cex=0.5,
                col= "#998ec3")
         dev.off()
      }
   
   }

   if ("deg" %in% names(dat)) {
      for (cell12 in names(genes.deg)) {
         if (!cell12%in% names(dat$deg)) {
            dat$deg[[cell12]] <- genes.deg[[cell12]] } }
   } else {
      dat$deg <- genes.deg
   }

   return(dat)

}



setSpecs <- function(){
# CR: note the single cell routein uses inferstate's scoreSCsig.setPaths

#   library(readxl)

   specs<-list(
      baseDir           = c("~/data/projects/exfoxp3")
   )

   specs$thresh<-list()
   specs$thresh$exprGenes.minTPM <- 10
   specs$thresh$deGenes.FC <- 1.5
   specs$thresh$sleuth.qval <- 0.05
   specs$thresh$maxGenePeakDist <- 50000

   specs$thresh$GOstats.pvalue <- 0.001

   specs$outDir <- paste0(specs$baseDir, "/analysis/20181004.paperFigs")

   specs$outExprFn <- paste0(specs$outDir, "/geneExpr.txt.gz")


   specs$geneStrxDir <- paste0(
      specs$baseDir,"/data/gene_annotation.GRCm38.ENSEMBL82/")
   specs$geneStrxTSS      <- paste0(specs$geneStrxDir,"/mm10_tss_exact.bed")
   specs$geneStrxGeneBody <- paste0(specs$geneStrxDir, "/mm10_genebodies.bed")
   specs$transcriptInfo <- read.table(paste0(specs$baseDir,
      "/data/txInfo/GRCm38.82.withpatch.ERCC.transcript_info.txt"),
      quote = "", header = TRUE, sep = "\t", as.is = TRUE)
   specs$geneInfo <- unique(specs$transcriptInfo[, c("gene_id", "gene_name")])


   specs$kallistoBaseDir <- paste0(specs$baseDir,
      "/results/RNAseq/kallisto/")

   specs$kallistoBaseDir.stranded <- paste0(specs$baseDir,
      "/results/RNAseq/kallisto_stranded/")

   specs$sleuthBaseDir <- paste0(specs$baseDir, "/results/RNAseq/sleuth/")

   specs$bin$bedtools <- "/usr/local/apps/bedtools/2.25.0/bin/bedtools"
   specs$bin$samtools <- "/usr/local/apps/samtools/1.2/bin/samtools"
   specs$bin$ucsctools <- "/usr/local/apps/ucsc/365/bin/x86_64"

   print("LOADING RNA SAMPLES")
   specs$rnaSamples <- read.table(paste0(specs$baseDir,
      "/metadata/exfoxp3_rna_samples.txt"), header=TRUE, sep="\t", as.is=TRUE)
   print(" DONE!")

#   specs$rnaSamples$cellType <- gsub(".si[0-9]+$", "", specs$rnaSamples$sampleName)


# Remove outliers
   specs$rnaSamples.outliers <- c()
   if (length(specs$rnaSamples.outliers) > 0) {
      print(paste("Skipping RNA outliers: ",
                  paste0(specs$rnaSamples.outliers, collapse=", ")))

      specs$rnaSamples <- specs$rnaSamples[!(specs$rnaSamples$sampleName %in%
                                           specs$rnaSamples.outliers),]
   }


# Remove chromium entries
   specs$rnaSamples <- specs$rnaSamples[specs$rnaSamples$libraryType != "chromium",]

   specs$chipSamples <- read.table(paste0(specs$baseDir,
      "/metadata/exfoxp3_chip_samples.txt"), header=TRUE, sep="\t", as.is=TRUE)

   specs$transcriptInfo <- read.table(paste(specs$baseDir,
      "/data/txInfo/GRCm38.82.FPtags.ERCC.transcript_info.txt",sep=""),
      quote="", header=TRUE,sep="\t",as.is=TRUE)


# Chromium 10X setup
   specs$sc <- list()
   specs$sc$dataDir <- paste0(specs$baseDir,"/run/",
                              "20171120.cellranger_aggr/trth2_sc201711")
   specs$sc$fn_tx_gene_map <- paste0(specs$baseDir,
                                     "/data/kallisto_files/tx_gene_name_map.txt")

   specs$sc$txGeneMap <- read.table(specs$sc$fn_tx_gene_map,header=FALSE,sep=" ")
   colnames(specs$sc$txGeneMap) <- c("gene_id", "tx_id", "gene_name")

   return(specs)

}


plotHmap.bulkReporter <- function(dat,
                                  nlMode = "minMax",
                                  geneList,
                                  show_rownames = TRUE,
                                  clusterCols = FALSE,
                                  clusterRows = TRUE,
                                  rowFontSize = 2,
                                  curHeight = 3,
                                  curWidth = 1.5,
                                  rowGaps,
                                  repAvg = FALSE,
                                  figName = "bulkReporterRNA_lung_hmap",
                                  dataset = "reporter", #nascent
                                  decontamEx = FALSE,
                                  decontamEx.f_c = 0.3
                                 ) {

# PURPOSE: make heatmap of bulk reporter RNA-seq data from lung

   curFig <- figName

   if (dataset == "reporter") {
      curCellTypes <- c("lung.Foxp3p", "lung.exFoxp3", "lung.nonFoxp3")
      labelCellTypes <- c("Foxp3p", "ex-Foxp3", "non-Foxp3")
      dataType2 <- "bulk"

      degTypes <- c("lung.nonFoxp3_vs_lung.exFoxp3",
                    "lung.Foxp3p_vs_lung.exFoxp3",
                    "lung.nonFoxp3_vs_lung.Foxp3p")

   } else if (dataset == "nascent") {

      curCellTypes <- c("iTreg", "TrTh2", "Th2")
      labelCellTypes <- c("iTreg", "eF-Treg", "Th2")
      dataType2 <- "nascent"

      degTypes <-  c( "Th2_vs_TrTh2", "iTreg_vs_TrTh2")
   }


   curSampleList <- dat$specs$rnaSamples[dat$specs$rnaSamples$dataType2 == dataType2 &
                                      dat$specs$rnaSamples$cellType %in% curCellTypes,]

   curSamples <- curSampleList$sampleName
   tpmCols <- paste0("tpm.", curSamples)


   if (missing(geneList)) {
      geneList <- c()
      for (comparison in degTypes) {
         geneList <- c(geneList, dat$deg[[comparison]]$upGenes.minExpr,
                                 dat$deg[[comparison]]$downGenes.minExpr)
      }
      geneList <- unique(geneList)
   }



   origExpr <- dat$edat$geneExpr[,c("gene_name","gene_id",tpmCols)]

   hMat <- as.data.frame(origExpr[
      match(geneList[geneList %in% origExpr$gene_name], origExpr$gene_name),
      c("gene_name", "gene_id", tpmCols)])
   rownames(hMat) <- paste0(hMat$gene_name, ":", hMat$gene_id)
   rowLabels <- hMat$gene_name
   hMat$gene_name <- NULL
   hMat$gene_id <- NULL


   if (repAvg) {
      newTpmCols <- c()
      for (celltype in curCellTypes) {
         curTypeSamples <- dat$specs$rnaSamples$sampleName[
            dat$specs$rnaSamples$dataType2 == dataType2 &
            dat$specs$rnaSamples$cellType == celltype]

         newTpmCol <- paste0("meantpm.",celltype)
         newTpmCols <- c(newTpmCols, newTpmCol)

         hMat[,newTpmCol] <- rowMeans(hMat[,paste0("tpm.",curTypeSamples)])
      }
      hMat <- hMat[,newTpmCols]

      if (dataset == "reporter" & decontamEx) {
         newTpmCol <- "meantpm.lung.exFoxp3.dec"
         hMat[,newTpmCol] <- ( hMat[,"meantpm.lung.exFoxp3"] -
                               decontamEx.f_c * hMat[,"meantpm.lung.Foxp3p"]) /
                             (1 - decontamEx.f_c)
         hMat[hMat[,newTpmCol] < 0, newTpmCol] <- 0

         curCellTypes <- c("lung.Foxp3p", "lung.exFoxp3", "lung.exFoxp3.dec", "lung.nonFoxp3")
         labelCellTypes <- c("Foxp3p", "ex-Foxp3", "ex-Foxp3-dec", "non-Foxp3")

         hMat[,"meantpm.lung.exFoxp3"] <- hMat[,"meantpm.lung.exFoxp3.dec"]
         curCellTypes <- c("lung.Foxp3p", "lung.exFoxp3", "lung.nonFoxp3")
         labelCellTypes <- c("Foxp3p", "ex-Foxp3", "non-Foxp3")
         hMat <- hMat[,paste0("meantpm.", curCellTypes)]
      }
   }

   hMat <- as.matrix(hMat)

   hMat <- (1 + hMat)


   if (nlMode == "logMean") {
      bottomCap <- -1.5
      topCap <- 1.5

      hMat <- log2(hMat / apply(hMat, 1, mean))
      hMat[hMat < bottomCap] <- bottomCap; hMat[hMat > topCap] <- topCap;
      curBreaks <- seq(bottomCap, topCap, length.out = 101)
      plotTitle  <- "Rel. expression: log2 (TPM/mean)"
      plotColors <- colorRampPalette(c("steelBlue2", "white", "darkOrange2"))(100)
   } else if (nlMode == "minMax") {
      hMin <-  apply(hMat, 1, min)
      hMax <-  apply(hMat, 1, max)
      hMat <- (hMat - hMin) / (hMax - hMin)
      plotTitle  <- "(TPM - min)/(max - min))"
      curBreaks <- seq(0, 1, length.out = 101)
      plotColors <- colorRampPalette(c("steelBlue2", "white", "darkOrange2"))(100)
   } else if (nlMode == "fracMax") {
      hMax <-  apply(hMat, 1, max)
      hMat <- hMat / hMax
      plotTitle  <- "TPM / max"
      curBreaks <- seq(0, 1, length.out = 101)
      plotColors <- colorRampPalette(c("white", "red"))(100)
   }


   
   colAnn <- data.frame(
                      celltype = rep("X", ncol(hMat)),
                      stringsAsFactors = FALSE)
   rownames(colAnn) <- colnames(hMat)
   print(paste0("COLUMN NAMES FOR ANNITATION"))
   print(rownames(colAnn))

   if (dataset == "reporter") {

      if (! repAvg) {
         colAnn$celltype[grep("RFPp_GFPn", rownames(colAnn))] <- c("ex-Foxp3")
         colAnn$celltype[grep("RFPp_GFPp", rownames(colAnn))] <- c("Foxp3")
         colAnn$celltype[grep("RFPn_GFPn", rownames(colAnn))] <- c("non-Foxp3")
      } else if (repAvg) {
         colAnn$celltype[grep("lung.Foxp3p", rownames(colAnn))] <- c("Foxp3")
         colAnn$celltype[grep("lung.exFoxp3", rownames(colAnn))] <- c("ex-Foxp3")
         colAnn$celltype[grep("lung.nonFoxp3", rownames(colAnn))] <- c("non-Foxp3")
      } 

   } else {

      colAnn$celltype[grep("TrTh2", rownames(colAnn))] <- c("TrTh2")
      colAnn$celltype[grep("iTreg", rownames(colAnn))] <- c("iTreg")
      colAnn$celltype[grep("Th2", rownames(colAnn))] <- c("Th2")

   }

   if (0) {
   if (nrow(hMat) > 500) {
      rowFontSize <- 1.5
      curHeight <- 20
      curWidth <- 5
   } else if (nrow(hMat) > 70)  {
      rowFontSize <- 3
      curHeight <- 7
      curWidth <- 5
   } else if (nrow(hMat) > 50)  {
      rowFontSize <- 3.5
      curHeight <- 7
      curWidth <- 5
   } else if (nrow(hMat) > 20)  {
      rowFontSize <- 5
      curHeight <- 7
      curWidth <- 5
   } else {
      rowFontSize <- 7
      curHeight <- 7
      curWidth <- 5
   }
   }

#   show_rownames <- TRUE
#   curHeight <- 5
#   curWidth <- 7
#   show_rownames <- FALSE


   gaps_row <- c()
   if (!missing(rowGaps)) {
      gaps_row <- rowGaps
   }

   hMat <- na.omit(hMat)

   pheatmap.options <- list(hMat,
            show_rownames  = show_rownames,
            show_colnames  = FALSE,
            cluster_cols   = clusterCols,
            cluster_rows   = clusterRows,
            labels_row = rowLabels,
            treeheight_row = 0,
            breaks         = curBreaks,
            gaps_row       = gaps_row,
            color          = colorRampPalette(c("steelBlue2", "white",
                                                "darkOrange2"))(100),
            border_color = NA,
            fontsize       = 9,
            fontsize_row   = rowFontSize,
            main = "")
#            main           = plotTitle)

   if (!repAvg) {
      pheatmap.options$annotation_col <- colAnn
      pheatmap.options$annotation_names_row <- TRUE

      if (dataset == "reporter") {
         pheatmap.options$annotation_colors <- list(
            celltype = c("ex-Foxp3"  = "darkOrange2",
                        "Foxp3" = "steelBlue2",
                        "non-Foxp3"     = "gray80"))
      } else if (dataset == "nascent") {
         pheatmap.options$annotation_colors <- list(
            celltype = c("TrTh2-Foxp3"  = "darkOrange2",
                        "iTreg" = "steelBlue2",
                        "Th2"     = "gray80"))
      }
      pheatmap.options$annotation_legend <- TRUE
      pheatmap.options$annotation_names_col <- TRUE
   } else {
      pheatmap.options$show_colnames <- TRUE
      pheatmap.options$labels_col <- labelCellTypes
   }

   print(" BOUT TO TSART PHEATMAP PDF")
   outFn <- paste0(dat$specs$outDir, "/msFig", curFig, "_",nlMode,".pdf")
   pdf(outFn,
       onefile = FALSE,
       family = "ArialMT",
       height = curHeight,
       width  = curWidth)
   do.call(pheatmap, pheatmap.options)
   dev.off()
   print(" ---> GOT BACK")

   curSampleList <- unique(curSampleList)
   figDescrFn <- paste0(dat$specs$outDir, "/", curFig, "_samples.txt")
   write.table(curSampleList, figDescrFn,
               row.names = FALSE, col.names = FALSE, quote = FALSE)

}



makeGeoExpressionTable <- function(dat, tabName = "geoTab1") {

   outSamples <- dat$specs$rnaSamples
   sampleOrder <- order(outSamples$geoName)

   oldNames <- outSamples$sampleName[sampleOrder]
   newNames <- outSamples$geoName[sampleOrder]

   oldCols <-paste0("tpm.",oldNames)
   newCols <-paste0("tpm.",newNames)

   outFn <- paste0(dat$specs$outDir, "/", tabName, "_bulk_expression_table.tsv.gz")
   outDat <- dat$edat$geneExpr.clean[,c("gene_id","gene_name",oldCols)]
   colnames(outDat) <- c("gene_id","gene_name",newCols)

   outDat[,newCols] <- round(outDat[,newCols],digits=2)

   gz1<-gzfile(outFn,"w")

   write.table(outDat,file=gz1, row.names=FALSE,col.names=TRUE,quote=FALSE,sep="\t")
   close(gz1)

}


scoreSCsig.setPaths <- function(dat,
                                baseDir,
                                kallistoDir,
                                sleuthDir,
                                rnaSamplesFn,
                                scDataDir,
                                scDataFn,
                                outDir       = "./") {

   if (missing(dat)) dat <- list()

#   library(readxl)

   paths<-list()
   paths$baseDir <- baseDir
   paths$outDir <- outDir

   if (! missing(kallistoDir)) {
      paths$kallistoBaseDir <- paste0(paths$baseDir,"/",kallistoDir)
      paths$sleuthBaseDir <- paste0(paths$baseDir, "/",sleuthDir)
   }

   if (! missing(rnaSamplesFn)) {
      paths$rnaSamplesFn <- paste0(paths$baseDir, "/",rnaSamplesFn)
   }


# Chromium 10X setup
   if (!missing(scDataFn)) {
      paths$sc.dataFn <- paste0(paths$baseDir,"/", scDataFn)
   } else {
      paths$sc.dataDir <- paste0(paths$baseDir,"/", scDataDir)
   }

   paths$transcriptInfoFn <- paste0(paths$baseDir,
         "/data/txInfo/GRCm38.82.FPtags.ERCC.transcript_info.txt")
   paths$sc.fn_tx_gene_map <- paste0(specs$baseDir,
                                     "/data/kallisto_files/tx_gene_name_map.txt")
   paths$sc.geneInfoFn  <- paste0(paths$baseDir,"/data/cellranger_files.v2.0.0/txInfo/genes.tsv")

   dat$paths <- paths
   return(dat)

}


scoreSCsig.setSpecs <- function(dat,
                                deg.celltypePairs,
                                deg.geneLists,
                                sc.clusterType,
                                manualLibNames,
                                experimentName,
                                keep.sc.samples,
                                skip.sc.samples){

#   library(readxl)

   specs <- list()
   specs$thresh<-list()
   specs$thresh$exprGenes.minTPM <- 10
   specs$thresh$deGenes.FC <- 1.5
   specs$thresh$sleuth.qval <- 0.05
   specs$thresh$maxGenePeakDist <- 50000

   specs$thresh$GOstats.pvalue <- 0.001

   specs$transcriptInfo <- read.table(dat$paths$transcriptInfoFn,
      quote = "", header = TRUE, sep = "\t", as.is = TRUE)
   specs$geneInfo <- unique(specs$transcriptInfo[, c("gene_id", "gene_name")])


   if ("rnaSamplesFn" %in% names(dat$paths)) {
      print("LOADING RNA SAMPLES")
      specs$rnaSamples <- read.table( dat$paths$rnaSamplesFn,
                                   header=TRUE, sep="\t", as.is=TRUE )
      print(" DONE!")
# Remove chromium entries
      specs$rnaSamples <- specs$rnaSamples[specs$rnaSamples$libraryType != "chromium",]
   }

   if (!missing(deg.celltypePairs)) { specs$deg.celltypePairs <- deg.celltypePairs }
   if (!missing(deg.geneLists)) { specs$deg.geneLists <- deg.geneLists }

#   specs$rnaSamples$cellType <- gsub(".si[0-9]+$", "", specs$rnaSamples$sampleName)



   if (!missing(manualLibNames)) { #to assign names by X to BARCODE-X cells
      specs$manualLibNames <- manualLibNames }

   if (!missing(keep.sc.samples)) {
      specs$keep.sc.samples <- keep.sc.samples }

   if (!missing(skip.sc.samples)) {
      specs$skip.sc.samples <- skip.sc.samples }

   if (!missing(sc.clusterType)) {
      specs$sc.clusterType <- sc.clusterType }

   if (!missing(experimentName)) {
      specs$experimentName <- experimentName
   } else {
      specs$experimentName <- "alldata"
   }


# Chromium 10X setup

   specs$sc <- list()
   if ("sc.dataDir" %in% names(dat$paths)) {
      specs$sc$dataType <- "cellranger"
   } else {
      specs$sc$dataType <- "seurat"
   }

   if (file.exists(dat$paths$sc.fn_tx_gene_map)) {
   specs$sc$txGeneMap <- read.table(dat$paths$sc.fn_tx_gene_map,
                                    header=FALSE,sep=" ")

   colnames(specs$sc$txGeneMap) <- c("gene_id", "tx_id", "gene_name")
   }

   dat$specs <- specs
   return(dat)

}



readSleuthOutput <- function(comparisonName, resultsDir){

   print(paste0("Reading SLEUTH results from ",resultsDir))

   comparisonName<-scan(file=paste0(resultsDir,"/comparison.txt"),
                        what="character")

   so<-NA
   soFn<-paste0(resultsDir,"/so.rds")
   if (file.exists(soFn)) { so<-readRDS(file=soFn) }

   sleuthResults<-readRDS(file=paste0(resultsDir,"/sleuthResults.rds"))
   designMat<-readRDS(file=paste0(resultsDir,"/designMat.rds"))
   upGenes<-scan(file=paste0(resultsDir,"/upGenes.txt"),what="character")
   downGenes<-scan(file=paste0(resultsDir,"/downGenes.txt"),what="character")

   return(list(
      comparisonName= comparisonName,
      so            = so,
      designMat     = designMat,
      sleuthResults = sleuthResults,
      upGenes       = upGenes,
      downGenes     = downGenes
   ))

}



writeSleuthOutput <- function(comparisonName,
                              so=NULL,
                              sleuthResults,
                              designMat,
                              upGenes=NULL,
                              downGenes=NULL,
                              outDir){
# Purpose: Given SLEUTH return object, write DE list to file
#          ALT: save R object dump
#             name with MD5 hash of string of samples compared

   if (!file.exists(outDir)){
      dir.create(outDir,recursive=TRUE) }

   print(paste0("Saving SLEUTH results to ",outDir))

# Write comparison string to disk
   fileConn<-file(paste0(outDir,"/comparison.txt"))
   writeLines(comparisonName, fileConn)
   close(fileConn)

# Dump Sleuth objects to disk
   if (!is.null(so)){ saveRDS(so, file=paste0(outDir,"/so.rds")) }
   saveRDS(sleuthResults, file=paste0(outDir,"/sleuthResults.rds"))
   saveRDS(designMat, file=paste0(outDir,"/designMat.rds"))

# Write gene list to disk
   if (!(is.null(upGenes))){
      upFn<-paste0(outDir,"/upGenes.txt")
      write.table(upGenes, file=upFn, quote=FALSE, row.names=FALSE)
   }

   if (!(is.null(downGenes))){
      downFn<-paste0(outDir,"/downGenes.txt")
      write.table(downGenes, file=downFn,quote=FALSE, row.names=FALSE)
   }

# Put file down saying run done
   fileConn<-file(paste0(outDir,"/done.txt"))
   writeLines("1", fileConn)
   close(fileConn)

   return(1)
}

mapColor <- function(x,
                     colorPal = colorRampPalette(c("steelBlue3","white","darkOrange3"))(100),
                     minCap,
                     maxCap,
                     minScore,
                     maxScore,
                     alpha,
                     symmetric=TRUE){


   if (!missing(minCap)) { x[x < minCap] <- minCap; }
   if (!missing(maxCap)) { x[x > maxCap] <- maxCap; }

   if (missing(minScore) | missing(maxScore)) {
   if (!symmetric) {
      minScore <- min(x)
      maxScore <- max(x)
   } else if (symmetric) {
      maxScore <- max(abs(x))
      minScore <- maxScore * -1
   }
   }

   print(paste0("minScore=",minScore))
   print(paste0("maxScore=",maxScore))
   cols <- colorPal[findInterval(x,seq(minScore,maxScore,
                                       length.out=length(colorPal)+1),
                                 all.inside=TRUE)]

   if (!missing(alpha)) {cols <- addAlpha(cols,alpha=alpha)}
   return(cols)
}
