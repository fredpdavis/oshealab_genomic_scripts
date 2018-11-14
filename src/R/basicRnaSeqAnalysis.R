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
      print("Loading differentially expressed genes")
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
       plotHmap.deg(dat$dat.bulk,
                    figName    = "2D",
                    dataset    = "reporter",
                    nlMode     = "fracMax",
                    repAvg        = TRUE,
                    clusterRows   = TRUE,
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

#ORIG:   celltypePairs <- dat$specs$deg.celltypePairs

   genes.deg <- list()
   t2g <- data.frame(
      target_id = dat$specs$transcriptInfo$transcript_id,
      ens_gene  = dat$specs$transcriptInfo$gene_id,
      ext_gene  = dat$specs$transcriptInfo$gene_name,
      stringsAsFactors = FALSE
   )

   allSamples <- dat$specs$rnaSamples

#ORIG   for (cellPair in celltypePairs)
   for (i in 1:nrow(dat$specs$rnaComps)) {
      cellPair <- 

# HERENOW 181114_1242
      sampleGroup1 <- dat$specs$rnaComps$sampleGroup1[i]
   
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

   specs<-list(
      baseDir           = c("~/data/projects/cytokineX")
   )

   specs$thresh<-list()
   specs$thresh$exprGenes.minTPM <- 10
   specs$thresh$deGenes.FC <- 1.5
   specs$thresh$sleuth.qval <- 0.05
   specs$thresh$maxGenePeakDist <- 50000

   specs$thresh$GOstats.pvalue <- 0.001

   specs$outDir <- paste0(specs$baseDir, "/analysis/basicFigures")

   specs$outExprFn <- paste0(specs$outDir, "/geneExpr.txt.gz")

   specs$transcriptInfo <- read.table(paste0(specs$baseDir,
      "/data/txInfo/GRCm38.82.withpatch.ERCC.transcript_info.txt"),
      quote = "", header = TRUE, sep = "\t", as.is = TRUE)
   specs$geneInfo <- unique(specs$transcriptInfo[, c("gene_id", "gene_name")])


   specs$kallistoBaseDir <- paste0(specs$baseDir,
      "/results/RNAseq/kallisto.GRCm38.94/")

   specs$sleuthBaseDir <- paste0(specs$baseDir, "/results/RNAseq/sleuth/")


   print("Loading RNA-seq sample information")
   specs$rnaSamples <- read.table(paste0(specs$baseDir,
      "/metadata/rnaseq_samples.txt"), header=TRUE, sep="\t", as.is=TRUE)
   print(" DONE!")

   print("Loading pairs of RNA-seq samples to compare")
   specs$rnaComps <- read.table(paste0(specs$baseDir,
      "/metadata/rnaseq_comparisons.txt"), header=TRUE, sep="\t", as.is=TRUE)
   print(" DONE!")

#   specs$rnaSamples$cellType <- gsub(".si[0-9]+$", "", specs$rnaSamples$sampleName)

# Remove outliers
   if (sum(specs$rnaSamples$outlier == "yes")) {
      print(paste("Skipping RNA outliers: ",
         paste0(specs$rnaSamples$sampleName[specs$rnaSamples$outlier == "yes"],
                collapse=", ")))

      specs$rnaSamples <- specs$rnaSamples[specs$rnaSamples.outliers != "yes",]
   }


   specs$transcriptInfo <- read.table(paste(specs$baseDir,
      "/data/kallisto_files.GRCm38.94/transcript_info.GRCm38.94.txt",sep=""),
      quote="", header=TRUE,sep="\t",as.is=TRUE)

   return(specs)

}


plotHmap.deg <- function(dat,
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
                         figName = "deg_hmap"
                         ) {

# PURPOSE: make heatmap of bulk reporter RNA-seq data from lung

   curFig <- figName

   curCellTypes <- c("iTreg", "TrTh2", "Th2")
   labelCellTypes <- c("iTreg", "eF-Treg", "Th2")
   dataType2 <- "nascent"
   degTypes <-  c( "Th2_vs_TrTh2", "iTreg_vs_TrTh2")


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

   colAnn$celltype[grep("TrTh2", rownames(colAnn))] <- c("TrTh2")
   colAnn$celltype[grep("iTreg", rownames(colAnn))] <- c("iTreg")
   colAnn$celltype[grep("Th2", rownames(colAnn))] <- c("Th2")

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

   if (!repAvg) {
      pheatmap.options$annotation_col <- colAnn
      pheatmap.options$annotation_names_row <- TRUE

      pheatmap.options$annotation_colors <- list(
            celltype = c("TrTh2-Foxp3"  = "darkOrange2",
                         "iTreg" = "steelBlue2",
                         "Th2"     = "gray80"))
      pheatmap.options$annotation_legend <- TRUE
      pheatmap.options$annotation_names_col <- TRUE
   } else {
      pheatmap.options$show_colnames <- TRUE
      pheatmap.options$labels_col <- labelCellTypes
   }

   outFn <- paste0(dat$specs$outDir, "/msFig", curFig, "_",nlMode,".pdf")
   pdf(outFn,
       onefile = FALSE,
       family = "ArialMT",
       height = curHeight,
       width  = curWidth)
   do.call(pheatmap, pheatmap.options)
   dev.off()

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
