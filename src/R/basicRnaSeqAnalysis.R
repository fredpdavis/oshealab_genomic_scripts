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




main <- function(dat,
                 makeFigs = TRUE,
                 returnData = TRUE, ...) {

   startup()

# 0. Load Data

   if (missing(dat))                    dat <- list()

   if (! "specs" %in% names(dat))       dat$specs <- setSpecs(...)

   if (! "edat" %in% names(dat)) {
      print("Loading expression")
      dat$edat <- loadExpr(specs       = dat$specs,
                           printTables = TRUE)
      if (returnData)    return(dat)
   }

   if (! "deg" %in% names(dat)) {
      print("Identifying/Loading differentially expressed genes")
      dat <- defineDEG(dat)
      if (returnData)    return(dat)
   }

   if (makeFigs) {
      print("Making figures")
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

# Heatmap of DEG genes
   if (any(c("all","1") %in% figList)) {
      print("Plot: DEG heatmap")
      plotHmap.deg(dat,
                    figName    = "1",
                    show_rownames = FALSE)
   }

# Heatmap of pairwise sample correlation
   if (any(c("all","2") %in% figList)) {
      print("Plot: Sample correlation heatmap")
      makeCorrHeatmap(dat)
   }

# Heatmap of variable genes
   if (any(c("all","3") %in% figList)) {
      print("Plot: Variable gene heatmap")
      makeVarGeneHeatmap(dat)
   }

# DEG scatterplots
   if (any(c("all","3") %in% figList)) {
      print("Plot: DEG scatterplot and fold change cumulative distribution plots")
      plotDEGscatters(dat)
   }

# GEO supp expression table
   if (any(c("all","GEO_table") %in% figList)) {
      print("Table: Gene expresion table")
      makeGeoExpressionTable(dat)
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

      curAbund<-read.table(
         paste(kallistoDir,"/",
               specs$rnaSamples$sampleID[i],"/",
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
                      plotDEG=TRUE) {

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

   for (i in 1:nrow(dat$specs$rnaComps)) {
      cell1 <- dat$specs$rnaComps[i,"group1name"]
      cell2 <- dat$specs$rnaComps[i,"group2name"]
      cell12 <- paste0(cell1,"_vs_",cell2)
      print(paste0("comparing ",cell1," to ",cell2))

      groupSamples <- dat$specs$degSamples[[cell12]]
      samples1   <- groupSamples[[1]]
      samples2   <- groupSamples[[2]]

      print(paste0("samples 1: ", paste0(samples1, collapse = ", ")))
      print(paste0("samples 2: ", paste0(samples2, collapse = ", ")))
   
      curSamples <- c(samples1, samples2)
         
      curSampleInfo <- allSamples[match(curSamples, allSamples$sampleName), ]
      print(curSampleInfo)


# Weirdo sleuth access business; coefficient exists in sleuth_wt() for the 
# condition with alphanumercally lesser name

      cond1 <- paste0(1,".",cell1)
      cond2 <- paste0(2,".",cell2)
         
      curSampleInfo$condition <- c(rep(cond1, length(samples1)),
                                   rep(cond2, length(samples2)))
         
      curSampleInfo$path <- paste0(dat$specs$kallistoBaseDir, "/",
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



setSpecs <- function(baseDir){

   if (missing(baseDir)) { baseDir <- "~/data/projects/cytokineX" }

   specs<-list(baseDir = baseDir)
   print(paste0("basedir = ",baseDir))

   specs$thresh<-list()
   specs$thresh$exprGenes.minTPM <- 10
   specs$thresh$deGenes.FC <- 1.5
   specs$thresh$sleuth.qval <- 0.05
   specs$thresh$maxGenePeakDist <- 50000

   specs$thresh$GOstats.pvalue <- 0.001

   specs$outDir <- paste0(specs$baseDir, "/analysis/basicFigures")
   if (!file.exists(specs$outDir)){
      dir.create(specs$outDir,recursive=TRUE) }

   specs$outExprFn <- paste0(specs$outDir, "/geneExpr.txt.gz")

   specs$transcriptInfo <- read.table(paste0(specs$baseDir,
      "/data/kallisto_files.GRCm38.94/transcript_info.GRCm38.94.txt"),
      quote = "", header = TRUE, sep = "\t", as.is = TRUE)
   specs$geneInfo <- unique(specs$transcriptInfo[, c("gene_id", "gene_name")])


   specs$kallistoBaseDir <- paste0(specs$baseDir,
      "/results/RNAseq/kallisto.GRCm38.94/")

   specs$sleuthBaseDir <- paste0(specs$baseDir, "/results/RNAseq/sleuth/")


   print("Loading RNA-seq sample information")
   print(paste0("LOOKING IN: ",specs$baseDir,"/metadata/rnaseq_samples.txt"))
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


# Figure out what samples to compare for DEG

   specs$degSamples <- list()
   for (i in 1:nrow(specs$rnaComps)) {

      cell1 <- specs$rnaComps[i,"group1name"]
      cell2 <- specs$rnaComps[i,"group2name"]
      cell12 <- paste0(cell1,"_vs_",cell2)
      print(paste0("comparing ",cell1," to ",cell2))

      groupSamples <- list()
      for (j in c(1,2)) {
         curCriteria  <- specs$rnaComps[i,paste0("group",j,"criteria")]

         criteriaBits <- unlist(strsplit(curCriteria,split=";"))
         orSamples <- c()
         for (k in 1:length(criteriaBits)) {
            curBits <- unlist(strsplit(criteriaBits[k],split=","))

            andSamples <- c()
            for (l in 1:length(curBits)) {
               featValPair <- unlist(strsplit(curBits[l], split="="))
               curSamples <- specs$rnaSamples$sampleName[
                              specs$rnaSamples[,featValPair[1]] == featValPair[2]]
               if (l == 1) {
                  andSamples <- curSamples
               } else {
                  andSamples <- intersect(andSamples, curSamples)
               }
            }
            orSamples <- union(orSamples, andSamples)
         }
         groupSamples[[j]] <- orSamples
         print(paste0("Critera: ",curCriteria))
         print(paste0("Matching samples: ", paste0(orSamples,collapse=", ")))
      }
      specs$degSamples[[cell12]] <- groupSamples

   }

   return(specs)

}


makeCorrHeatmap <- function(dat,
                            sampleAnnotate, #properties to show as sample annotation
                            figName = "corrHeatmap") {

   eMat <- dat$edat$geneExpr[,paste0("tpm.",dat$specs$rnaSamples$sampleName)]

   eMat <- eMat[apply(eMat,1,max) >= dat$specs$thresh$exprGenes.minTPM,]

   eMat <- log2(1 + eMat)
   corMat <- cor(eMat)
   colnames(corMat) <- gsub("tpm.","",colnames(corMat))
   rownames(corMat) <- gsub("tpm.","",rownames(corMat))

   if (!missing(sampleAnnotate)) {
      colAnn <- data.frame( sampleNames = colnames(corMat),
                            stringsAsFactors = FALSE)
      rownames(colAnn) <- colnames(corMat)
      print(paste0("ROWNAMES COLANN=",rownames(colAnn)))
   
      if (!missing(sampleAnnotate)) {
         for (curProp in sampleAnnotate) {
            print(paste0("NOW ON ",curProp))
            colAnn[,curProp] <- dat$specs$rnaSamples[
               match(rownames(colAnn),
                     dat$specs$rnaSamples$sampleName),
               curProp]
         }
      }
      colAnn$sampleNames <- NULL
   }
   print("TEST: SampleAnnotates")
   print(colAnn)

   pheatmap.options <- list(corMat, fontsize_row = 4, fontsize_col = 4,
            main = "Pearson correlation of log2(TPM + 1)",
            border_color = NA,
            show_rownames = TRUE,
            show_colnames = TRUE,
            color = colorRampPalette(c("white", "#008500"))(50)
            )

   if (!missing(sampleAnnotate)) {
      pheatmap.options$annotation_col <-  colAnn
      pheatmap.options$annotation_row <- colAnn
      pheatmap.options$annotation_legend <- TRUE
      pheatmap.options$annotation_names_col <- TRUE
#      pheatmap.options$show_rownames <- FALSE,
#      pheatmap.options$show_colnames <- FALSE
   }

   pdf(paste0(dat$specs$outDir, "/",
           figName, "_correlation_heatmap.pdf"),
       onefile = FALSE,
       height = 5,
       width  = 6)
   do.call(pheatmap, pheatmap.options)
   dev.off()
}


makeVarGeneHeatmap <- function(dat,
                            sampleAnnotate, #properties to show as sample annotation
                            nlMode = "logMean",
                            figName = "variableGeneHeatmap") {

   eMat <- dat$edat$geneExpr[,paste0("tpm.",dat$specs$rnaSamples$sampleName)]
   eMat <- data.frame(eMat)
   rownames(eMat) <- paste0(dat$edat$geneExpr$gene_name,":",
                            dat$edat$geneExpr$gene_id)

   eMat <- eMat[apply(eMat,1,max) >= dat$specs$thresh$exprGenes.minTPM,]
   eMat <- eMat[apply(eMat,1,max) >= dat$specs$thresh$deGenes.FC *
                                     apply(eMat,1,min),]

   if (nlMode == "logMean") {
      bottomCap <- -1.5
      topCap <- 1.5

      eMat <- log2(eMat / apply(eMat, 1, mean))
      eMat[eMat < bottomCap] <- bottomCap; eMat[eMat > topCap] <- topCap;
      curBreaks <- seq(bottomCap, topCap, length.out = 101)
      plotTitle  <- "Rel. expression: log2 (TPM/mean)"
      plotColors <- colorRampPalette(c("steelBlue2", "white", "darkOrange2"))(100)
   } else if (nlMode == "minMax") {
      hMin <-  apply(eMat, 1, min)
      hMax <-  apply(eMat, 1, max)
      eMat <- (eMat - hMin) / (hMax - hMin)
      plotTitle  <- "(TPM - min)/(max - min))"
      curBreaks <- seq(0, 1, length.out = 101)
      plotColors <- colorRampPalette(c("steelBlue2", "white", "darkOrange2"))(100)
   } else if (nlMode == "fracMax") {
      hMax <-  apply(eMat, 1, max)
      eMat <- eMat / hMax
      plotTitle  <- "TPM / max"
      curBreaks <- seq(0, 1, length.out = 101)
      plotColors <- colorRampPalette(c("white", "red"))(100)
   }

   colnames(eMat) <- gsub("tpm.","",colnames(eMat))

   if (!missing(sampleAnnotate)) {
      colAnn <- data.frame( sampleNames = colnames(eMat),
                            stringsAsFactors = FALSE)
      rownames(colAnn) <- colnames(eMat)
      print(paste0("ROWNAMES COLANN=",rownames(colAnn)))
   
      if (!missing(sampleAnnotate)) {
         for (curProp in sampleAnnotate) {
            print(paste0("NOW ON ",curProp))
            colAnn[,curProp] <- dat$specs$rnaSamples[
               match(rownames(colAnn),
                     dat$specs$rnaSamples$sampleName),
               curProp]
         }
      }
      colAnn$sampleNames <- NULL
   }

   pheatmap.options <- list(eMat, fontsize_row = 4, fontsize_col = 4,
            main = plotTitle,
            cluster_rows = TRUE,
            cluster_cols = TRUE,
            border_color = NA,
            show_rownames = FALSE,
            show_colnames = TRUE,
            breaks = curBreaks,
            color = plotColors
            )

   if (!missing(sampleAnnotate)) {
#      print(colAnn)
      pheatmap.options$annotation_col <- colAnn
#      pheatmap.options$annotation_row <- colAnn
#      pheatmap.options$annotation_legend <- TRUE
      pheatmap.options$annotation_legend <- FALSE
      pheatmap.options$annotation_names_col <- TRUE
   }

   pdf(paste0(dat$specs$outDir, "/",
           figName, "_varGene_heatmap.pdf"),
       onefile = FALSE,
       height = 5,
       width  = 6)
   do.call(pheatmap, pheatmap.options)
   dev.off()
}


plotHmap.deg <- function(dat,
                         nlMode = "logMean",
                         geneList, #used if specified; otherwise picks DEGs
                         degTypes, #if not specified, uses all DEG's
                         sampleList, #used if specified; otherwise uses DEG samples
                         sampleAnnotate, #properties to show as sample annotation
                         show_rownames = TRUE,
                         clusterCols = FALSE,
                         clusterRows = TRUE,
                         rowFontSize = 2,
                         curHeight = 5,
                         curWidth = 5,
                         rowGaps,
                         repAvg, #named list of sample groups, if specified, show average across groups, rather than individual samples
                         figName = "deg_hmap"
                         ) {

# PURPOSE: make heatmap of bulk reporter RNA-seq data from lung

# By Default: show DEG separately for all pairwise comparisons

   curFig <- figName

# Genes: if geneList specified, use it;
#        else if degTypes specified, use only specified DEG's
#        else, use all DEGs
#
# Samples: if sampleList specified, use it;
#          else, if degTypes specified, use only those samples
#          else, use all DEG compared samples

   if (missing(geneList)) {
      if (missing(degTypes)) {
         degTypes <- names(dat$specs$degSamples)
      }
      geneList <- c()
      for (degType in degTypes) {
         geneList <- c(geneList, dat$deg[[degType]]$upGenes.minExpr,
                                 dat$deg[[degType]]$downGenes.minExpr)
      }
      geneList <- unique(geneList)
   }

   if (missing(sampleList)) {
      if (missing(degTypes)) {
         degTypes <- names(dat$specs$degSamples)
      }

      sampleList <- c()
      for (degType in degTypes) {
         sampleList <- c(sampleList, dat$specs$degSamples[[degType]][[1]],
                                     dat$specs$degSamples[[degType]][[2]])
      }
      sampleList <- sort(unique(sampleList))
   }
   tpmCols <- paste0("tpm.", sampleList)


   origExpr <- dat$edat$geneExpr[,c("gene_name","gene_id",tpmCols)]
   print(head(origExpr))


   print(geneList)

   hMat <- as.data.frame(origExpr[
      match(geneList[geneList %in% origExpr$gene_name], origExpr$gene_name),
      c("gene_name", "gene_id", tpmCols)])
   rownames(hMat) <- paste0(hMat$gene_name, ":", hMat$gene_id)

   rowLabels <- hMat$gene_name
   hMat$gene_name <- NULL
   hMat$gene_id <- NULL

   if (!missing(repAvg)) {
      newTpmCols <- c()
      for (celltype in names(repAvg)) {
         curTypeSamples <- repAvg[[celltype]]
         newTpmCol <- paste0("meantpm.",celltype)
         newTpmCols <- c(newTpmCols, newTpmCol)

         hMat[,newTpmCol] <- rowMeans(hMat[,paste0("tpm.",curTypeSamples)])
      }
      hMat <- hMat[,newTpmCols]
      print(newTpmCols)
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

   colnames(hMat) <- gsub("tpm.","",colnames(hMat))


   if (!missing(sampleAnnotate)) {
      colAnn <- data.frame( sampleNames = colnames(hMat),
                            stringsAsFactors = FALSE)
      rownames(colAnn) <- colnames(hMat)
      print("rownames(colAnn):")
      print(rownames(colAnn))

      for (curProp in sampleAnnotate) {
         print(paste0("NOW ON ",curProp))
         colAnn[,curProp] <- dat$specs$rnaSamples[
            match(rownames(colAnn),
                  dat$specs$rnaSamples$sampleName),
            curProp]
      }
      colAnn$sampleNames <- NULL
      print("Sample Annotations:")
      print(colAnn)
   }

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
            main = nlMode)

   if (missing(repAvg) & !missing(sampleAnnotate)) {
      pheatmap.options$annotation_col <- colAnn
      pheatmap.options$annotation_names_row <- TRUE

#      pheatmap.options$annotation_colors <- list(
#            celltype = c("TrTh2-Foxp3"  = "darkOrange2",
#                         "iTreg" = "steelBlue2",
#                         "Th2"     = "gray80"))
      pheatmap.options$annotation_legend <- TRUE
      pheatmap.options$annotation_names_col <- TRUE
   } else {
      pheatmap.options$show_colnames <- TRUE
   }

   outFn <- paste0(dat$specs$outDir, "/msFig", curFig, "_",nlMode,".pdf")
   pdf(outFn,
       onefile = FALSE,
       family = "ArialMT",
       height = curHeight,
       width  = curWidth)
   do.call(pheatmap, pheatmap.options)
   dev.off()

   curSampleInfo <- dat$specs$rnaSamples[dat$specs$rnaSamples$sampleName %in%
                                         sampleList,]
   figDescrFn <- paste0(dat$specs$outDir, "/", curFig, "_samples.txt")
   write.table(curSampleInfo, figDescrFn,
               row.names = FALSE, col.names = FALSE, quote = FALSE)

}



makeGeoExpressionTable <- function(dat, tabName = "geoTab1") {

   outSamples <- dat$specs$rnaSamples
   if (!"geoName" %in% colnames(outSamples)) {
      outSamples$geoName <- outSamples$sampleName
   }
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
      write.table(upGenes, file=upFn, quote=FALSE, row.names=FALSE,col.names=FALSE)
   }

   if (!(is.null(downGenes))){
      downFn<-paste0(outDir,"/downGenes.txt")
      write.table(downGenes, file=downFn,quote=FALSE, row.names=FALSE,col.names=FALSE)
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


plotDEGscatters <- function(dat, scatterPlot=TRUE) {

   for (i in 1:nrow(dat$specs$rnaComps)) {

      cell1 <- dat$specs$rnaComps[i,"group1name"]
      cell2 <- dat$specs$rnaComps[i,"group2name"]
      cell12 <- paste0(cell1,"_vs_",cell2)
      print(paste0("comparing ",cell1," to ",cell2))

      groupSamples <- dat$specs$degSamples[[cell12]]
      samples1   <- groupSamples[[1]]
      samples2   <- groupSamples[[2]]

      comparisonFC <- data.frame(
         gene_name      = dat$edat$geneExpr$gene_name,
         tpm.celltype2  = apply(dat$edat$geneExpr[, paste0("tpm.", samples2)],
                                   1, mean),
         tpm.celltype1  = apply(dat$edat$geneExpr[, paste0("tpm.", samples1)],
                                   1, mean),
         stringsAsFactors = FALSE
      )
      comparisonFC$log2fc.celltype1_vs_celltype2 <- round(
         (1 + comparisonFC$tpm.celltype1) / (1 + comparisonFC$tpm.celltype2),
         digits = 2 )
   
      exprRange <- range(1 + comparisonFC$tpm.celltype1,
                         1 + comparisonFC$tpm.celltype2)
   
      if (scatterPlot) {
      tmppngfn <- tempfile()
      png(file = tmppngfn, height = 3.1, width = 3.1, units = "in", res = 300,
          family = "ArialMT")
      par(mar = c(0, 0,0, 0))
      plot(1 + comparisonFC$tpm.celltype2,
           1 + comparisonFC$tpm.celltype1,
           ann=FALSE,axes=FALSE,
           pch=20, cex=0.5, log="xy", col="grey",las=1,
           xlab = paste0(cell2," (TPM + 1)"),
           ylab = paste0(cell1," (TPM + 1)"),
           xlim=exprRange,
           ylim=exprRange,
           main="")
      dev.off()
      pngbg <- readPNG(tmppngfn)
      pngbg <- as.raster(pngbg)
      
      outFn <- paste0(dat$specs$outDir,"/DEGscatter.",cell12,".pdf")
      pdf(outFn, height=3.5,width=3.5)
   
      par(mar = c(3.75, 3.75, 0.5, 0.5),
          mgp = c(2, 0.6, 0),
          cex = 1, cex.axis = 0.8, cex.lab = 1)
   
      plot(1 + comparisonFC$tpm.celltype2,
           1 + comparisonFC$tpm.celltype1,
           pch=20,
           cex=0.5,
           log="xy",
           type="n",
           las=1,
           xlab = paste0(cell2," (TPM + 1)"),
           ylab = paste0(cell1," (TPM + 1)"),
           xlim=exprRange,
           ylim=exprRange,
           main="")
      lim<-par()
      rasterImage(pngbg, 10^lim$usr[1], 10^lim$usr[3],
                         10^lim$usr[2], 10^lim$usr[4])
   
      degGenes.up <- comparisonFC[comparisonFC$gene_name %in%
                                  dat$deg[[cell12]]$upGenes.minExpr,]
   
      degGenes.down  <- comparisonFC[comparisonFC$gene_name %in%
                                     dat$deg[[cell12]]$downGenes.minExpr,]
                              
      points(1 + degGenes.up$tpm.celltype2,
             1 + degGenes.up$tpm.celltype1,
             pch=20,
             cex=0.5,
             col= "#998ec3")
   
      points(1 + degGenes.down$tpm.celltype2,
             1 + degGenes.down$tpm.celltype1,
             pch=20,
             cex=0.5,
             col="#f1a340")
   
   
      outliers.down <- head(degGenes.up[order(degGenes.up$log2fc.celltype1_vs_celltype2,decreasing=TRUE),], n = 10)
      outliers.up <- head(degGenes.down[order(degGenes.down$log2fc.celltype1_vs_celltype2),], n = 10)
   
      lab_x <- c(1 + outliers.down$tpm.celltype2, 1 + outliers.up$tpm.celltype2)
      lab_y <- c(1 + outliers.down$tpm.celltype1, 1 + outliers.up$tpm.celltype1)
      
      lab_text <- c(outliers.down$gene_name, outliers.up$gene_name)
      lab_textadj <- c(rep(1, nrow(outliers.down)), rep(0, nrow(outliers.up)))
      lab_textx <- c(rep(2.5, nrow(outliers.down)), rep(6500, nrow(outliers.up)))
      laborder <- order(lab_y)
      
      texty_logo <- ( 0.15 * log(max(exprRange), base = 10))
      texty_loginc <- ((0.75 * log(max(exprRange), base = 10)) / length(lab_y))
   
      for (j in 1:length(laborder)) {
         k <- laborder[j]
         lab_texty <- (10**( (j - 1) * texty_loginc + texty_logo))
         text(lab_textx[k], lab_texty, lab_text[k], cex = 0.5, adj = lab_textadj[k])
         segments(lab_textx[k], lab_texty, lab_x[k], lab_y[k], lwd = 1, col = "grey")
      }
      legend("topleft",
             paste0("n=",nrow(degGenes.up)," genes"),
             text.col="#998ec3",
             cex=0.75, bty="n")
      legend("bottomright",
             paste0("n=",nrow(degGenes.down)),
             text.col="#f1a340",
             cex=0.75, bty="n")
   
   
      dev.off()
      }
   
   }

}


plotDEGvolcanos <- function(dat) {

   for (i in 1:nrow(dat$specs$rnaComps)) {

      cell1 <- dat$specs$rnaComps[i,"group1name"]
      cell2 <- dat$specs$rnaComps[i,"group2name"]
      cell12 <- paste0(cell1,"_vs_",cell2)
      print(paste0("comparing ",cell1," to ",cell2))

      groupSamples <- dat$specs$degSamples[[cell12]]
      samples1   <- groupSamples[[1]]
      samples2   <- groupSamples[[2]]

      comparisonFC <- data.frame(
         transcript_id  = dat$edat$txExpr$transcript_id,
         gene_name      = dat$edat$txExpr$gene_name,
         tpm.celltype2  = apply(dat$edat$txExpr[, paste0("tpm.", samples2)],
                                   1, mean),
         tpm.celltype1  = apply(dat$edat$txExpr[, paste0("tpm.", samples1)],
                                   1, mean),
         stringsAsFactors = FALSE
      )

#      comparisonFC$fc.celltype1_vs_celltype2 <- round(
#         (1 + comparisonFC$tpm.celltype1) / (1 + comparisonFC$tpm.celltype2),
#         digits = 2 )

      comparisonFC$fc.celltype1_vs_celltype2 <- 
         (1 + comparisonFC$tpm.celltype1) / (1 + comparisonFC$tpm.celltype2)

#      return(list(comparisonFC=comparisonFC, sr=dat$deg[[cell12]]$sleuthResults))
      print("getting qvalue")
      comparisonFC$qval <- dat$deg[[cell12]]$sleuthResults[
                              match(comparisonFC$transcript_id,
                                    dat$deg[[cell12]]$sleuthResults$target_id),
                              "qval"]
      comparisonFC$qval.orig <- comparisonFC$qval
      comparisonFC$qval <- -1 * log10(comparisonFC$qval)
      print("->done")
      comparisonFC <- na.omit(comparisonFC)

      print(paste0("Restricting to genes at least ",dat$specs$thresh$exprGenes.minTPM," TPM"))
      comparisonFC  <- comparisonFC[comparisonFC$tpm.celltype1 >= 
                                    dat$specs$thresh$exprGenes.minTPM |
                                    comparisonFC$tpm.celltype2 >= 
                                    dat$specs$thresh$exprGenes.minTPM,]

      fcRange <- range(comparisonFC$fc.celltype1_vs_celltype2)
      qvalRange <- range(comparisonFC$qval)
      print("GOT HERE 1")
   
      tmppngfn <- tempfile()
      png(file = tmppngfn, height = 3.1, width = 3.1, units = "in", res = 300,
          family = "ArialMT")
      par(mar = c(0, 0,0, 0))
      plot(comparisonFC$fc.celltype1_vs_celltype2,
           comparisonFC$qval,
           ann=FALSE,axes=FALSE,
           pch=20, cex=0.5, log="x", col="grey",las=1,
           xlab = "",
           ylab = "-log10 q-value",
           xlim=fcRange,
           ylim=qvalRange,
           main="")
      dev.off()
      pngbg <- readPNG(tmppngfn)
      pngbg <- as.raster(pngbg)
      print("GOT HERE 2")
      
      outFn <- paste0(dat$specs$outDir,"/DEGvolcano.",cell12,".pdf")
      pdf(outFn, height=3.5,width=3.5)
   
      par(mar = c(3.75, 3.75, 0.5, 0.5),
          mgp = c(2, 0.6, 0),
          cex = 1, cex.axis = 0.8, cex.lab = 1)
   
      plot(comparisonFC$fc.celltype1_vs_celltype2,
           comparisonFC$qval,
           pch=20,
           cex=0.5,
           log="x",
           type="n",
           las=1,
#           xlab = paste0(cell12," (FC)"),
           xlab = "fold change",
           ylab = "-log10 q-value",
           xlim=fcRange,
           ylim=qvalRange,
           main="")
      lim<-par()
      rasterImage(pngbg, 10^lim$usr[1], lim$usr[3],
                         10^lim$usr[2], lim$usr[4])
      print("GOT HERE 3")
   
      if (1) {
      degGenes.up <- comparisonFC[comparisonFC$qval.orig <= dat$specs$thresh$sleuth.qval &
                                  comparisonFC$tpm.celltype1 > comparisonFC$tpm.celltype2,]
      degGenes.down <- comparisonFC[comparisonFC$qval.orig <= dat$specs$thresh$sleuth.qval &
                                    comparisonFC$tpm.celltype2 > comparisonFC$tpm.celltype1,]
                              
      points(degGenes.up$fc.celltype1_vs_celltype2,
             degGenes.up$qval,
             pch=20,
             cex=0.5,
             col= "#998ec3")
   
      points(degGenes.down$fc.celltype1_vs_celltype2,
             degGenes.down$qval,
             pch=20,
             cex=0.5,
             col="#f1a340")
   
   
      outliers.down <- head(degGenes.up[order(degGenes.up$qval,decreasing=TRUE),], n = 10)
      outliers.up <- head(degGenes.down[order(degGenes.down$qval,decreasing=TRUE),], n = 10)
   
      lab_x <- c(outliers.down$fc.celltype1_vs_celltype2,
                 outliers.up$fc.celltype1_vs_celltype2)
      lab_y <- c(outliers.down$qval, outliers.up$qval)
      
      lab_text <- c(outliers.down$gene_name, outliers.up$gene_name)
      lab_textadj <- c(rep(0, nrow(outliers.down)), rep(1, nrow(outliers.up)))
      lab_textx <- c(rep(25, nrow(outliers.down)), rep(0.04, nrow(outliers.up)))
      laborder <- order(lab_y)
      
      texty_logo <- ( 0.15 * max(qvalRange))
      texty_loginc <- ((0.8 * max(qvalRange)) / length(lab_y))
   
      for (j in 1:length(laborder)) {
         k <- laborder[j]
         lab_texty <- ( (j - 1) * texty_loginc + texty_logo)
         text(lab_textx[k], lab_texty, lab_text[k], cex = 0.5, adj = lab_textadj[k])
         segments(lab_textx[k], lab_texty, lab_x[k], lab_y[k], lwd = 1, col = "grey")
      }

      legend("topright",
             paste0("n=",nrow(degGenes.up)," transcripts"),
             text.col="#998ec3",
             cex=0.75,bty="n")

      legend("topleft",
             paste0("n=",nrow(degGenes.down)),
             text.col="#f1a340",
             cex=0.75,bty="n")
      }

      mtext(1,line=2,at=fcRange[1],adj=0,text=cell2,col="#f1a340")
      mtext(1,line=2,at=fcRange[2],adj=1,text=cell1,col="#998ec3")
   
      dev.off()
   
   }

}
