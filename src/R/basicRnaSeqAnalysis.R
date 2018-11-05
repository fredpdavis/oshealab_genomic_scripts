################################################################################
#
# basicRnaSeqAnalysis.R calls differentially expressed genes and makes figures
#
# Author:       Fred P. Davis, NIAMS/NIH. fred.davis@nih.gov
#
# To run the analyses, start R and run these commands:
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


paperRun <- function(dat.sc, dat.bulk, returnData = FALSE) {

   startup()

   if (missing(dat.sc)) {
      dat.sc <- absc.main();
      dat.sc <- absc.main(dat.sc) ;
   }

   if (missing(dat.bulk)) {
      dat.bulk <- bulk.main(returnData=TRUE);
      dat.bulk <- bulk.main(dat.bulk, returnData=TRUE);
      dat.bulk <- bulk.main(dat.bulk, returnData=TRUE);
      dat.bulk$specs$outDir <- dat.sc$paths$outDir
   }

   if (returnData) {
      return(list(dat.sc        = dat.sc,
                  dat.bulk      = dat.bulk))
   }

   makePaperFigures(list(dat.sc=dat.sc, dat.bulk=dat.bulk))

}

makePaperFigures <- function(dat, figList) {

   if (missing(figList)) {figList <- "all"}

# Fig     Desc                                            Data
# ------ ----------------------------------------------   --------------------
# 1A.     tsne. cell type labeled, colored                SC
   if (any(c("all","1A") %in% figList)) {
      absc.makeFig.tsneClusterLabels(dat$dat.sc, figName = "1A")
   }

# 1B.     heatmap. sc cluster-averaged marker levels      SC
   if (any(c("all","1B") %in% figList)) {
      absc.makeFig.clusterMarkerHeatmap(dat$dat.sc, figName="1B")
   }

# 1C.     violin. ctla4, tnfrsf18, foxp3                  SC
   if (any(c("all","1C") %in% figList)) {
      absc.makeFig.Violin( dat$dat.sc,
                           figName  = "1C",
                           geneList = c("Foxp3", "Ctla4", "Tnfrsf18"),
                           ident    = "celltype.lib",
                           coordFlip = TRUE,
                           addCellNum = FALSE,
                           point.size.use = FALSE,
                           cols.use = rev(brewer.pal(6,"Paired"))[c(1,3,5,2,4,6)],
                           ident.include = c("lung act Treg_Th.ova",
                                             "lung naive Treg_Th.ova",
                                             "lung naive_Th.ova",
                                             "lung act Treg_Th.pbs",
                                             "lung naive Treg_Th.pbs",
                                             "lung naive_Th.pbs"),
                           groupLabels = c("active Treg OVA",
                                           "naive Treg OVA",
                                           "naive OVA",
                                           "active Treg PBS",
                                           "naive Treg PBS",
                                           "naive PBS"),
                           height=2, width=2,
                           font.size=5 )
   }

# S1A.    tsne. library colored                           SC
   if (any(c("all","S1A") %in% figList)) {
      absc.makeFig.tsneAllLibs(dat$dat.sc, figName = "S1A")
   }

# S1B     tsne. gene Sell, Cd44, Foxp3                    SC
   if (any(c("all","S1B") %in% figList)) {
      absc.makeFig.tsneGeneMaps(dat$dat.sc, figName="S1B",
                             geneList = c("Sell","Cd44","Foxp3"))
   }

# 2D      heatmap. bulk reporter DEGs                     bulkRNA
# BACKBURN: col annotate part of the two signatures or not?
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


# S4      heatmap. nascent DEG                            nascentRNA
   if (any(c("all","S4") %in% figList)) {
       plotHmap.bulkReporter(dat$dat.bulk,
                             figName    = "S4",
                             nlMode     = "fracMax",
                             dataset    = "nascent",
                             repAvg        = TRUE,
                             clusterRows   = TRUE,
                             show_rownames = FALSE)
   }

# GEO supp expression table
   if (any(c("all","GEO_table") %in% figList)) {
      makeGeoExpressionTable(dat$dat.bulk)
   }

}

bulk.main <- function(dat, returnData = FALSE){

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

   if (! "sc" %in% names(dat)) {
      print("Loading SC data")
      dat <- scidr.loadExprSC(dat)
      if (returnData)    return(dat)
   }


   return(dat)

}



loadExpr <- function(specs, printTables = TRUE){

#   library(magrittr)
#   library(dplyr)

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



defineDEG <- function(dat){

# consistent DEG pairwise SLEUTH comparison

   celltypePairs <- list( list(dataType2="nascent",
                               cellType1 = "iTreg",
                               cellType2 = "TrTh2" ),
                          list(dataType2="nascent",
                               cellType1 = "Th2",
                               cellType2 = "TrTh2" ),
                          list(dataType2="nascent",
                               cellType1 = "Th2",
                               cellType2 = "iTreg"),
                          list(dataType2="bulk",
                               cellType1 = "lung.exFoxp3",
                               cellType2="lung.Foxp3p"),
                            list(dataType2 = "bulk",
                                 cellType1   = "lung.Foxp3p",
                                 cellType2   = "lung.nonFoxp3"),
                            list(dataType2 = "bulk",
                                 cellType1   = "lung.nonFoxp3",
                                 cellType2   = "lung.exFoxp3")
                         )

   genes.deg <- list()
   t2g <- data.frame(
      target_id = dat$specs$transcriptInfo$transcript_id,
      ens_gene  = dat$specs$transcriptInfo$gene_id,
      ext_gene  = dat$specs$transcriptInfo$gene_name,
      stringsAsFactors = FALSE
   )

   allSamples <- dat$specs$rnaSamples

   for (cellPair in celltypePairs) {
   
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
         
      if (cellPair$dataType2 == "nascent") {

         curSampleInfo$path <- paste0(dat$specs$kallistoBaseDir.stranded, "/",
                                      curSampleInfo$runID, "/",
                                      curSampleInfo$sampleName)
      } else {

         curSampleInfo$path <- paste0(dat$specs$kallistoBaseDir, "/",
                                      curSampleInfo$runID, "/",
                                      curSampleInfo$sampleName)

      }

# sleuth info https://rawgit.com/pachterlab/sleuth/master/inst/doc/intro.html
   
      print(paste0("comparing ",cellPair$cellType1," to ",cellPair$cellType2))
      cell12 <- paste0(cellPair$cellType1,"_vs_",cellPair$cellType2)
   
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
   
         sleuthData <- scoreSCsig.readSleuthOutput(resultsDir = curSleuthPath)
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
   
         scoreSCsig.writeSleuthOutput(
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



absc.main <- function(dat) {

   startup()

   curRdsFn <- "exfoxp3_sc_insurg.rds"

# 1. load full cellranger results
   if (missing(dat)) {
      dat <- absc.loadData(curRdsFn = curRdsFn)
      return(dat)
   }

# 2. run seurat to make all- and treg-only seurat objects
   if (!("all.so" %in% names(dat) &
         "treg.so" %in% names(dat))) {
      dat <- absc.processScData(dat)
      return(dat)
   }

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


absc.makeFig.cellTypeAbundance <- function(dat, figName = "4D_bar") {

# grouped barplots of cell type abundances
   ident <- "celltype.lib"

   goalOrder <- rev(c("lung naive",
                  "lung effector",
                  "lung naive Treg",
                  "lung act Treg",
                  "other"))

   plotLabels <- rev(c("naive",
                   "effector",
                   "naive Treg",
                   "active Treg",
                   "other"))

# Th.pbs
   nTable  <- table(dat$all.so@meta.data[[ident]])
   deleteCols <- c()
   for (libType in c("pbs", "ova")) {
      collapseCols <- c( paste0("proliferating_Th.",libType),
                         paste0("contam_Th.",libType),
                         paste0("mLN effector_Th.",libType),
                         paste0("mLN naive_Th.",libType),
                         paste0("mLN naive Treg_Th.",libType))

      nTable[paste0("other_Th.",libType)] <- sum(nTable[collapseCols])
      deleteCols <- c(deleteCols, collapseCols)
   }
   nTable <- nTable[!names(nTable) %in% deleteCols]

   pbs.nTable <- nTable[grepl("[^N]_Th.pbs", names(nTable))]
   ova.nTable <- nTable[grepl("[^N]_Th.ova", names(nTable))]

   pbs.nTable <- 100 * pbs.nTable / sum(pbs.nTable)
   ova.nTable <- 100 * ova.nTable / sum(ova.nTable)

# Th.ova
   plotDat <- rbind( data.frame(condition = rep("PBS", length(pbs.nTable)),
                                celltype = factor(gsub("_.*","", names(pbs.nTable))),
                                percCell = as.vector(pbs.nTable)),
                     data.frame(condition = rep("OVA", length(ova.nTable)),
                                celltype = factor(gsub("_.*","", names(ova.nTable))),
                                percCell = as.vector(ova.nTable)))
   print("BEFORE")
   print(plotDat)

   plotDat$celltype <- plyr::mapvalues(x        = plotDat$celltype,
                                       from     = goalOrder,
                                       to       = plotLabels)

   plotDat$celltype <- factor(plotDat$celltype,
                              levels = plotLabels)
   print("AFTER")
   print(plotDat)

   outFn <- paste0(dat$paths$outDir,"/msFig",figName,"_celltype_frequency.pdf")
   tx <- ggplot(plotDat,
                aes(x=celltype,
                    y=percCell,
                    fontsize=2,
                    fill=factor(condition)))+
         geom_bar(stat="identity",position="dodge",width=0.75)+
         scale_fill_discrete(name="",
                             breaks=c("OVA", "PBS"),
                             labels=c("OVA", "PBS"))+
         xlab("Cell type")+ylab("Relative abundance (%)") +
         scale_y_log10() +
         coord_flip()


   ggsave(filename=outFn,
          plot=tx,
          width=4,
          height=3,
          units="in",
          device=cairo_pdf)
}



absc.processScData <- function(dat) {

   library(Seurat)
   library(dplyr)
   library(Matrix)

   seuratRobjFn.all <- paste0(dat$paths$outDir,"/seurat_all.Robj")
   seuratRobjFn.treg <- paste0(dat$paths$outDir,"/seurat_treg.Robj")

   if (file.exists(seuratRobjFn.all)) {

      print("Loading seurat object (all libraries) from disk")
      load(seuratRobjFn.all)

   } else {
   
      print("Read10X()")
      all.so <- Read10X(data.dir =
                            paste0(dat$paths$sc.dataDir,
                                   "/outs/filtered_gene_bc_matrices_mex/mm10_egfp/"))
   
      print("CreateSeuratObject()")
      all.so <- CreateSeuratObject(raw.data = all.so,
                                   min.cells = 3,
                                   min.genes = 200,
                                   scale.factor = 10000)
   
      print(paste0("Number of cells in all libraries: ", length(all.so@cell.names)))
   
      print("*************** NOTE! skipping all-library tsne for testing *********************")
      print("NormalizeData()")
      all.so <- NormalizeData(object = all.so,
                              normalization.method = "LogNormalize",
                              scale.factor = 10000)
   
      print("ScaleData()")
      all.so <- ScaleData(object = all.so, vars.to.regress = c("nUMI"))
   
      print("FindVariableGenes()")
      all.so <- FindVariableGenes(object = all.so, mean.function=ExpMean,
                                  do.plot = FALSE,
                                  dispersion.function=LogVMR,
                                  x.low.cutoff = 0.0125,
                                  x.high.cutoff = 3,
                                  y.cutoff = 0.5)
   
      print("RunPCA()")
      all.so <- RunPCA(object = all.so,
                       pc.genes = all.so@var.genes,
                       do.print = TRUE, pcs.print = 1:10,
                       genes.print = 10)
   
      print("FindClusters()")
      all.so <- FindClusters(object = all.so, reduction.type = "pca", dims.use = 1:10, 
                             resolution = c(1.2),
                             print.output = 0, save.SNN = TRUE)
      PrintFindClustersParams(object = all.so)

# Relabel clusters
      # Remove library name suffix
      print("SetAllIdent()")
      all.so@meta.data$clean_ident <- gsub("_.*","",as.character(all.so@ident))
      all.so <- SetAllIdent(all.so, id="clean_ident")

#       0,3,1,8: blue – lung naive
#       10, 5, 2, 14: light blue – mLN naive
#       11, 4: orange – lung effector
#       12: light orange – mLN effector
#       15: pink - proliferating
#       6: dark green –lung naïve Treg
#       9: dark yellow-green – mLN naïve Treg
#       7: light green – lung activated Treg

      print("Manually labeling clusters")
      current.cluster.ids <- 0:16
      new.cluster.ids <- c("lung naive", "lung naive", "mLN naive", "lung naive",
                           "lung effector", "mLN naive", "lung naive Treg", "lung act Treg",
                           "lung naive", "mLN naive Treg", "mLN naive", "lung effector",
                           "mLN effector", "lung naive", "mLN naive", "proliferating", "contam")
      all.so@meta.data$orig_cluster_id <- all.so@ident
      all.so@ident <- plyr::mapvalues(x = all.so@ident, from = current.cluster.ids, to = new.cluster.ids)
      all.so@meta.data$orig.ident <- as.factor(gsub(".*-","",all.so@cell.names))

      all.so@meta.data$libraryName <- dat$sc$libNames[as.numeric(gsub(".*-","",all.so@cell.names))]

      print("Additional new ident: celltype.lib")
      all.so@meta.data$celltype.lib <- paste0(all.so@ident,"_",all.so@meta.data$libraryName)
      print(paste0("NEW IDENT: ", head(all.so@meta.data$celltype.lib)))
      all.so <- StashIdent(all.so, save.name = "celltype")
      all.so <- SetAllIdent(all.so, id = "celltype.lib")

   
      print("RunTSNE()")
      all.so <- RunTSNE(object = all.so, dims.use=1:10, do.fast=TRUE)
      pdf(paste0(dat$paths$outDir,"/seurat_alllib_tsne.pdf"))
      TSNEPlot(object = all.so, do.label=TRUE)
      dev.off()
   
      print("Saving seurat object (all) to disk")
      save(all.so, file=seuratRobjFn.all)
   
   }


   if (file.exists(seuratRobjFn.treg)) {

      print("Loading seurat object (Treg only) from disk")
      load(seuratRobjFn.treg)

   } else {

      print("Process just Treg libraries")
   
      print("SubsetData()")
      treg.cells <- all.so@cell.names[gsub(".*-","", all.so@cell.names) %in%
                                      which(grepl("Treg",dat$sc$libNames))]
   
      treg.so <- SubsetData(all.so, cells.use = treg.cells, do.clean=TRUE)
      print(paste0("Number of cells in all libraries: ", length(treg.so@cell.names)))

      print("NormalizeData()")
      treg.so <- NormalizeData(object = treg.so,
                              normalization.method = "LogNormalize",
                              scale.factor = 10000)
   
      print("ScaleData()")
      treg.so <- ScaleData(object = treg.so, vars.to.regress = c("nUMI"))
   
   
      print("FindVariableGenes()")
      treg.so <- FindVariableGenes(object = treg.so, mean.function=ExpMean,
                                  do.plot = FALSE,
                                  dispersion.function=LogVMR,
                                  x.low.cutoff = 0.0125,
                                  x.high.cutoff = 3,
                                  y.cutoff = 0.5)
   
      print("RunPCA()")
      treg.so <- RunPCA(object = treg.so,
                       pc.genes = treg.so@var.genes,
                       do.print = TRUE, pcs.print = 1:10,
                       genes.print = 10)
   
      print("FindClusters()")
      treg.so <- FindClusters(object = treg.so, reduction.type = "pca", dims.use = 1:10, 
                             resolution = c(1),
                             print.output = 0, save.SNN = TRUE)
      PrintFindClustersParams(object = treg.so)

# Relabel clusters
      current.cluster.ids <- 0:8
      new.cluster.ids <- c("naive",
                           "activated",
                           "naive",
                           "naive",
                           "activated",
                           "activated",
                           "activated",
                           "proliferating",
                           "contam")

      treg.so@meta.data$orig_cluster_id <- treg.so@ident
      treg.so@ident <- plyr::mapvalues(x = treg.so@ident, from = current.cluster.ids, to = new.cluster.ids)
   
      print("RunTSNE()")
      treg.so <- RunTSNE(object = treg.so, dims.use=1:10, do.fast=TRUE)
      pdf(paste0(dat$paths$outDir,"/seurat_treg_tsne.pdf"))
      TSNEPlot(object = treg.so, do.label=TRUE)
      dev.off()
   
      print("Saving seurat object (Treg only) to disk")
      save(treg.so, file=seuratRobjFn.treg)
   }

   dat$all.so <- all.so
   dat$treg.so <- treg.so

#add foxp3.plus.egfp entry

   tx2 <- rbind(dat$all.so@data, foxp3.plus.egfp=dat$all.so@data["Foxp3",] + dat$all.so@data["egfp",])
   dat$all.so@data <- tx2

   tx2 <- rbind(dat$all.so@scale.data, foxp3.plus.egfp=dat$all.so@scale.data["Foxp3",] + dat$all.so@scale.data["egfp",])
   dat$all.so@scale.data <- tx2



   return(dat)

}


absc.makeFig.tsneAllLibs <- function(dat, figName = "S1A") {

   print("Plotting labeled tSNE")
   for (noLegend in c(TRUE, FALSE)) {
      pdf(paste0(dat$paths$outDir,"/msFig",figName,
               "_seurat_alllib_tsne_library_nolegend",
                  as.character(noLegend),".pdf"),
         height=3,width=3)
      TSNEPlot(object = dat$all.so, group.by="libraryName",
               no.legend=noLegend,
               pt.size=0.5, vector.friendly=TRUE)
      dev.off()
   }

}

absc.subCluster <- function(dat,
                            clusterName = "lung effector",
                            clusterRes = 1.5,
                            ident = "celltype",
                            nMarker = 10,
                            figName="5A") {

# given an existing cluster (eg, lung effectors), recluster at
#   to define subclusters

   goalCells <- dat$all.so@cell.names[dat$all.so@meta.data[[ident]]== clusterName]
   cur.so <- SubsetData(dat$all.so, cells.use=goalCells, do.clean=TRUE)
   print(paste0("Number of selected cells in cluster ", clusterName,
                ": ", length(cur.so@cell.names)))

   print("NormalizeData()")
   cur.so <- NormalizeData(object = cur.so,
                           normalization.method = "LogNormalize",
                           scale.factor = 10000)

   print("ScaleData()")
   cur.so <- ScaleData(object = cur.so, vars.to.regress = c("nUMI"))


   print("FindVariableGenes()")
   outFn <- paste0(dat$paths$outDir,"/msFig",figName,
                   "_",clusterName,"_vargenes.pdf")
   outFn <- gsub(" ", "_", outFn)
   pdf(outFn)
   cur.so <- FindVariableGenes(object = cur.so, mean.function=ExpMean,
#                               do.plot = FALSE,
                               dispersion.function=LogVMR,
                               x.low.cutoff = 0.0125,
                               x.high.cutoff = 3,
                               sort.results=TRUE,
                               y.cutoff = 0.5)
   dev.off()
   return(cur.so)

   print("RunPCA()")
   cur.so <- RunPCA(object = cur.so,
                    pc.genes = cur.so@var.genes,
                    do.print = TRUE, pcs.print = 1:10,
                    genes.print = 10)

   print("FindClusters()")
   cur.so <- FindClusters(object = cur.so, reduction.type = "pca", dims.use = 1:10, 
                          resolution = clusterRes,
                          print.output = 0, save.SNN = TRUE)
   PrintFindClustersParams(object = cur.so)


   print("Assigning new subcluster identities")
# Paint new subclusters on original all cell TSNE
   newIdents <- dat$all.so@meta.data[[ident]]
   newIdents[newIdents != clusterName] <- "other"
   newIdents[match(names(cur.so@ident), dat$all.so@cell.names)] <- cur.so@ident
   dat$all.so@meta.data$subcluster <- newIdents

   print("TSNEPlot()")
   outFn <- paste0(dat$paths$outDir,"/msFig",figName,
                   "_subcluster_res_",clusterRes,
                   "_",clusterName,"_tsne.pdf")
   outFn <- gsub(" ", "_", outFn)
   pdf(outFn)
   try(TSNEPlot(object = dat$all.so, group.by="subcluster", pt.size=0.25, do.label=TRUE))
   dev.off()


# Heatmap of subcluster markers

   print("FindMarkers()")
   markers <- FindAllMarkers(object = cur.so,
                             only.pos = TRUE,
                             min.pct = 0.25,
                             return.thresh = 0.05, #5% FDR
                             thresh.use = 0.25)
   top10 <- markers %>% group_by(cluster) %>% top_n(nMarker, avg_logFC)
   geneList <- top10$gene
   
   outFn <- paste0(dat$paths$outDir,"/msFig",figName,
                   "_subcluster_res_",clusterRes,
                   "_",clusterName,"_markers.pdf")
   outFn <- gsub(" ", "_", outFn)
   hmap <- DoHeatmap(object = cur.so,
             genes.use = geneList,
             slim.col.label = TRUE, 
             group.label.loc="top",
             remove.key = TRUE, cex.row=4,
             do.plot=FALSE)
#   dev.off()

    try(ggplot2::ggsave(filename=outFn,
                    plot = hmap,
                    width = 8,
                    height = 5,
                    units = "in",
                    device=cairo_pdf))

   return(list(subset.so = cur.so))

# eg, given lung effector Th cells, find subclusters (and viz on tsne)

}


absc.makeFig.tsneGeneMaps <- function(dat, figName = "1B",
                                 geneList) {
   print("Plotting genes on tsne plot")

   if (missing(geneList)) {
      geneList <- c("Cd28", "Gata3","Il4", "Il5", "Il13")
#      geneList <- c("Anxa2", "Lgals1")
#      geneList <- c("Il2ra", "Il2rb", "Icos", "Nt5e")
#                  c( "Sell", "Cd44", "Foxp3", "Glrx", "egfp",
#                     "Dusp10", "Dusp2", "Mki67", "Igfbp4")
#                     "foxp3.plus.egfp" )
   }

   for (gene in geneList) {
      pdf(paste0(dat$paths$outDir,"/msFig",figName,"_seurat_alllib_tsne_gene_",gene,".pdf"),
          height=2.5,width=2.5)
      FeaturePlot(object = dat$all.so,
                  features.plot = gene,
                  pt.size=1,
                  cols.use = c("grey", "blue"), 
                  reduction.use = "tsne",
                  vector.friendly=TRUE)
      dev.off()
   }
}




absc.makeFig.tsneClusterLabels <- function(dat, figName = "1A") {

   print("Plotting labeled tSNE")
   for (noLegend in c(TRUE, FALSE)) {
      pdf(paste0(dat$paths$outDir,"/msFig",figName,
                 "_seurat_alllib_tsne_cluster_nolegend",
                 as.character(noLegend),".pdf"),
         height=3,width=3)
      TSNEPlot(object = dat$all.so, group.by="celltype",
               no.legend=noLegend, pt.size=0.5, vector.friendly=TRUE)
      dev.off()
   }

   print("Plotting labeled tSNE (original cluster)")
   pdf(paste0(dat$paths$outDir,"/msFig",figName,"_seurat_alllib_tsne_original_cluster_label.pdf"),
       height=3,width=3)
   TSNEPlot(object = dat$all.so, group.by="orig_cluster_id", pt.size=0.5, vector.friendly=TRUE)
   dev.off()

}


absc.makeFig.clusterMarkerHeatmap <- function(dat,
                                         figName = "S1F",
                                         geneList,
                                         plot=TRUE,
                                         nMarker = 10) {

   cur.so <- SetAllIdent(dat$all.so, id="celltype")
   origClusterOrder <- names(table(cur.so@ident))
   goalClusterOrder <- c("lung naive",
                         "mLN naive",
                         "lung effector",
                         "lung naive Treg",
                         "lung act Treg",
                         "mLN naive Treg",
                         "mLN effector",
                         "proliferating")

   if (missing(geneList)) {

      allMarkersFn <- paste0(dat$paths$outDir,"/",figName,"_marker_genes_celltype_all.txt")
      markersFn <- paste0(dat$paths$outDir,"/",figName,"_marker_genes_celltype_top10.txt")

      if (file.exists(allMarkersFn)) {
         markers <- read.table(allMarkersFn, header=TRUE,sep="\t", stringsAsFactors=FALSE,as.is=TRUE)
      } else {
         print("FindMarkers()")
         markers <- FindAllMarkers(object = cur.so,
                                   only.pos = TRUE,
                                   min.pct = 0.25,
                                   return.thresh = 0.05, #5% FDR
                                   thresh.use = 0.25)

         write.table(markers, file=allMarkersFn, row.names=FALSE,col.names=TRUE,quote=FALSE,sep="\t")
      }

      markers <- markers[!markers$gene %in% "Foxp3",]

      top10 <- markers %>% group_by(cluster) %>% top_n(nMarker, avg_logFC)
      clusterGeneList <- top10[,c("cluster","gene")]
      geneList <- top10$gene
      print("-> Saving marker gene list")
      write.table(clusterGeneList,
                     file=markersFn,
                     row.names=FALSE,
                     col.names=TRUE,
                     quote=FALSE,sep="\t")
   }

   rowGaps <- c()
   rowTypes <- c()
   if (any(goalClusterOrder != origClusterOrder)) {
      geneList <- c()
      for (cluster in goalClusterOrder) {
         geneList <- c(geneList, clusterGeneList$gene[clusterGeneList$cluster == cluster])
         rowTypes <- c(rowTypes, rep(cluster, sum(clusterGeneList$cluster == cluster)))
         rowGaps <- c(rowGaps, length(geneList))
         print(geneList)
      }
   }


   if (0 & plot) {
   outFn <- paste0(dat$paths$outDir,"/msFig",figName,
                   "_seurat_alllib_markerHeatmap.pdf")
#   pdf(outFn, width=11,height=8.5)
   hmap <- DoHeatmap(object = dat$all.so,
             genes.use = geneList,
             slim.col.label = TRUE, 
             group.by = "celltype",
             group.order=c("lung naive",
                           "mLN naive",
                           "lung effector",
                           "lung naive Treg",
                           "lung act Treg",
                           "mLN naive Treg",
                           "mLN effector",
                           "proliferating",
                           "contam"),
#             group.by = "orig_cluster_id",
             group.label.loc="top",
             group.label.rot=TRUE,
             remove.key = TRUE, cex.row=6,
             do.plot=FALSE)
#   dev.off()

    ggplot2::ggsave(filename=outFn,
                    plot = hmap,
                    width = 8,
                    height = 5,
                    units = "in",

                    device=cairo_pdf)
   }


   if (1) {
      library(pheatmap)
      library(RColorBrewer)
      outFn <- paste0(dat$paths$outDir,"/msFig",figName,
                     "_seurat_alllib_markerHeatmap_clusteravg.pdf")
      tx <- AverageExpression(cur.so, use.scale=TRUE)

      tx <- tx[c("lung naive",
                 "mLN naive",
                 "lung effector",
                 "lung naive Treg",
                 "lung act Treg",
                 "mLN naive Treg",
                 "mLN effector",
                 "proliferating") ]

      if (1) {
         rowLabels <- geneList
         rowLabels[rowLabels == "foxp3.plus.egfp"] <- "Foxp3 (plus eGFP)"
         duplLabel <- unique(rowLabels[duplicated(rowLabels)])
         rowLabels[rowLabels %in% duplLabel] <- paste0(rowLabels[rowLabels %in% duplLabel]," *")
         plotMat <- tx[geneList,]
         rownames(plotMat) <- 1:nrow(plotMat)


         colAnn <- data.frame( celltype=colnames(plotMat))
         rownames(colAnn) <- colnames(plotMat)

         rowAnn <- data.frame( celltype=rowTypes)
         rownames(rowAnn) <- 1:nrow(rowAnn)

         annColors <- rev(brewer.pal(n = 8, name = "Dark2"))
         names(annColors) <- colnames(plotMat)

         pdf(outFn,height=5,width=2.5)
         pheatmap(plotMat,
                  show_rownames=TRUE,
                  show_colnames=TRUE,
                  cluster_cols=FALSE,
                  cluster_rows=FALSE,
                  annotation_col = colAnn,
                  annotation_row = rowAnn,
                  annotation_names_row = FALSE,
                  annotation_names_col = FALSE,
                  annotation_colors = list(celltype = annColors),
                  annotation_legend = FALSE,
                  labels_row = rowLabels,
                  gaps_row = rowGaps,
                  scale="row",
                  treeheight_row = 0,
                  treeheight_col = 0,
                  color          = colorRampPalette(c("steelBlue2", "white",
                                                      "darkOrange2"))(100),
                  fontsize_row = 3,
                  fontsize_col = 6,
                  border_color = NA,
                  main = "")
         dev.off()
         return(plotMat)
      }


      if (0) {
         hmap <- DoHeatmap(cur.so,
                           data.use = tx,
                           genes.use = geneList,
                           cells.use = colnames(tx),
             group.label.loc="top",
             group.label.rot=TRUE,
                           slim.col.label=TRUE,
                           remove.key = TRUE, cex.row=4,
                           do.plot=FALSE)

    ggplot2::ggsave(filename=outFn,
                    plot = hmap,
                    width = 4,
                    height = 4,
                    units = "in",

                    device=cairo_pdf)
      }
   }

   return(geneList)
}

absc.makeFig.Violin <- function(dat, figName="1Y",
                                geneList = c("foxp3.plus.egfp"),
                                ident = "celltype",
                                height = 4,
                                width = 4,
                                addCellNum = TRUE,
                                font.size = 1,
                                point.size.use=0.025,
                                cols.use = NULL,
                                coordFlip = FALSE,
                                xLabel,
                                groupLabels,
                                ident.include.order = TRUE, #use ident.include order
                                ident.include) {

   if (missing(xLabel)) {xLabel <- "Expression level"}

   if (coordFlip & !missing(ident.include)) {
      ident.include <- rev(ident.include)
      if (!missing(groupLabels)) { groupLabels <- rev(groupLabels)}
   }

   if (coordFlip & !is.null(cols.use)) {
      cols.use <- rev(cols.use) }

   cur.so <- SetAllIdent(dat$all.so, id=ident)

   if (!missing(ident.include) & ident.include.order) {
      cur.so@ident <- factor(cur.so@ident,
                             levels=c(ident.include,
                                      setdiff(levels(cur.so@ident), ident.include)))

      if (!missing(groupLabels)) {
         cur.so@ident <- plyr::mapvalues(x      = cur.so@ident,
                                         from   = ident.include,
                                         to     = groupLabels)
      }
      ident.include <- groupLabels

   }

   countTable <- table(cur.so@ident)
   if (missing(ident.include)) { ident.include <- names(table(cur.so@ident)) }
   countTable <- as.data.frame(countTable[ident.include])
   colnames(countTable) <- c("cluster","n")

   if (addCellNum) {
   if (!coordFlip) {
      countTable$n[1] <- paste0("n=",countTable$n[1])
   } else {
      countTable$n[nrow(countTable)] <- paste0("n=",countTable$n[nrow(countTable)])
   }
   }


   for (gene in geneList) {
      outFn <- paste0(dat$paths$outDir,"/msFig",figName,
                     "_seurat_alllib_violin_",ident,"_",gene,".pdf")

      VlnPlot.options <- list(cur.so,
                    feature     = gene,
                    do.return   = TRUE,
#                    x.lab.rot   = TRUE,
                    point.size.use=point.size.use,
                    ident.include=ident.include,
                    size.title.use=font.size,
                    size.y.use  = font.size,
                    size.x.use  = font.size)

      if (!missing(cols.use)) {
         VlnPlot.options$cols.use <- cols.use
      }
      tx <- do.call(VlnPlot, VlnPlot.options)
#      tx <- tx + coord_cartesian(clip='off')
      if (coordFlip) {

         if (addCellNum) {yMaxMult <- 1.25} else {yMaxMult <- 1}

         tx <- tx + coord_flip(clip='off',
                               ylim = c(
                                 min(ggplot_build(tx)$layout$panel_scales_y[[1]]$range$range),
                      yMaxMult * max(ggplot_build(tx)$layout$panel_scales_y[[1]]$range$range)))
         vjust <- 0.5
         hjust <- 1
      } else {
         vjust <- 1
         hjust <- 0.5
      }
      
      if (addCellNum) {
      tx <- tx + geom_text(data = countTable,
                           aes(x=cluster,y=Inf,label=n),
                           color="black",
                           size=font.size / 3,
                           vjust=vjust,
                           hjust=hjust)
      }

      tx <- tx + theme(text = element_text(size=font.size),
                       axis.title  = element_text(size=font.size),
                       axis.title.y  = element_blank(),
#                       axis.title.x  = element_blank(),
                       axis.title.x  = element_text(color="black"),
                       axis.text.x = element_text(size=font.size),
                       axis.text.y = element_text(size=font.size)) +
                 ylab(xLabel)
#                       plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "in"))

      try(ggplot2::ggsave(filename=outFn,
                    plot = tx,
                    width = width,
                    height = height,
                    units = "in",
                    device=cairo_pdf))
      print("GOT HERE")
   }

}


absc.makeFig.allSigScore <- function(dat, figName = "4Y",
                                     sigScore = "iTreg_vs_TrTh2",
                                     flipScore = FALSE,
                                     scoreLabel
                                     ) {

# Load Treg seurat object with loadData() to calculate and display sig scores
   dat.seurat <- absc.loadData(curRdsFn = "seurat_alllib.rds",
                          seuratFn = paste0("analysis/20180904.fisher_bulkDEG",
                                            "/seurat_all.Robj"),
                          libNames = c("1" = "Treg.ova",
                                       "2" = "Treg.pbs",
                                       "3" = "Th.ova",
                                       "4" = "Th.pbs",
                                       "5" = "Th.mLN"))

   tx <- absc.cleanDEG.exFoxp3(dat.seurat, f_c = 0.3)

   if (missing(scoreLabel)) { scoreLabel <- sigScore }

   axisType <- "sleuth"
   if (sigScore == "pantreg_mathis") {
      axisType <- "geneList"
      tx <- readDispiritoTregSig(baseDir = dat$paths$baseDir)
      dat.seurat$specs$deg.geneLists <- list()
      dat.seurat$specs$deg.geneLists[[sigScore]] <- list(
         "upGenes"   = tx[["pan-tissue"]]$up,
         "downGenes" = tx[["pan-tissue"]]$down)
   } else if ( sigScore == "lung.nonFoxp3_vs_lung.exFoxp3.decontam" ) {
      dat.seurat$deg[["lung.nonFoxp3_vs_lung.exFoxp3.decontam"]] <- 
         dat.seurat$deg[["lung.nonFoxp3_vs_lung.exFoxp3"]]

      dat.seurat$deg[["lung.nonFoxp3_vs_lung.exFoxp3.decontam"]]$upGenes.minExpr <- tx$newUpGenes
      dat.seurat$deg[["lung.nonFoxp3_vs_lung.exFoxp3.decontam"]]$downGenes.minExpr <- tx$newDownGenes
   }       

   tx <- scoreSCsig.scoreSConBulkDEG(dat.seurat,
                                     axis1 = sigScore,
                                     axisTypes = axisType,
                                     figName = figName,
                                     score = "zscore",
                                     tsne.minCap = -5, tsne.maxCap = 5,
                                     minTPMthresh = 150,
                                     printModule = TRUE,
                                     tsne.ptCex = 0.5,
                                     separateClusters = TRUE);

   curDat <- list(all.so = dat.seurat$sc$origSO,
                  paths = dat$paths)
   curDat$all.so@meta.data[[scoreLabel]] <- tx$degScores[[1]][match(names(curDat$all.so@ident),
                                                                  names(tx$degScores[[1]]))]
   names(curDat$all.so@meta.data[[scoreLabel]]) <- curDat$all.so@cell.names

   if (flipScore) {
      curDat$all.so@meta.data[[scoreLabel]] <- -1 *
         curDat$all.so@meta.data[[scoreLabel]] }

   seuratCols <- rebalanceColorMap(
      scoreRange = range(curDat$all.so@meta.data[[scoreLabel]]),
      minCap = -5,
      maxCap = 5,
      colorScale = c("#2c7bb6", "#ffffbf", "#d7191c")
#      colorScale = c("steelBlue2", "white", "darkOrange2")
      )

   for (libName in c("Th.ova","Th.pbs")) {
      print("PLOTTING SCORE TSNE USING SEURAT INTERNAL!!")
      outFn <- paste0(dat$paths$outDir,"/msFig",figName,
                      "_seurat_score_tsne_",scoreLabel,
                      "__library_",libName,".pdf")
      outFn <- gsub(" ", "_", outFn)
      pdf(outFn, height=2.5, width=2.5)
      FeaturePlot(object=curDat$all.so,
                  cells.use = curDat$all.so@cell.names[curDat$all.so@meta.data$libraryName == libName], 
                  features.plot=scoreLabel,
                  pt.size=2,
   #               no.legend=FALSE,
                  cols.use= seuratCols$colors,
   #               cols.use= colorRampPalette(c("steelBlue2", "white",
   #                                            "darkOrange2"))(100),
                  reduction.use="tsne",
                  vector.friendly=TRUE)
      dev.off()

      outFn <- paste0(dat$paths$outDir,"/msFig",figName,
                      "_seurat_score_tsne_",scoreLabel,
                      "__library_",libName,".legend.pdf")
      pdf(outFn, height=2.5,width=2.5)
      par(mar=c(1,4,0,2))
      plot(c(0,1),
           c(seuratCols$scaleMin,seuratCols$scaleMax),
           col="white",
           bty="n",
#           type = 'n',
#           axes = F,
           xlab = '', ylab = '', main = "")
#      axis(3)
#      text(x=1.5,
#           y=seq(0,seurat,l=5),
#           labels=seq(seuratCols$scaleMin,seuratCols$scaleMax,l=5),
#           cex=0.5)
      rasterImage(as.raster(matrix(rev(seuratCols$colors),ncol=1)),
                  0,seuratCols$scaleMin,
                  0.5,seuratCols$scaleMax)
      dev.off()
   }

#   return(list(curDat = curDat,
#               tx = tx))

# Make violin plot showing score along celltype.lib axes
   absc.makeFig.Violin(curDat,
#                       geneList = "Foxp3",
                       figName = figName,
                       geneList = scoreLabel,
                       xLabel   = "Signature score (Z)",
                       addCellNum = FALSE,
                       ident    ="celltype.lib",
                       ident.include = c("lung act Treg_Th.ova",
                                         "lung act Treg_Th.pbs",
                                         "lung naive Treg_Th.ova",
                                         "lung naive Treg_Th.pbs",
                                         "lung effector_Th.ova",
                                         "lung effector_Th.pbs"),
#                                         "lung naive_Th.ova",
#                                         "lung naive_Th.pbs"),
                       groupLabels   = c("active Treg OVA",
                                         "active Treg PBS",
                                         "naive Treg OVA",
                                         "naive Treg PBS",
                                         "effector OVA",
                                         "effector PBS"),
#                                         "naive OVA",
#                                         "naive PBS"),
                       coordFlip=TRUE,
                       cols.use = rev(brewer.pal(6,"Paired")),
                       point.size.use = FALSE,
                       height=2, width=2,
                       font.size=5)
#   return(curDat)
#
   return(list(dat.seurat = dat.seurat, degScores = tx$degScores))

}


absc.loadData <- function(curRdsFn,
                     seuratFn,
                     libNames,
                     experimentName) {

# Purpose: loads either cellranger or seurat datasets

   if (!missing(curRdsFn) & file.exists(curRdsFn)) {
      dat <- readRDS(curRdsFn)
      return(dat)
   }

   if (missing(experimentName)) {
      experimentName <- "all_Th"
   }

   deg.celltypePairs = list(list(libraryType = "bulk",
                                 cellType1   = "iTreg",
                                 cellType2   = "TrTh2" ),
                            list(libraryType = "bulk",
                                 cellType1   = "Th2",
                                 cellType2   = "TrTh2" ),
                            list(libraryType = "bulk",
                                 cellType1   = "iTreg",
                                 cellType2   = "Th2" ),
                            list(libraryType = "bulk",
                                 cellType1   = "TrTh2",
                                 cellType2   = "iTreg"),
                            list(libraryType = "bulk",
                                 cellType1   = "lung.Foxp3p",
                                 cellType2   = "lung.nonFoxp3"),
                            list(libraryType = "bulk",
                                 cellType1   = "lung.Foxp3p",
                                 cellType2   = "lung.exFoxp3"),
                            list(libraryType = "bulk",
                                 cellType1   = "lung.nonFoxp3",
                                 cellType2   = "lung.exFoxp3")
                            )

   if (missing(libNames)) {
      manualLibNames <- list("S015_WT_lung_Treg_ova_set2" = "Treg.ova",
                             "S018_WT_lung_Treg_pbs_set2" = "Treg.pbs",
                             "WT_lung_Th_ova" = "Th.ova",
                             "WT_lung_Th_PBS" = "Th.pbs",
                             "WT_mLN_Th_PBS" = "mLN_Th.pbs")
   } else {
      manualLibNames <- libNames
   }
   dat <- list() ;

   if (missing(seuratFn)) {
      dat <- scoreSCsig.setPaths(
         baseDir      = "/data/davisfp/projects/exfoxp3",
         kallistoDir  = "results/RNAseq/kallisto/",
         sleuthDir    = "results/RNAseq/sleuth/",
         rnaSamplesFn = "metadata/exfoxp3_rna_samples.txt",
         scDataDir    = "run/20180708.reaggr_bulk_treg/trth2_sc201807" ,
         outDir       = "./"
      ) ;
   } else {
      dat <- scoreSCsig.setPaths(
         baseDir      = "/data/davisfp/projects/exfoxp3",
         kallistoDir  = "results/RNAseq/kallisto/",
         sleuthDir    = "results/RNAseq/sleuth/",
         rnaSamplesFn = "metadata/exfoxp3_rna_samples.txt",
         scDataFn     = seuratFn,
         outDir       = "./"
      ) ;
   }
   
   dat <- scoreSCsig.setSpecs(dat,
                              experimentName = experimentName,
                              deg.celltypePairs = deg.celltypePairs,
                              manualLibNames = manualLibNames) ;
   
   dat <- scoreSCsig.loadData(dat, returnData=TRUE) ;

   if (!missing(curRdsFn)) {
      saveRDS(dat,curRdsFn)
   }

   return(dat)

}


# Load 10X data into Seurat object
# run Seurat tSNE
# run Seurat clustering at several res; find markers
# DE genes: beween Treg assigned to different clusters
# DE genes: between active Treg in OVA v PBS
# DE genes: between active Th in OVA v PBS
# DE genes: between naive Th in lung v mLN

absc.demix.exFoxp3 <- function(dat,
                          geneList = c("Foxp3", "Ctla4", "Izumo1r", "Ikzf3", "eGFP")) {

# PURPOSE: Fancy, in progress, decontamination of exFoxp3 reporter profile
#          Implements linear model using STAN

   library(rstan)
   stanFn <- paste0(dat$specs$baseDir,"/src/STAN/decontaminate.stan")

# Purpose: call STAN routine to de-contaminate ex-Foxp3 reporter profile

   allCols <- colnames(dat$edat$geneExpr.clean)
   contamCols <- allCols[grep("tpm.*lung.*RFPp.*GFPp", allCols)]
   mixCols <- allCols[grep("tpm.*lung.*RFPp.*GFPn", allCols)]
   print(contamCols)
   print(mixCols)

   eMat <- dat$edat$geneExpr.clean
   eMat <- eMat[,c("gene_name",contamCols, mixCols)]
   eMat <- eMat %>% group_by(gene_name) %>% summarise_all(funs(sum))
   print("NOT SURE WHY: IGNORING NA gene name entry!!! *********** ****** ~~~~~~~ !!!!!!!!!!!!!")
   eMat <- na.omit(eMat)
   logEmat <- data.frame(log1p(eMat[,c(contamCols,mixCols)]))

   rownames(logEmat) <- eMat$gene_name
   logEmat$gene_name <- NULL

   colnames(logEmat) <- c(contamCols, mixCols)

   contamCols <- allCols[grep("tpm.*lung.*RFPp.*GFPp", allCols)]
   mixCols <- allCols[grep("tpm.*lung.*RFPp.*GFPn", allCols)]

   sm <- NA
   stanControlList <- list() #list(adapt_delta = 0.95)

   inData <- list(nMix = length(mixCols),
                  nContam = length(contamCols),
                  f_c = 0.3,
                  logO_mix = logEmat[geneList,mixCols],
                  logO_contam = logEmat[geneList,contamCols],
                  nGenes = length(geneList)
                  )

   stanFit <- rstan::stan(file = stanFn,
                           data = inData,
                           iter = 1000,
                           control = stanControlList,
                           chains = 4)

   return(list(stanFit = stanFit, logEmat = logEmat, geneList = geneList))

   if (0) {
   stanFits <- list()
   for (gene in geneList) {

      inData <- list(nMix = length(mixCols),
                     nContam = length(contamCols),
                     f_c = 0.3,
                     logO_mix = unlist(logEmat[gene,mixCols]),
                     logO_contam = unlist(logEmat[gene,contamCols]),
                     nGenes = 1
                     )

      if (is.na(sm)) {
         stanFits[[gene]] <- rstan::stan(file = stanFn,
                                         data = inData,
                                         iter = 1000,
                                         control = stanControlList,
                                         chains = 4)
         sm <- stanFits[[gene]]
      } else {
         stanFits[[gene]] <- rstan::stan(fit= sm,
                                         data = inData,
                                         iter = 1000,
                                         control = stanControlList,
                                         chains = 4)
      }

   }
   return(stanFits)
   }

}

absc.cleanDEG.exFoxp3 <- function(dat,
                                  tpmOnly = TRUE,
                                  f_c = 0.3          # contaminating Treg fraction
                                 ) {

# Purpose: Given FACS-estimated Treg fraction in Rfp+Gfp- ex-Foxp3 reporter,
#          subtract weighted Rfp+Gfp+ contribution from observed Rfp+Gfp- 
#          and then use this to filter out sleuth-defined ex-Foxp3 DEG that
#          are no longer > FC cutoff.


   eMat <- dat$edat$geneExpr.clean
   allCols <- colnames(dat$edat$geneExpr.clean)
   nonCols <- allCols[grep("tpm.*lung.*RFPn.*GFPn", allCols)]
   contamCols <- allCols[grep("tpm.*lung.*RFPp.*GFPp", allCols)]
   mixCols <- allCols[grep("tpm.*lung.*RFPp.*GFPn", allCols)]
   print(nonCols)
   print(mixCols)
   print(contamCols)

   eMat <- eMat[,c("gene_name",contamCols, mixCols,nonCols)]
   eMat <- eMat %>% group_by(gene_name) %>% summarise_all(funs(sum))
   print("NOT SURE WHY: IGNORING NA gene name entry!!! *********** ****** ~~~~~~~ !!!!!!!!!!!!!")
   eMat <- as.data.frame(na.omit(eMat))
   rownames(eMat) <- eMat$gene_name
   eMat$gene_name <- NULL

   print(contamCols)

   decontam.exFoxp3 <- ( apply(eMat[,mixCols],1,mean) -
                         f_c * apply(eMat[,contamCols],1,mean) ) /
                       ( 1 - f_c )

   newFC.non_vs_ex <- (apply(eMat[,nonCols],1,mean) + 1) / (decontam.exFoxp3 + 1)

   upGenes.non_vs_ex <- dat$deg$lung.nonFoxp3_vs_lung.exFoxp3$upGenes.minExpr
   downGenes.non_vs_ex <- dat$deg$lung.nonFoxp3_vs_lung.exFoxp3$downGenes.minExpr


   newMat <- cbind(exFoxp3 = decontam.exFoxp3,
                   nonFoxp3 = apply(eMat[,nonCols],1,mean))
# Remove genes that no longer meet minTPM or FC thresholds

   newUpGenes <- intersect(upGenes.non_vs_ex,
                           rownames(newMat)[newMat[,"nonFoxp3"] >=
                              dat$specs$thresh$deGenes.FC * newMat[,"exFoxp3"]])
                              print("GOT HERE")

   newDownGenes <- intersect(downGenes.non_vs_ex,
                           rownames(newMat)[newMat[,"exFoxp3"] >=
                              dat$specs$thresh$deGenes.FC * newMat[,"nonFoxp3"]])
   
   return(list(
               newUpGenes = newUpGenes,
               newDownGenes = newDownGenes,
               origUpGenes = upGenes.non_vs_ex,
               origDownGenes = downGenes.non_vs_ex))
}

rebalanceColorMap <- function( scoreRange, scoreMid = 0,
                               minCap, maxCap, colorScale,
                               nBins = 100) {

   scaleMin <- max(minCap,scoreRange[1])
   scaleMax <- min(maxCap,scoreRange[2])

# scaleMin is assigned colorScale[1]
# scaleMax gets colorScale[3]
# now define the midpoint -- which bin (of nBins) contains scoreMid

   midBin <- floor( nBins * (scoreMid - scaleMin) / (scaleMax - scaleMin))

   colors <- c( colorRampPalette(c(colorScale[1], colorScale[2]))(midBin),
                colorRampPalette(c(colorScale[2], colorScale[3]))(nBins - midBin - 1))

   return(list(colors=colors,scaleMin=scaleMin,scaleMax=scaleMax))
}


readDispiritoTregSig <- function(baseDir) {
# Purpose: read Treg signatures from DiSpirito et al., Sci Immunol 2018

   # Read signature name legend
   tabFn <- paste0(baseDir,
                   "/data/DiSpirito2018/DiSpirito2018_aat5861_Table_S1.xlsx")
   tx2 <- read_excel(tabFn, range="D2:D16")
   tx2 <- strsplit(as.data.frame(tx2)[,1], split=" = ")

   sigCodes2Name <- list()
   sigGenes <- list()
   for (i in 1:length(tx2)) {
      sigCode <- tx2[[i]][1]
      origName <- tx2[[i]][2]
      sigDir <- gsub(".* ", "", origName)
      sigName <- gsub(" down", "", origName)
      sigName <- gsub(" up", "", sigName)

      print(paste0("sigCode=",sigCode,"; sigDir=",sigDir,"; sigName=",sigName))

      sigCodes2Name[[sigCode]] <- c(sigName, sigDir)
      if (! sigName %in% names(sigGenes)) {
         sigGenes[[sigName]] <- list()
      }
      sigGenes[[sigName]][[sigDir]] <- c()
   }

   sigTab <- as.data.frame(read_excel(tabFn, range="A2:B1911"))
   for (i in 1:nrow(sigTab)) {
      sigName <- sigCodes2Name[[sigTab[i,2]]][1]
      sigDir <-sigCodes2Name[[sigTab[i,2]]][2]
      sigGenes[[sigName]][[sigDir]] <- c( sigGenes[[sigName]][[sigDir]],
                                          sigTab[i,1])
   }

   return(sigGenes)

}


scoreSCsig.loadData <- function(dat, returnData = FALSE, plotDEG=FALSE){
#   scoreSCsig.startup()

# 0. Load Data

   if (missing(dat))                    dat <- list()

   if (! "paths" %in% names(dat))       dat <- scoreSCsig.setPaths(dat)
   if (! "specs" %in% names(dat))       dat <- scoreSCsig.setSpecs(dat)

   if (! "sc" %in% names(dat)) {
      print("Loading single cell RNA-seq data")
      dat <- scoreSCsig.scidr.loadExpr(dat)
   }

   if (! "edat" %in% names(dat) &
       "kallistoBaseDir" %in% names(dat$paths)) {
      print("Loading bulk RNA-seq expression")
      dat$edat <- scoreSCsig.loadExpr(dat)
   }

   if (! "deg" %in% names(dat) &
       "sleuthBaseDir" %in% names(dat$paths)) {
      print("Loading bulk differentially expressed genes")
      dat <- scoreSCsig.defineDEG(dat, plotDEG=plotDEG)
   }

   return(dat)

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


scoreSCsig.scidr.loadExpr <- function(dat) {


   print("Loading SC expression")
   if (dat$specs$sc$dataType == "cellranger") {
      dat <- scoreSCsig.loadExprSC.cellranger(dat)
   } else if (dat$specs$sc$dataType == "seurat") {
      dat <- scoreSCsig.loadExprSC.seurat(dat)
   }

   return(dat)

}

scoreSCsig.loadExprSC.cellranger <- function( dat,
                                              loadMat         = TRUE,
                                              loadTsne        = TRUE,
                                              loadClusters    = TRUE) {

   dataDir <- dat$paths$sc.dataDir

   fn.invoc <- paste0(dataDir,"/_invocation")

   runMode <- "count"
   fn.mexString <- ""
   if (any(grepl('aggregation_csv', readLines(fn.invoc)))) {

       runMode <- "aggr"
       fn.mexString <- "_mex"

       fn.aggrCSV <- paste0(dataDir,"/SC_RNA_AGGREGATOR_CS/PARSE_CSV/fork0/files/aggregation_csv.csv")

       print(paste0("FN: ", fn.aggrCSV))
   }
   print(paste0("Loading cellranger ",runMode," output"))

   fn.eMat <- Sys.glob(paste0(dataDir,'/outs/filtered_gene_bc_matrices',fn.mexString,'/*/matrix.mtx'))
   fn.eMat.genes <- Sys.glob(paste0(dataDir,'/outs/filtered_gene_bc_matrices',fn.mexString,'/*/genes.tsv'))
   fn.eMat.cells <- Sys.glob(paste0(dataDir,'/outs/filtered_gene_bc_matrices',fn.mexString,'/*/barcodes.tsv'))

   if (missing(dat)) { dat <- list() }

   if (runMode == "aggr") {
      x <- read.csv(fn.aggrCSV,header=TRUE)
      dat$sc$fullLibInfo <- x
      if ("Group" %in% colnames(x)) {
         dat$sc$libNames <- as.character(x$Group)
      } else {
         dat$sc$libNames <- as.character(x$library_id)
      }
   }

   if ("manualLibNames" %in% names(dat$specs)) {
      print("Renaming libraries")
      oldNames <- dat$sc$libNames
      newNames <- oldNames
      for (origName in names(dat$specs$manualLibNames)) {
         newNames[dat$sc$libNames == origName] <- dat$specs$manualLibNames[[origName]]
      }
      dat$sc$libNames <- newNames
   }

   print("LIBRARIES:")
   print(dat$sc$libNames)


   skipCells <- c()
   if (loadMat) {
      print(paste0("Loading eppression matrix: ", fn.eMat))

      print(paste0("- read expression matrix: ",fn.eMat))
      dat$sc$eMat <- readMM(fn.eMat)
      print("-> DONE")

      print(paste0("- read gene names: ",fn.eMat.genes))
      dat$sc$geneInfo <- read.table(fn.eMat.genes, sep="\t",header=FALSE)
      colnames(dat$sc$geneInfo) <- c("geneID","geneName")
      print("-> DONE")

      print(paste0("- read cell names: ",fn.eMat.cells))
      cellNames <- read.table(fn.eMat.cells, sep="\t", header=FALSE)
      print("-> DONE")

      rownames(dat$sc$eMat) <- dat$sc$geneInfo$geneID
      colnames(dat$sc$eMat) <- cellNames$V1
      dat$sc$cell2lib <- gsub(".*-","", colnames(dat$sc$eMat))
      dat$sc$cell2lib <- dat$sc$libNames[as.numeric(dat$sc$cell2lib)]
      names(dat$sc$cell2lib) <- colnames(dat$sc$eMat)

      origCells <- names(dat$sc$cell2lib)
      if ("skip.sc.samples" %in% names(dat$specs)) {

         dat$sc$libNames <- setdiff(dat$sc$libNames, dat$specs$skip.sc.samples)
         skipCells <- names(dat$sc$cell2lib)[dat$sc$cell2lib %in% skipLibraries]

         print(paste0("Keeping ", length(dat$sc$libNames)," samples"))
         print(paste0("Skipping ", length(dat$specs$skip.sc.samples)," samples"))

      } else if ("keep.sc.samples" %in% names(dat$specs)) {

         keepLibraries <- dat$specs$keep.sc.samples
         skipLibraries <- setdiff(dat$sc$libNames, keepLibraries)
         skipCells <- names(dat$sc$cell2lib)[dat$sc$cell2lib %in% skipLibraries]
         dat$sc$libNames <- keepLibraries

         print(paste0("Keeping ", length(dat$sc$libNames)," samples"))
         print(paste0("Skipping ", length(skipLibraries)," samples"))
      }

      if (length(skipCells) > 0) {
         print(paste0("SKIPPING n=",length(skipCells)," cells"))

         print(paste0("orig eMat n=",ncol(dat$sc$eMat)," cells"))
         dat$sc$eMat <- dat$sc$eMat[,! colnames(dat$sc$eMat) %in% skipCells]

         print(paste0("->trim to n=",ncol(dat$sc$eMat)," cells"))
         dat$sc$cell2lib <- dat$sc$cell2lib[! names(dat$sc$cell2lib) %in% skipCells]
      }

   }


   if (loadClusters) {
      print("Loading clusters")
      clusterTypes <- list.files(paste0(dataDir,"/outs/analysis/clustering"))
      print(paste0(clusterTypes,collapse=", "))
      dat$sc$clusters <- list()
      for (clusterType in clusterTypes) {
         print(paste0("-> ",clusterType))
         fn.clusters <- paste0(dataDir,"/outs/analysis/clustering/",
                               clusterType,"/clusters.csv")

         x <- read.csv(fn.clusters, header=TRUE)
         x$Barcode <- as.character(x$Barcode)
         x$libraryID <- as.numeric(gsub(".*-","",x$Barcode))

         x <- x[!x$Barcode %in% skipCells,]

         nLib <- length(unique(x$libraryID))
         nCluster <- max(x$Cluster)

         if (runMode == "aggr") { x$libraryName <- dat$sc$libNames[x$libraryID]}

         dat$sc$clusters[[clusterType]] <- list(details = x,
                                                nCluster = nCluster,
                                                nLib = nLib)
      }

      if ("sc.clusterType" %in% names(dat$specs)) {
         dat$sc$clusterType <- dat$specs$sc.clusterType
      } else {
         dat$sc$clusterType <- "graphclust"
      }
   }

   if (loadTsne) {
      print("Loading tSNE results")
      tsneFn <- paste0(dataDir,"/outs/analysis/tsne/2_components/projection.csv")
      x <- read.csv(tsneFn, header=TRUE,stringsAsFactors=FALSE)
      dat$sc$tsneCoords <- x
      colnames(dat$sc$tsneCoords) <- c("Barcode","tsne1","tsne2")
      dat$sc$tsneCoords <- dat$sc$tsneCoords[! names(dat$sc$tsneCoords) %in%
                                             skipCells,]
   }

   return(dat)

}


scoreSCsig.loadExprSC.seurat <- function( dat,
                                          loadMat = TRUE,
                                          loadTsne = TRUE,
                                          loadClusters = TRUE) {

# PURPOSE: load Seurat output

   dataFn <- dat$paths$sc.dataFn

   print(paste0("Loading seurat object: ",dataFn))
   tempenv <- new.env()
   load(dataFn,envir=tempenv)
   print(paste0("- picking out right bit"))
   objName <- names(tempenv)[1]

# so = seurat object
   so <- tempenv[[names(tempenv)[1]]]

   if (missing(dat)) { dat <- list() }

   if ("manualLibNames" %in% names(dat$specs)) {
      print("Renaming libraries")
      for (libNum in names(dat$specs$manualLibNames)) {
         curName <- dat$specs$manualLibNames[[libNum]]
         print(paste0(" ",libNum," -> ", curName))
         so@meta.data[grep(paste0("-",libNum),
                           colnames(so@data)),"sample"] <- curName
      }
      so@meta.data$orig.ident <- as.factor(so@meta.data$sample)
   }

   dat$sc$libNames <- unique(so@meta.data$orig.ident)
   print(paste0("ALL LIBRARIES: ",dat$sc$libNames))

   skipCells <- c()
   if (loadMat) {
      print("Loading expression matrix")
      dat$sc$eMat <- so@data
      dat$sc$eMat.raw <- so@raw.data # adding eMat.raw to hold UMI data

      dat$sc$geneInfo <- read.table(dat$paths$sc.geneInfoFn,
                                    sep="\t", header=FALSE)
      colnames(dat$sc$geneInfo) <- c("geneID","geneName")

# convert rownames to geneID to match cellranger pipeline

      curGeneNames <- intersect(rownames(dat$sc$eMat), dat$sc$geneInfo$geneName)
      if (length(curGeneNames) < nrow(dat$sc$eMat)) {
         print(paste0("WARNING! skipped n=",
                      length(setdiff(rownames(dat$sc$eMat), curGeneNames)),
                      "genes in SC matrix that we don't have info for"))
      }
      dat$sc$eMat <- dat$sc$eMat[curGeneNames,]
      dat$sc$eMat.raw <- dat$sc$eMat.raw[curGeneNames,]

      curGeneID <- dat$sc$geneInfo$geneID[match(rownames(dat$sc$eMat),
                                                dat$sc$geneInfo$geneName)]

      rownames(dat$sc$eMat) <- curGeneID
      rownames(dat$sc$eMat.raw) <- curGeneID

      # rownames is gene name, not id.
      dat$sc$cell2lib <- so@meta.data$orig.ident
      names(dat$sc$cell2lib) <- colnames(dat$sc$eMat)

      origCells <- names(dat$sc$cell2lib)
      if ("skip.sc.samples" %in% names(dat$specs)) {

         dat$sc$libNames <- setdiff(dat$sc$libNames, dat$specs$skip.sc.samples)
         skipCells <- names(dat$sc$cell2lib)[dat$sc$cell2lib %in% skipLibraries]

         print(paste0("Keeping ", length(dat$sc$libNames)," samples"))
         print(paste0("Skipping ", length(dat$specs$skip.sc.samples)," samples"))

      } else if ("keep.sc.samples" %in% names(dat$specs)) {

         keepLibraries <- dat$specs$keep.sc.samples
         skipLibraries <- setdiff(dat$sc$libNames, keepLibraries)
         skipCells <- names(dat$sc$cell2lib)[dat$sc$cell2lib %in% skipLibraries]
         dat$sc$libNames <- keepLibraries

         print(paste0("Keeping ", length(dat$sc$libNames)," samples"))
         print(paste0("Skipping ", length(skipLibraries)," samples"))

      }

      if (length(skipCells) > 0) {
         print(paste0("SKIPPING n=",length(skipCells)," cells"))

         print(paste0("orig eMat n=",ncol(dat$sc$eMat)," cells"))
         dat$sc$eMat <- dat$sc$eMat[,! colnames(dat$sc$eMat) %in% skipCells]
         dat$sc$eMat.raw <- dat$sc$eMat.raw[,! colnames(dat$sc$eMat.raw) %in% skipCells]

         print(paste0("->trim to n=",ncol(dat$sc$eMat)," cells"))
         dat$sc$cell2lib <- dat$sc$cell2lib[! names(dat$sc$cell2lib) %in% skipCells]
      }
   }


   if (loadClusters) {
      print("Loading clusters")

      clusterTypes <- names(so@meta.data)[grepl("^res", names(so@meta.data))]
      print(paste0(clusterTypes,collapse=", "))
      dat$sc$clusters <- list()
      for (clusterType in clusterTypes) {
         print(paste0("-> ",clusterType))

         x <- data.frame(Barcode        = colnames(dat$sc$eMat),
                         Cluster        = as.numeric(so@meta.data[[clusterType]][!origCells %in% skipCells]) + 1,
                         libraryName    = so@meta.data$orig.ident[!origCells %in% skipCells],
                         libraryID      = so@meta.data$orig.ident[!origCells %in% skipCells],
                         stringsAsFactors = FALSE)

         nLib <- length(unique(x$libraryID))
         nCluster <- max(x$Cluster)

         dat$sc$clusters[[clusterType]] <- list(details = x,
                                                nCluster = nCluster,
                                                nLib = nLib)
      }
      if ("sc.clusterType" %in% names(dat$specs)) {
         dat$sc$clusterType <- dat$specs$sc.clusterType
      } else {
         dat$sc$clusterType <- names(dat$sc$clusters)[1]
      }
   }

   if (loadTsne & "tsne" %in% names(so@dr)) {
      print("Loading tSNE results")
      dat$sc$tsneCoords <- data.frame(so@dr$tsne@cell.embeddings)
      colnames(dat$sc$tsneCoords) <- c("tsne1", "tsne2")
      dat$sc$tsneCoords$Barcode <- rownames(dat$sc$tsneCoords)
      dat$sc$tsneCoords <- dat$sc$tsneCoords[,c("Barcode","tsne1","tsne2")]

      dat$sc$tsneCoords <- dat$sc$tsneCoords[! dat$sc$tsneCoords$Barcode %in% 
                                             skipCells,]
   }

   dat$sc$origSO <- so

   return(dat)

}


scoreSCsig.defineDEG <- function(dat,
                                 plotDEG=TRUE ){

# consistent DEG pairwise SLEUTH comparison

# loaded in startup()
#   library(sleuth)
#   library(digest)

   genes.deg <- list()
   t2g <- data.frame(
      target_id = dat$edat$txExpr$orig_transcript_id,
      ens_gene  = dat$edat$txExpr$gene_id,
      ext_gene  = dat$edat$txExpr$gene_name,
      stringsAsFactors = FALSE
   )

#   t2g <- data.frame(
#      target_id = dat$specs$transcriptInfo$transcript_id,
#      ens_gene  = dat$specs$transcriptInfo$gene_id,
#      ext_gene  = dat$specs$transcriptInfo$gene_name,
#      stringsAsFactors = FALSE
#   )

   allSamples <- dat$specs$rnaSamples

   for (cellPair in dat$specs$deg.celltypePairs) {

      print(paste0("comparing ",cellPair$cellType1," to ",cellPair$cellType2))
      cell12 <- paste0(cellPair$cellType1,"_vs_",cellPair$cellType2)
      cell1 <- cellPair$cellType1
      cell2 <- cellPair$cellType2

      samples1   <- allSamples$sampleName[allSamples$cellType == cellPair$cellType1 &
                                          allSamples$libraryType == cellPair$libraryType]

      samples2   <- allSamples$sampleName[allSamples$cellType == cellPair$cellType2 &
                                          allSamples$libraryType== cellPair$libraryType]

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
         
      curSampleInfo$path <- paste0(dat$paths$kallistoBaseDir, "/",
                                      curSampleInfo$runID, "/",
                                      curSampleInfo$sampleName)

# sleuth info https://rawgit.com/pachterlab/sleuth/master/inst/doc/intro.html


      comparisonString <- paste0( paste0(sort(samples1), collapse = "_"),
                                  "_vs_",
                                  paste0(sort(samples2), collapse = "_") )

      comparisonMD5 <- digest(comparisonString, algo = "md5")
      curSleuthPath <- paste0(dat$paths$sleuthBaseDir, "/", comparisonMD5)

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

         sleuthData <- scoreSCsig.readSleuthOutput(resultsDir = curSleuthPath)
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

         scoreSCsig.writeSleuthOutput(
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

scoreSCsig.readSleuthOutput <- function(comparisonName, resultsDir){

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



scoreSCsig.writeSleuthOutput <- function(comparisonName,
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
