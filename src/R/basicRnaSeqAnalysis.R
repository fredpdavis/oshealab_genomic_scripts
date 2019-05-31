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
      print("Plot: DEG scatterplot")
      plotDEGscatters(dat)
   }

   if (any(c("all","4") %in% figList)) {
      print("Plot: fold change cumulative distribution plots")
      plotSmale(dat)
   }

# GEO supp expression table
   if (any(c("all","GEO_table") %in% figList)) {
      print("Table: Gene expresion table")
      makeGeoExpressionTable(dat)
   }

# Tables of DEG genes
   if (any(c("all","DEG_table") %in% figList)) {
      print("Table: DEG table")
      writeDEG(dat)
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
         gene_id        = dat$edat$geneExpr$gene_id,
         gene_name      = dat$edat$geneExpr$gene_name,
         tpm.celltype2  = apply(dat$edat$geneExpr[, paste0("tpm.", samples2)],
                                1, mean),
         tpm.celltype1  = apply(dat$edat$geneExpr[, paste0("tpm.", samples1)],
                                1, mean),
         stringsAsFactors = FALSE
      )
      comparisonFC$log2fc <- log2(1 + comparisonFC$tpm.celltype2) -
                             log2(1 + comparisonFC$tpm.celltype1)


      upGenes <- comparisonFC$gene_id[
         comparisonFC$tpm.celltype2 >= comparisonFC$tpm.celltype1 *
             dat$specs$thresh$deGenes.FC]

      print("head(upGenes)")
      print(head(upGenes))
   
      downGenes <- comparisonFC$gene_id[
         comparisonFC$tpm.celltype1 >= comparisonFC$tpm.celltype2 *
             dat$specs$thresh$deGenes.FC]
  
      print("head(downGenes)")
      print(head(downGenes))

      exprGenes.up <- intersect(upGenes, 
                                comparisonFC$gene_id[
         comparisonFC$tpm.celltype2 >= dat$specs$thresh$exprGenes.minTPM])

      print("head(exprGenes.up)")
      print(head(exprGenes.up))

      exprGenes.down <- intersect(downGenes, 
                                comparisonFC$gene_id[
         comparisonFC$tpm.celltype1 >= dat$specs$thresh$exprGenes.minTPM])
      print("head(exprGenes.down)")
      print(head(exprGenes.down))
  

      if (file.exists(paste0(curSleuthPath, "/done.txt"))) {
   
         sleuthData <- readSleuthOutput(resultsDir = curSleuthPath)
         genes.deg[[cell12]]$degGenes           <- sleuthData$degGenes
         genes.deg[[cell12]]$upGenes            <- sleuthData$upGenes
         genes.deg[[cell12]]$downGenes          <- sleuthData$downGenes
         genes.deg[[cell12]]$sleuthResults.gene <- sleuthData$sleuthResults.gene
         genes.deg[[cell12]]$sleuthResults.tx   <- sleuthData$sleuthResults.tx
         genes.deg[[cell12]]$designMat          <- sleuthData$designMat
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
                                               aggregation_column       = 'ext_gene',
                                               extra_bootstrap_summary  = TRUE,
                                               target_mapping = t2g)
      
         print(paste("sleuth_fit()", date()))
         genes.deg[[cell12]]$so <- sleuth_fit(genes.deg[[cell12]]$so)
      
         print(paste("sleuth_wt()", date()))
         genes.deg[[cell12]]$so <- sleuth_wt(genes.deg[[cell12]]$so,
                                            paste0('condition',cond2))
      
         print(paste("sleuth_results() gene level", date()))

         genes.deg[[cell12]]$sleuthResults.gene <- sleuth_results(
            genes.deg[[cell12]]$so, paste0('condition',cond2),
            pval_aggregate=TRUE)

         print("head(genes.deg[[cell12]]$sleuthResults.gene)")
         print(head(genes.deg[[cell12]]$sleuthResults.gene))


         print(paste("sleuth_results() tx level", date()))

         genes.deg[[cell12]]$sleuthResults.tx <- sleuth_results(
            genes.deg[[cell12]]$so, paste0('condition',cond2),
            pval_aggregate=FALSE)

         print(paste("DONE", date()))
   
         genes.deg[[cell12]]$degGenes <- unique(na.omit(
           genes.deg[[cell12]]$sleuthResults.gene$ens_gene[
               genes.deg[[cell12]]$sleuthResults.gene$qval <=
               dat$specs$thresh$sleuth.qval]))

         print("head(genes.deg[[cell12]]$degGenes)")
         print(head(genes.deg[[cell12]]$degGenes))

         genes.deg[[cell12]]$upGenes <- intersect(genes.deg[[cell12]]$degGenes,
                                                  upGenes)

         genes.deg[[cell12]]$downGenes <- intersect(genes.deg[[cell12]]$degGenes,
                                                  downGenes)

         writeSleuthOutput(
            sleuthResults.gene  = genes.deg[[cell12]]$sleuthResults.gene,
            sleuthResults.tx    = genes.deg[[cell12]]$sleuthResults.tx,
            degGenes            = genes.deg[[cell12]]$degGenes,
            upGenes             = genes.deg[[cell12]]$upGenes,
            downGenes           = genes.deg[[cell12]]$downGenes,
            designMat           = genes.deg[[cell12]]$designMat,
            outDir              = curSleuthPath,
            comparisonName      = comparisonString )
     }
   
      genes.deg[[cell12]]$fc.deg <- comparisonFC[ comparisonFC$gene_id %in%
         unique(c(genes.deg[[cell12]]$upGenes,
                  genes.deg[[cell12]]$downGenes)),]


      genes.deg[[cell12]]$fc.all <- comparisonFC

      genes.deg[[cell12]]$fc.all <- merge(
         genes.deg[[cell12]]$fc.all,
         genes.deg[[cell12]]$sleuthResults.gene,
         by.x="gene_id",
#         by.x="gene_name",
         by.y="target_id",
         all.x=TRUE,all.y=FALSE)

      genes.deg[[cell12]]$fc.all$pval[is.na(genes.deg[[cell12]]$fc.all$pval)] <- 1
      genes.deg[[cell12]]$fc.all$qval[is.na(genes.deg[[cell12]]$fc.all$qval)] <- 1
      genes.deg[[cell12]]$fc.all$pi <- genes.deg[[cell12]]$fc.all$log2fc *
                                       -1 * log10(genes.deg[[cell12]]$fc.all$qval)
   
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

   specs$ensembl_release<- "94"

   specs$genomeVer <- "GRCm38"
   specs$txVer <- paste0(specs$genomeVer,".",specs$ensembl_release)

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
      "/data/kallisto_files.",specs$txVer,"/transcript_info.",specs$txVer,".txt"),
      quote = "", header = TRUE, sep = "\t", as.is = TRUE)
#   specs$geneInfo <- unique(specs$transcriptInfo[, c("gene_id", "gene_name")])
   specs$geneInfo <- specs$transcriptInfo %>% group_by(gene_id,gene_name,chr) %>% summarize(minStart = min(start), maxEnd = max(end))



   specs$kallistoBaseDir <- paste0(specs$baseDir,
      "/results/RNAseq/kallisto.",specs$txVer,"/")

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
                            sampleNames, #sampleName to include in plot?
                            height = 5,
                            width = 6,
                            figName = "corrHeatmap") {

   if (missing(sampleNames)) {
      sampleNames <- dat$specs$rnaSamples$sampleName
   }

   eMat <- dat$edat$geneExpr[,paste0("tpm.",sampleNames)]

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
      print("TEST: SampleAnnotates")
      print(colAnn)
   }

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
       height = height ,
       width  = width)
   do.call(pheatmap, pheatmap.options)
   dev.off()
}


makeVarGeneHeatmap <- function(dat,
                            sampleAnnotate, #properties to show as sample annotation
                            sampleNames, #sampleName to include in plot?
                            pheatmap.gaps_col = NULL,
                            nlMode = "logMean", #"minMax", "fracMax"
                            height = 5,
                            width = 6,
                            justNumbers = FALSE,
                            pheatmap.cluster_rows = TRUE,
                            pheatmap.cluster_cols = TRUE,
                            figName = "variableGeneHeatmap") {

   if (missing(sampleNames)) {
      sampleNames <- dat$specs$rnaSamples$sampleName
   }

   eMat <- dat$edat$geneExpr[,paste0("tpm.",sampleNames)]
   eMat <- data.frame(eMat)
   rownames(eMat) <- paste0(dat$edat$geneExpr$gene_name,":",
                            dat$edat$geneExpr$gene_id)

   eMat <- eMat[apply(eMat,1,max) >= dat$specs$thresh$exprGenes.minTPM,]
   eMat <- eMat[apply(eMat,1,max) >= dat$specs$thresh$deGenes.FC *
                                     apply(eMat,1,min),]

   if (justNumbers) {
      print(paste0("NUMBER OF VARIABLE GENES IN HEATMAP:",nrow(eMat)))
      return(1) ;
   }


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
            cluster_rows = pheatmap.cluster_rows,
            cluster_cols = pheatmap.cluster_cols,
            gaps_col = pheatmap.gaps_col,
            treeheight_row = 0,
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
       height = height,
       width  = width)
   do.call(pheatmap, pheatmap.options)
   dev.off()
}


plotHmap.deg <- function(dat,
                         nlMode = "logMean",
                         geneList, #used if specified; otherwise picks DEGs
                         geneList.type = "gene_name", # or gene_id
                         degTypes, #if not specified, uses all DEG's
                         sampleList, #used if specified; otherwise uses DEG samples
                         sampleAnnotate, #properties to show as sample annotation
                         overlayTPM = FALSE,
                         show_rownames = TRUE,
                         clusterCols = FALSE,
                         clusterRows = TRUE,
                         rowFontSize = 2,
                         colFontSize = 6,
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
      geneList.type <- "gene_id"
      if (missing(degTypes)) {
         degTypes <- names(dat$specs$degSamples)
      }
      geneList <- c()
      for (degType in degTypes) {
         geneList <- c(geneList, dat$deg[[degType]]$upGenes.minExpr,
                                 dat$deg[[degType]]$downGenes.minExpr)
      }
      geneList <- unique(geneList)
   } else {
      if (geneList.type == "gene_name") {
         geneList.type <- "gene_id"
         geneList.orig <- geneList
         geneList <- c()
         for (gene in geneList.orig) {
            geneList <- c(geneList,
                          dat$specs$geneInfo$gene_id[dat$specs$geneInfo$gene_name == gene])
         }
      }
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
      match(geneList[geneList %in% origExpr$gene_id], origExpr$gene_id),
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
   origTPMmat <- hMat

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
   origTPMmat <- origTPMmat[rownames(hMat),]
   origTPMmat <- round(origTPMmat,digits=0)

   plotTitle <- nlMode
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
            fontsize_col   = colFontSize,
            fontsize_row   = rowFontSize,
            main = plotTitle)

   if (overlayTPM) {
      pheatmap.options$display_numbers <- origTPMmat
      pheatmap.options$main <- paste0("color=",pheatmap.options$main,
                                      "; #=TPM")
   }

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

   so <- NA
   soFn <- paste0(resultsDir,"/so.rds")
   if (file.exists(soFn)) { so <- readRDS(file=soFn) }

   sleuthResults.gene   <- readRDS(file=paste0(resultsDir,
                                               "/sleuthResults.gene.rds"))

   sleuthResults.tx     <- readRDS(file=paste0(resultsDir,
                                               "/sleuthResults.tx.rds"))

   designMat            <- readRDS(file=paste0(resultsDir,
                                               "/designMat.rds"))

   upGenes              <- scan(file=paste0(resultsDir,"/upGenes.txt"),
                                what="character")

   downGenes            <- scan(file=paste0(resultsDir,"/downGenes.txt"),
                                what="character")

   degGenes             <- scan(file=paste0(resultsDir,"/degGenes.txt"),
                                what="character")

   return(list(
      comparisonName     = comparisonName,
      so                 = so,
      designMat          = designMat,
      sleuthResults.gene = sleuthResults.gene,
      sleuthResults.tx   = sleuthResults.tx,
      degGenes           = degGenes,
      upGenes            = upGenes,
      downGenes          = downGenes
   ))

}



writeSleuthOutput <- function(comparisonName,
                              so=NULL,
                              sleuthResults.gene,
                              sleuthResults.tx,
                              designMat,
                              upGenes=NULL,
                              downGenes=NULL,
                              degGenes=NULL,
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

   saveRDS(sleuthResults.gene, file=paste0(outDir,"/sleuthResults.gene.rds"))
   saveRDS(sleuthResults.tx, file=paste0(outDir,"/sleuthResults.tx.rds"))
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

   if (!(is.null(degGenes))){
      upFn<-paste0(outDir,"/degGenes.txt")
      write.table(degGenes, file=upFn, quote=FALSE, row.names=FALSE,col.names=FALSE)
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


plotSmale <- function(dat, pseudocount = 1) {

   for (i in 1:nrow(dat$specs$rnaComps)) {

      cell1 <- dat$specs$rnaComps[i,"group1name"]
      cell2 <- dat$specs$rnaComps[i,"group2name"]
      cell12 <- paste0(cell1,"_vs_",cell2)

      if (grepl("unstim",cell2)){next;}

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
      comparisonFC$fc1v2 <- (pseudocount + comparisonFC$tpm.celltype1) /
                            (pseudocount + comparisonFC$tpm.celltype2)


      degGenes <- list()
      degGenes$down <- comparisonFC[comparisonFC$gene_name %in%
                                    dat$deg[[cell12]]$upGenes.minExpr,]
   
      degGenes$up    <- comparisonFC[comparisonFC$gene_name %in%
                                     dat$deg[[cell12]]$downGenes.minExpr,]

      for (direction in c("up","down")) {
   
         curMat <- degGenes[[direction]]

         hlColor <- "steelBlue2"
         labDir <- "repressed"
         if (direction == "down") {
            curMat$fc1v2 <- 1 / curMat$fc1v2
            hlColor <- "darkOrange2"
            labDir <- "induced"
         }
#         curMat$rank <- rank(-1 * curMat$fc1v2)
         curMat$rank <- rank(curMat$fc1v2)
         curMat$rank <- 100 * curMat$rank / length(curMat$rank)

         fcRange <- range(curMat$fc1v2)
         fcRange[1] <- 1


         tmppngfn <- tempfile()
         png(file = tmppngfn, height = 3.1, width = 3.1, units = "in", res = 300,
             family = "ArialMT")
         par(mar = c(0, 0,0, 0))
         plot(curMat$rank,
              curMat$fc1v2,
              ann=FALSE,axes=FALSE,
              pch=20, cex=0.5,
              log="y", col="grey",las=1,
              xlab = paste0(cell2," (TPM + 1)"),
              ylab = paste0(cell1," (TPM + 1)"),
              ylim=fcRange,
              main="")
         dev.off()
         pngbg <- readPNG(tmppngfn)
         pngbg <- as.raster(pngbg)
   
         outFn <- paste0(dat$specs$outDir,"/smalePlot.",cell12,".",direction,
                         ".pseudcount",pseudocount,".pdf")
         pdf(outFn, height=3.5,width=3.5)
      
         par(mar = c(3.75, 3.75, 0.5, 0.5),
             mgp = c(2, 0.6, 0),
             cex = 1, cex.axis = 0.8, cex.lab = 1)
      
         plot(curMat$rank,
              curMat$fc1v2,
              pch = 20,
              cex = 0.5,
              log = "y",
              type= "n",
#              axes=FALSE,
              las = 1,
              xlab = paste0("Percentile of ",labDir," genes"),
              ylab = paste0(cell1," vs ",cell2," (FC)"),
              ylim=fcRange,
              main="")
         lim<-par()
         rasterImage(pngbg, lim$usr[1], 10^lim$usr[3],
                            lim$usr[2], 10^lim$usr[4])
     
#         axis(1, at=-1 * seq(0,100,20), label=seq(0,100,20))
#         axis(2,las=1)
      
#         outliers <- curMat[curMat$rank <= 10,]
         outliers <- head(curMat[order(curMat$rank,decreasing=TRUE),],n=10)
   
         if (1) {
   
         points(outliers$rank,
                outliers$fc1v2,
                pch=20,
                cex=0.5,
                col= hlColor)
      
         lab_x <- c(outliers$rank)
         lab_y <- c(outliers$fc1v2)
         
         lab_text <- c(outliers$gene_name)
         lab_textadj <- c(rep(1, nrow(outliers)))
#         lab_textx <- c(rep(nrow(curMat)/2, nrow(outliers)))
         lab_textx <- c(rep(50, nrow(outliers)))
         laborder <- order(lab_y)
         
#         texty_logo <- ( 0.15 * log(fcRange[2], base = 10))
         texty_logo <- (log(fcRange[1], base = 10))
         texty_logo <- (0.15 * log(fcRange[2], base = 10))
         texty_loginc <- ((0.75 * log(fcRange[2], base = 10)) / length(lab_y))

#         return(fcRange = fcRange,
#                lab_x = lab_x,
#                lab_y = lab_y,
#                outliers=outliers)

      
         for (j in 1:length(laborder)) {
            k <- laborder[j]
            lab_texty <- (10**( (j - 1) * texty_loginc + texty_logo))
            text(lab_textx[k], lab_texty, lab_text[k], cex = 0.5, adj = lab_textadj[k])
            segments(lab_textx[k], lab_texty, lab_x[k], lab_y[k], lwd = 1, col = "grey")
         }
         legend("topleft",
                paste0("n=",nrow(curMat)," genes"),
                text.col=hlColor,
                cex=0.75, bty="n")
         }
      
      
         dev.off()
      }
   
   }

}




generic.plotDEGscatters <- function(exprMat, # data frame with expression levels
                                    outFn = "scatter.pdf",
                                    upGenes = c(),   # points to color above diagonal
                                    downGenes = c(), # points to color below diagonal
                                    labeledUpGenes = c(), # genes to label above diagonal
                                    labeledDownGenes = c(), # genes to label below diagonal
                                    color.upGenes = "darkOrange2",
                                    color.downGenes = "steelBlue2",
                                    colname.tpm1, #if not specified, assumes column 1
                                    colname.tpm2, #if not specified, assumes column 2
                                    colname.gene, #if not specified, assumes row name=gene name
                                    label.tpm1 = "tpm1", #condition 1 name
                                    label.tpm2 = "tpm2", #condition 2 name
                                    exprUnits = "TPM",
                                    axisMax #optional max TPM for plot
                                   ) {
# GOAL: Generic drop-in routine for others non-YARP code


   if (missing(colname.gene)) {
      exprMat$added.gene.name <- rownames(exprMat)
      colname.gene <- "added.gene.name"
   }

   if (missing(colname.tpm1)) { colname.tpm1 <- 1 }
   if (missing(colname.tpm2)) { colname.tpm2 <- 2 }

   if (missing(axisMax)) {
      axisMax <- 1 + max(exprMat[[colname.tpm1]],
                         exprMat[[colname.tpm2]])
   }
   exprRange <- c(1,axisMax)

   tmppngfn <- tempfile()
   png(file = tmppngfn, height = 3.1, width = 3.1, units = "in", res = 300,
       family = "ArialMT")
   par(mar = c(0, 0,0, 0))
   plot(1 + exprMat[[colname.tpm2]],
        1 + exprMat[[colname.tpm1]],
        ann=FALSE,axes=FALSE,
        pch=20, cex=0.5, log="xy", col="grey",las=1,
        xlab = paste0(label.tpm2," (",exprUnits," + 1)"),
        ylab = paste0(label.tpm1," (",exprUnits," + 1)"),
        xlim=exprRange,
        ylim=exprRange,
        main="")
   dev.off()
   pngbg <- readPNG(tmppngfn)
   pngbg <- as.raster(pngbg)
   
   pdf(outFn, height=3.5,width=3.5)

   par(mar = c(3.75, 3.75, 1, 1),
       mgp = c(2, 0.6, 0),
       cex = 1, cex.axis = 0.8, cex.lab = 1)

   plot(1 + exprMat[[colname.tpm2]],
        1 + exprMat[[colname.tpm1]],
        pch=20,
        cex=0.5,
        log="xy",
        type="n",
        las=1,
        xlab = paste0(label.tpm2," (",exprUnits," + 1)"),
        ylab = paste0(label.tpm1," (",exprUnits," + 1)"),
        xlim=exprRange,
        ylim=exprRange,
        main="")
   lim<-par()
   rasterImage(pngbg, 10^lim$usr[1], 10^lim$usr[3],
                      10^lim$usr[2], 10^lim$usr[4])

# Color upGenes and downGenes points
   degGenes.up <- exprMat[exprMat[[colname.gene]] %in% upGenes,]
   degGenes.down <- exprMat[exprMat[[colname.gene]] %in% downGenes,]

   if (nrow(degGenes.up) > 0) {
      points(1 + degGenes.up[[colname.tpm2]],
             1 + degGenes.up[[colname.tpm1]],
             pch=20,
             cex=0.5,
             col= color.upGenes)
   }

   if (nrow(degGenes.down) > 0) {
      points(1 + degGenes.down[[colname.tpm2]],
             1 + degGenes.down[[colname.tpm1]],
             pch=20,
             cex=0.5,
             col= color.downGenes)
   }

# Label labeledUpGenes and labeledDownGenes points
   outliers.down <- exprMat[exprMat[[colname.gene]] %in% labeledUpGenes,]
   outliers.up <- exprMat[exprMat[[colname.gene]] %in% labeledDownGenes,]

   if (nrow(outliers.down) + nrow(outliers.up) > 0) {

   lab_x <- c(1 + outliers.down[[colname.tpm2]],
              1 + outliers.up[[colname.tpm2]])

   lab_y <- c(1 + outliers.down[[colname.tpm1]],
              1 + outliers.up[[colname.tpm1]])

   rightLabelX <- exp(0.7 * log(axisMax))
   
   lab_text <- c(outliers.down[[colname.gene]], outliers.up[[colname.gene]])
   lab_textadj <- c(rep(1, nrow(outliers.down)),
                    rep(0, nrow(outliers.up)))
   lab_textx <- c(rep(2.5, nrow(outliers.down)),
                  rep(rightLabelX, nrow(outliers.up)))
   laborder <- order(lab_y)
   
   texty_logo <- ( 0.15 * log(axisMax, base = 10))
   texty_loginc <- ((0.75 * log(axisMax, base = 10)) / length(lab_y))

   for (j in 1:length(laborder)) {
      k <- laborder[j]
      lab_texty <- (10**( (j - 1) * texty_loginc + texty_logo))
      text(lab_textx[k], lab_texty, lab_text[k], cex = 0.5, adj = lab_textadj[k])
      segments(lab_textx[k], lab_texty, lab_x[k], lab_y[k], lwd = 1, col = "grey")
   }
   }



# Add legends for number of upGenes and downGenes

   legend("topleft",
          paste0("n=",nrow(degGenes.up)," genes"),
          text.col=color.upGenes,
          cex=0.75, bty="n")
   legend("bottomright",
          paste0("n=",nrow(degGenes.down)),
          text.col=color.downGenes,
          cex=0.75, bty="n")


   dev.off()
   
}


plotDEGscatters <- function(dat, scatterPlot=TRUE, nLabel=10,
                            colorUp = "darkOrange2",colorDown = "steelBlue2") {

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
         gene_id        = dat$edat$geneExpr$gene_id,
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
      plot(1 + comparisonFC$tpm.celltype1,
           1 + comparisonFC$tpm.celltype2,
           ann=FALSE,axes=FALSE,
           pch=20, cex=0.5, log="xy", col="grey",las=1,
           xlab = paste0(cell1," (TPM + 1)"),
           ylab = paste0(cell2," (TPM + 1)"),
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
   
      plot(1 + comparisonFC$tpm.celltype1,
           1 + comparisonFC$tpm.celltype2,
           pch=20,
           cex=0.5,
           log="xy",
           type="n",
           las=1,
           xlab = paste0(cell1," (TPM + 1)"),
           ylab = paste0(cell2," (TPM + 1)"),
           xlim=exprRange,
           ylim=exprRange,
           main="")
      lim<-par()
      rasterImage(pngbg, 10^lim$usr[1], 10^lim$usr[3],
                         10^lim$usr[2], 10^lim$usr[4])
  
#HERNEOW 190429_1642 - MAJORLY!  should keep DEGs in gene-id form. since multiple gene-id with same gene-name, fucking up DEG scatterpot.
      degGenes.up <- comparisonFC[comparisonFC$gene_id %in%
                                  dat$deg[[cell12]]$upGenes.minExpr,]
   
      degGenes.down  <- comparisonFC[comparisonFC$gene_id %in%
                                     dat$deg[[cell12]]$downGenes.minExpr,]
                              
      points(1 + degGenes.up$tpm.celltype1,
             1 + degGenes.up$tpm.celltype2,
             pch=20,
             cex=0.5,
             col= colorUp)
   
      points(1 + degGenes.down$tpm.celltype1,
             1 + degGenes.down$tpm.celltype2,
             pch=20,
             cex=0.5,
             col=colorDown)
   
   
      outliers.up   <- head(degGenes.up[order(degGenes.up$log2fc.celltype1_vs_celltype2,decreasing=FALSE),],
                            n = nLabel)
      outliers.down <- head(degGenes.down[order(degGenes.down$log2fc.celltype1_vs_celltype2,decreasing=TRUE),],
                          n = nLabel)
   
      lab_x <- c(1 + outliers.up$tpm.celltype1, 1 + outliers.down$tpm.celltype1)
      lab_y <- c(1 + outliers.up$tpm.celltype2, 1 + outliers.down$tpm.celltype2)
      
      lab_text <- c(outliers.up$gene_name, outliers.down$gene_name)
      lab_textadj <- c(rep(1, nrow(outliers.up)), rep(0, nrow(outliers.down)))
      lab_textx <- c(rep(2.5, nrow(outliers.up)), rep(6500, nrow(outliers.down)))
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
             text.col=colorUp,
             cex=0.75, bty="n")
      legend("bottomright",
             paste0("n=",nrow(degGenes.down)),
             text.col=colorDown,
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

runClusterProfiler <- function(dat) {

#   devtools::install_github(c("GuangchuangYu/clusterProfiler"))
#   install.packages(c("misgdbr"))


   library(clusterProfiler)
   library(enrichplot)
   library(ReactomePA)
   library(msigdbr)
   library(org.Mm.eg.db)
   library(plyr)
   library(dplyr)
   library(data.table)

   for (i in 1:nrow(dat$specs$rnaComps)) {
      cell1 <- dat$specs$rnaComps[i,"group1name"]
      cell2 <- dat$specs$rnaComps[i,"group2name"]
      cell12 <- paste0(cell1,"_vs_",cell2)
      if (cell12 == "WT.unstim_vs_Stat4KO.unstim") next;
      print(paste0("NOW ON ",cell12,"----------------------"))

      tx <- runClusterProfiler.meat(
               degGenes.up      = dat$deg[[cell12]]$upGenes.minExpr,
               degGenes.down    = dat$deg[[cell12]]$downGenes.minExpr,
               fullGeneInfo     = dat$deg[[cell12]]$fc.all,
               name             = cell12,
               outDir           = dat$specs$outDir
            )
   }

}

runClusterProfiler.meat <- function(degGenes.up,
                                    degGenes.down,
                                    fullGeneInfo,
                                    name,
                                    outDir) {

# Parameters
   nShow      <- 20

# Input: DEG list

   degGenes <- c(degGenes.up, degGenes.down)
   degDF      <- data.frame(SYMBOL = degGenes)

# Add additional gene identifiers (names, entrezid, KEGG codes, etc)
   tx.entrez  <- bitr(degDF$SYMBOL,
                      fromType="SYMBOL",
                      toType="ENTREZID",
                      OrgDb="org.Mm.eg.db")

   tx.kegg    <- bitr_kegg(tx.entrez$ENTREZID,
                           fromType='ncbi-geneid',
                           toType='kegg',
                           organism='mmu')

   tx.kegg    <- rename(tx.kegg, c("kegg" = "KEGG",
                                   "ncbi-geneid" = "ENTREZID"))

   degDF      <- merge(degDF, tx.entrez, by="SYMBOL")
   degDF      <- merge(degDF, tx.kegg, by="ENTREZID")
   degDF      <- degDF[,c("SYMBOL","ENTREZID","KEGG")]
   tx.entrez <- NULL; tx.kegg <- NULL


# Prep ranked gene data frame for GSEA analysis
   gseaDF     <- data.frame(SYMBOL = fullGeneInfo$gene_name,
                            GSEArank = fullGeneInfo$pi)
   gseaDF     <- gseaDF[order(gseaDF$GSEArank,decreasing=TRUE),]

   tx.entrez  <- bitr(gseaDF$SYMBOL,
                      fromType="SYMBOL",
                      toType="ENTREZID",
                      OrgDb="org.Mm.eg.db")

   tx.kegg    <- bitr_kegg(tx.entrez$ENTREZID,
                           fromType='ncbi-geneid',
                           toType='kegg',
                           organism='mmu')

   tx.kegg    <- rename(tx.kegg, c("kegg" = "KEGG",
                                   "ncbi-geneid" = "ENTREZID"))

   gseaDF      <- merge(gseaDF, tx.entrez, by="SYMBOL")
   gseaDF      <- merge(gseaDF, tx.kegg, by="ENTREZID")
   gseaDF      <- gseaDF[,c("SYMBOL","GSEArank","ENTREZID","KEGG")]

   gseaDF <- na.omit(gseaDF)
   tx <- as.character(gseaDF$ENTREZID)
   gseaDF <- gseaDF$GSEArank
   names(gseaDF) <- tx
   gseaDF <- sort(gseaDF,decreasing=TRUE)
   tx.entrez <- NULL; tx.kegg <- NULL
#   return(gseaDF)



   analysisOptions <- list(
      "enrichKEGG"      = list(gene = degDF$KEGG,
                               organism = "mmu",
                               pAdjustMethod = "BH"),

      "enrichGO"        = list(gene = degDF$SYMBOL,
                               OrgDb = org.Mm.eg.db,
                               keyType = 'SYMBOL'),

      "enrichPathway"   = list(gene = degDF$ENTREZID,
                               organism = "mouse",
                               pAdjustMethod = "BH",
                               readable=TRUE),

      "enricher"        = list(gene = degDF$ENTREZID,
                               pAdjustMethod = "BH"),

      "gseKEGG"         = list(geneList = gseaDF,
                               organism = "mmu",
                               nPerm = 10000,
                               pAdjustMethod = "BH",
                               verbose=FALSE),

      "gseGO"           = list(geneList = gseaDF,
                               OrgDB = org.Mm.eg.db,
                               minGSSize = 10,
                               nPerm = 10000,
                               pAdjustMethod = "BH",
                               verbose=FALSE),

      "gseGO"           = list(geneList = gseaDF,
                               OrgDB = org.Mm.eg.db,
                               minGSSize = 10,
                               nPerm = 10000,
                               pAdjustMethod = "BH",
                               verbose=FALSE),

      "gsePathway"      = list(geneList = gseaDF,
                               organism = "mouse",
                               nPerm = 10000,
                               pAdjustMethod = "BH",
                               verbose=FALSE)
   )

   setReadable.options <- list(
      "enrichKEGG" = list("org.Mm.eg.db","ENTREZID"),
      "enricher" = list('org.Mm.eg.db', 'ENTREZID'),
      "gseKEGG"  = list("org.Mm.eg.db", "ENTREZID"),
      "gseGO"    = list("org.Mm.eg.db", "ENTREZID"),
      "gsePathway"=list("org.Mm.eg.db", "ENTREZID"))

   analysisList <- list(
      "enrichKEGG" = list(
         "KEGG1_BHp20"  = list( pvalueCutoff = 0.2),

         "KEGG1_all"    = list( pvalueCutoff = 1,
                                qvalueCutoff = 1,
                                minGSSize = 1,
                                maxGSSize = 10000),

         "KEGG1_M_all"  = list( pvalueCutoff = 1,
                                qvalueCutoff = 1)
      ),

      "enrichGO" = list(
         "GO1_BP_BHp001_DF"     = list( ont = "BP",
                                        pvalueCutoff  = 0.001),

         "GO1_MF_BHp05"         = list( ont = "MF",
                                        pvalueCutoff  = 0.05)

      ),

      "enrichPathway" = list(
         "Reactome1_BHp05"      = list( pvalueCutoff = 0.05),
   
         "Reactome1_all"        = list(pvalueCutoff  = 1,
                                       qvalueCutoff = 1,
                                       minGSSize = 1,
                                       maxGSSize = 10000)
      ),

      "enricher" = list(
         "MsigDB_H_BHp20"       = list( pvalueCutoff = 0.2,
                                        TERM2GENE = msigdbr(
                                          species="Mus musculus",
                                          category="H")[,
                                             c("gs_name", "entrez_gene")]),
   
         "MsigDB_H_all"         = list( pvalueCutoff = 1,
                                        TERM2GENE = msigdbr(
                                          species="Mus musculus",
                                          category="H")[,
                                             c("gs_name", "entrez_gene")]),
   
         "MsigDB_C2_CGP_BHp0001"= list( pvalueCutoff = 0.0001,
                                        TERM2GENE = msigdbr(
                                          species="Mus musculus",
                                          category="C2",
                                          subcategory="CGP")[,
                                             c("gs_name", "entrez_gene")]),
   
         "MsigDB_C2_CP_BHp20"   = list( pvalueCutoff = 0.2,
                                        TERM2GENE = msigdbr(
                                          species="Mus musculus",
                                          category="C2",
                                          subcategory="CP")[,
                                             c("gs_name", "entrez_gene")]),
   
         "MsigDB_C2_CPb_BHp20"  = list( pvalueCutoff = 0.2,
                                        TERM2GENE = msigdbr(
                                          species="Mus musculus",
                                          category="C2",
                                          subcategory="CP:BIOCARTA")[,
                                             c("gs_name", "entrez_gene")]),

         "MsigDB_C3_TFT_BHp20"  = list( pvalueCutoff = 0.2,
                                        TERM2GENE = msigdbr(
                                          species="Mus musculus",
                                          category="C3",
                                          subcategory="TFT")[,
                                             c("gs_name", "entrez_gene")]),
   
         "MsigDB_C6_BHp20"      = list( pvalueCutoff = 0.2,
                                        TERM2GENE = msigdbr(
                                          species="Mus musculus",
                                          category="C6")[,
                                             c("gs_name", "entrez_gene")]),
   
         "MsigDB_C7_BHp0001"    = list( pvalueCutoff = 0.0001,
                                        TERM2GENE = msigdbr(
                                          species="Mus musculus",
                                          category="C7")[,
                                             c("gs_name", "entrez_gene")])
      ),

      "gseKEGG" = list(
         "gKEGG1_BHp10" = list(pvalueCutoff=0.1),
         "gKEGG1_all"   = list(pvalueCutoff=1,
                               minGSSize = 1,
                               maxGSSize = 10000)
      ),

#      "gseGO" = list(
#         "gGO1_BP_BHp05"= list(pvalueCutoff = 0.05,
#                               ont = "BP"),
#
#         "gGO1_BP_BHp05"= list(pvalueCutoff = 0.1,
#                               ont = "MF")
#      ),
#
      "gsePathway" = list(
         "gReactome1_BHp05" = list(pvalueCutoff = 0.05)
#         "gReactome1_all" = list(pvalueCutoff = 1,
#                                 minGSSize = 1,
#                                 maxGSSize = 10000)
      )
   )

   analysisList$GSEA <- analysisList$enricher
   names(analysisList$GSEA) <- paste0("g",names(analysisList$GSEA))
   for (type in names(analysisList$GSEA)) {
      analysisList$GSEA[[type]]$geneList <- gseaDF
   }

   res.list <- list() ;
   for (analysisType in names(analysisList)) {
      print(paste0("-> analyzing ",analysisType))

      if (grepl("^e",analysisType)) {next;}

      res.list[[analysisType]] <- list()
      for (analysis in names(analysisList[[analysisType]])) {
         print(paste0("  |-> ",analysis))

         curOptions <- analysisList[[analysisType]][[analysis]]
         for (option in names(analysisOptions[[analysisType]])) {
            curOptions[[option]] <- analysisOptions[[analysisType]][[option]]
         }

#         print(paste0("about to do.call on ", analysisType," with these options: ",curOptions))
         print(paste0("about to do.call on ", analysisType))
         print(paste0(names(curOptions),collapse=", "))
         
         res.list[[analysisType]][[analysis]] <- do.call(analysisType, curOptions)
         print(paste0("-> got back"))


         if (analysisType %in% names(setReadable.options)) {
            curOptions <- c(list(x=res.list[[analysisType]][[analysis]]),
                            setReadable.options[[analysisType]])

            res.list[[analysisType]][[analysis]] <- do.call(setReadable, curOptions)

            write.table(as.data.frame(res.list[[analysisType]][[analysis]]),
                        file=paste0(outDir,"/",analysisType,".",analysis,".txt"),
                        quote=FALSE,
                        sep="\t",
                        row.names=FALSE)

         }
      }
   }

# fpd 181218_1559 -- write out to summary file
#   write.table(hSummary$V1, file=paste(Output_file_prefix, "clusterProfiler_Summary.txt", sep="_"), sep='\t', quote=FALSE,row.names=FALSE)

   print("-> plotting: ")

   for (analysisType in names(res.list)) {

      for (analysis in names(res.list[[analysisType]])) {

         if (nrow(res.list[[analysisType]][[analysis]]) == 0) next;

         if (1) { #dotplot and emapplot for both DEG and GSEA
            outFile <- paste0(outDir,"/clusterProfiler_plots.",name,".",
                              analysisType,".",analysis,".dotplot.pdf")
            ggplot2::ggsave(filename=outFile,
                                plot = dotplot(
                                 res.list[[analysisType]][[analysis]],
                                 showCategory=nShow,
                                 font.size = 8,
                                 title = analysis),
                                width=8,
                                height=8,
                                units="in",
                                device=cairo_pdf)
   
            outFile <- paste0(outDir,"/clusterProfiler_plots.",name,".",
                              analysisType,".",analysis,".heatplot.pdf")
            ggplot2::ggsave(filename=outFile,
                                plot = heatplot(
                                 res.list[[analysisType]][[analysis]],
                                 showCategory=nShow)+ggplot2::ggtitle(analysis),
                                width=8,
                                height=8,
                                units="in",
                                device=cairo_pdf)

         }
   
         if (0 & grepl("^g|G", analysisType)) { # GSEA plots

            outFile <- paste0(outDir,"/clusterProfiler_plots.",name,".",
                              analysisType,".",analysis,".gseaplot2.pdf")
            ggplot2::ggsave(filename=outFile,
                                plot = gseaplot2(
                                 res.list[[analysisType]][[analysis]],
                                 title = analysis,
                                 geneSetID = 1:5,
                                 pvalue_table = TRUE,
                                 color = c("blue",
                                    "red",
                                    "green",
                                    "gold",
                                    "violet"),
                                 ES_geom = "line",
                                 base_size = 11,
                                 subplots = 1:3,
                                 rel_heights = c(2, 0.5, 1)
                                ),
                                width=8,
                                height=8,
                                units="in",
                                device=cairo_pdf)

         }
      }

   }

   return(res.list)

}


writeDEG <- function(dat) {

# Goal; write generic signature file

   for (degPair in names(dat$deg)) {
      outFn <- paste0(dat$specs$outDir,"/",degPair,"_signature.txt")
      headerLines <- c(paste0("# SIGNATURE: ",degPair),
                       "# AUTHOR: FRED P. DAVIS",
                       paste0("# CRITERIA: |FOLD CHANGE| >= ",dat$specs$thresh$deGenes.FC),
                       paste0("# CRITERIA: SLEUTH Q-VALUE <= ", dat$specs$thresh$sleuth.qval),
                       paste0("# CRITERIA: MINIMUM ",dat$specs$thresh$exprGenes.minTPM," TPM IN AT LEAST ONE SAMPLE")
                       )

      uniqCond <- unique(dat$deg[[degPair]]$designMat$condition)

      headerLines <- c(headerLines,
                       paste0("# LOG FOLD CHANGE: log2(" , gsub("^[0-9]\\.","",uniqCond[1]), " / ",
                                                           gsub("^[0-9]\\.","",uniqCond[2]), ")"))

      writeLines(headerLines, con=outFn, sep="\n")

      outMat <- dat$deg[[degPair]]$fc.deg[,c("gene_name","log2fc")]
      outMat <- outMat[order(outMat$gene_name),]
      outMat$log2fc <- round(outMat$log2fc, digits=2)
      write.table(outMat, outFn, append=TRUE,quote=FALSE,row.names=FALSE,col.names=TRUE)
   }

}
