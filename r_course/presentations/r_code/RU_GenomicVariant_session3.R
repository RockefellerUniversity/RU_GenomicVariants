params <-
list(isSlides = "no")

## ----setup_varManS3, include=FALSE--------------------------------------------
knitr::opts_chunk$set(echo = TRUE,cache = TRUE,cache.lazy = FALSE,message = FALSE,warning = FALSE, tidy = T)
# AsSlides <- TRUE
#
suppressPackageStartupMessages(library(VariantAnnotation))
suppressPackageStartupMessages(library(DT))
suppressPackageStartupMessages(library(BSgenome.Hsapiens.UCSC.hg19))
suppressPackageStartupMessages(library(TxDb.Hsapiens.UCSC.hg19.knownGene))
suppressPackageStartupMessages(library(SNPlocs.Hsapiens.dbSNP144.GRCh37))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(GenomicFeatures))
suppressPackageStartupMessages(library(maftools))
suppressPackageStartupMessages(library(NMF))


## ---- results='asis',include=TRUE,echo=FALSE----------------------------------
if(params$isSlides != "yes"){
  cat("# Genomic Variants (part 3)

---
"    
  )
  
}



## ---- results='asis',include=TRUE,echo=FALSE----------------------------------
if(params$isSlides == "yes"){
  cat("class: inverse, center, middle

# Working with MAF files

<html><div style='float:left'></div><hr color='#EB811B' size=1px width=720px></html> 

---
"    
  )
}else{
  cat("# Working with MAF files

---
"    
  )
  
}



## ----mult_mafIntro_advan------------------------------------------------------
laml_tab <- read.delim("data/tcga_laml.maf",sep="\t")
laml_tab[1:2,]


## ----mult_samInfo_advan,fig.align="center"------------------------------------
tbl <- table(laml_tab$Tumor_Sample_Barcode)
hist(tbl,breaks = 10, xlab = "Mutations")


## ----mafTools_intro,echo=FALSE,out.width = "50%",fig.align="center"-----------
knitr::include_graphics("imgs/vcfMan_fig7r.png")


## ----mult_mafT_advan----------------------------------------------------------
library(maftools)
laml <- read.maf("data/tcga_laml.maf.gz")
class(laml)


## ----mult_samSum_advan--------------------------------------------------------
sample_sum <- getSampleSummary(laml)
sample_sum[1:2,]


## ----mult_samSumPlot_advan1---------------------------------------------------
var_to <- sample_sum$total
names(var_to) <- sample_sum$Tumor_Sample_Barcode
sample_sum <- dplyr::select(sample_sum,-total)
melt_dat <- reshape2::melt(sample_sum, id="Tumor_Sample_Barcode")
melt_dat[1:3,]


## ----mult_samSumPlot_advan15--------------------------------------------------
melt_dat$totalVar <- var_to[match(melt_dat$Tumor_Sample_Barcode,names(var_to))]
melt_dat$prop <- melt_dat$value / melt_dat$totalVar
head(melt_dat)


## ----mult_samSumPlot_advan2,echo=TRUE,eval=FALSE,tidy=FALSE-------------------
## ggplot(melt_dat,aes(x=Tumor_Sample_Barcode,y=value,fill=variable))+
##   geom_bar(stat='identity',position = 'stack')+
##   labs(x="",y="Mutations",fill="")+
##   theme(axis.text.x=element_blank())


## ----mult_samSumPlot_advan3,echo=FALSE,eval=TRUE,tidy=TRUE,fig.align="center"----
ggplot(melt_dat,aes(x=Tumor_Sample_Barcode,y=value,fill=variable))+
  geom_bar(stat='identity',position = 'stack')+
  labs(x="",y="Mutations",fill="")+
  theme(axis.text.x=element_blank())


## ----mult_samSumPlotPres_advan4,echo=TRUE,eval=FALSE,tidy=FALSE---------------
## ggplot(melt_dat,aes(x=Tumor_Sample_Barcode,y=prop,fill=variable))+
##   geom_bar(stat='identity',position = 'stack')+
##   labs(x="",y="Proportion",fill="")+
##   theme(axis.text.x=element_blank())


## ----mult_samSumPlotPres_advan5,echo=FALSE,eval=TRUE,tidy=FALSE,fig.align="center"----
ggplot(melt_dat,aes(x=Tumor_Sample_Barcode,y=prop,fill=variable))+
  geom_bar(stat='identity',position = 'stack')+
  labs(x="",y="Proportion",fill="")+
  theme(axis.text.x=element_blank())


## ----mult_geneSum_advan-------------------------------------------------------
gene_sum <- getGeneSummary(laml)
gene_sum[1:2,]


## ---- results='asis',include=TRUE,echo=FALSE----------------------------------
if(params$isSlides == "yes"){
  cat("class: inverse, center, middle

# Functional analysis of mutations

<html><div style='float:left'></div><hr color='#EB811B' size=1px width=720px></html> 

---
"    
  )
}else{
  cat("# Functional analysis of mutations

---
"    
  )
  
}



## ----mult_lolli_advan_eval1,echo=TRUE,eval=FALSE,tidy=FALSE-------------------
## lollipopPlot(maf = laml,
##              gene = 'NRAS',
##              AACol = 'Protein_Change',
##              showMutationRate = TRUE,
##              labelPos = "all")


## ----mult_lolli_advan_eval2,echo=FALSE,eval=TRUE,tidy=FALSE,fig.align="center"----
lollipopPlot(maf = laml, gene = 'NRAS', AACol = 'Protein_Change', showMutationRate = TRUE,labelPos = "all")


## ----mult_oncoplot_advan,fig.align="center"-----------------------------------
oncoplot(maf=laml,top = 5)


## ----mult_pathPlotWF_advan1,eval=TRUE,echo=TRUE,tidy=FALSE,fig.show="hide"----
OncogenicPathways(maf = laml)


## ----mult_pathPlotWF_advan2,eval=TRUE,echo=FALSE,tidy=FALSE,results="hide",fig.align="center"----
OncogenicPathways(maf = laml)


## ----mult_pathPlotWF_advan3,echo=TRUE,eval=TRUE,tidy=FALSE,fig.align="center"----
PlotOncogenicPathways(maf = laml, pathways = "RTK-RAS")


## ---- results='asis',include=TRUE,echo=FALSE----------------------------------
if(params$isSlides == "yes"){
  cat("class: inverse, center, middle

# Mutation Signatures

<html><div style='float:left'></div><hr color='#EB811B' size=1px width=720px></html> 

---
"    
  )
}else{
  cat("# Mutation Signatures

---
"    
  )
  
}



## ----mult_mutSig_TiTv_advan1,eval=FALSE,echo=TRUE,tidy=FALSE------------------
## laml.titv = titv(maf = laml, plot = FALSE, useSyn = TRUE)
## plotTiTv(res = laml.titv)


## ----mult_mutSig_TiTv_advan2,eval=TRUE,echo=FALSE,tidy=FALSE,fig.align="center"----
laml.titv = titv(maf = laml, plot = FALSE, useSyn = TRUE)
plotTiTv(res = laml.titv)


## ----mult_mutSig_triMut_advan,tidy=FALSE--------------------------------------
library(BSgenome.Hsapiens.UCSC.hg19, quietly = TRUE)
laml.tnm = trinucleotideMatrix(maf = laml,
                               prefix = 'chr', 
                               add = TRUE, 
                               ref_genome = "BSgenome.Hsapiens.UCSC.hg19")


## ----mult_mutSig_triMutPres_advan---------------------------------------------
dim(laml.tnm$nmf_matrix)
laml.tnm$nmf_matrix[1,]


## ----mult_mutSig_triMutPat_advan1---------------------------------------------
tarSam_triNuc <- laml.tnm$nmf_matrix['TCGA-AB-3009',]
tarSam_triNuc[1:2]


## ----mult_mutSig_triMutPat_advan2---------------------------------------------
yd <- data.frame(triNuc=names(tarSam_triNuc),
                 count=tarSam_triNuc,
                 stringsAsFactors = FALSE)
yd$cat <- gsub("(.*)\\[(.*)\\](.*)","\\2",yd$triNuc)
yd$num <- seq(1,length(yd$triNuc))


## ----mult_mutSig_triMutPat_advan3,eval=FALSE,echo=TRUE,tidy=FALSE-------------
## ggplot(yd,aes(x=num,y=count,fill=cat))+
##   geom_bar(stat='identity')+
##   labs(x="",y="Counts",fill="")+
##   theme(axis.text.x=element_blank())


## ----mult_mutSig_triMutPat_advan4,eval=TRUE,echo=FALSE,tidy=FALSE,fig.align="center"----
ggplot(yd,aes(x=num,y=count,fill=cat))+geom_bar(stat='identity')+
  labs(x="",y="Counts",fill="")+theme(axis.text.x=element_blank())


## ----mult_mutSig_sigEst_advan1,eval=FALSE,echo=TRUE,tidy=FALSE----------------
## library('NMF')
## laml.sign <- estimateSignatures(mat = laml.tnm,
##                                 nTry = 6,
##                                 pConstant = 0.1,
##                                 parallel = 1)


## ----mult_mutSig_sigEst_advan2,eval=TRUE,echo=FALSE,include=FALSE-------------
library('NMF')
laml.sign <- estimateSignatures(mat = laml.tnm,
                                nTry = 6,
                                pConstant = 0.1,
                                parallel = 1)


## ----mult_mutSig_sigEst_advan3,eval=TRUE,echo=FALSE,include=TRUE,fig.align="center"----
plotCophenetic(laml.sign)


## ----mult_mutSig_sigExt_advan,tidy=FALSE--------------------------------------
laml.sig.ext <- extractSignatures(mat = laml.tnm, 
                                  n = 3,
                                  pConstant = 0.1,
                                  parallel = 1)
laml.sig.ext$signatures[1:5,] # use for mapping to mutational signature database


## ----mult_muSig_mapSig_advan,tidy=FALSE---------------------------------------
laml.og30.cosm = compareSignatures(nmfRes = laml.sig.ext,
                                   sig_db = "legacy")
laml.og30.cosm$cosine_similarities[,1:5]


## ----mult_muSig_mapSigPres_advan,fig.align="center"---------------------------
pheatmap::pheatmap(mat = laml.og30.cosm$cosine_similarities, cluster_rows = FALSE)


## ----mult_muSig_mapSigSBS_advan1,fig.align="center"---------------------------
laml.sign.sbs = compareSignatures(nmfRes = laml.sig.ext, sig_db = "SBS")
laml.sign.sbs$cosine_similarities[,1:5]


## ----mult_muSig_mapSigSBS_advan2,fig.align="center"---------------------------
pheatmap::pheatmap(mat = laml.sign.sbs$cosine_similarities, cluster_rows = FALSE)


## ----mult_muSig_plotSigCOS_advan1,eval=FALSE,echo=TRUE,tidy=FALSE-------------
## plotSignatures(nmfRes = laml.sig.ext,
##                title_size = 1.2,
##                contributions = FALSE,
##                show_title = TRUE,
##                sig_db = 'legacy')


## ----mult_muSig_plotSigCOS_advan2,eval=TRUE,echo=FALSE,tidy=FALSE,fig.align="center"----
plotSignatures(nmfRes = laml.sig.ext, 
               title_size = 1.2,
               contributions = FALSE,
               show_title = TRUE,
               sig_db = 'legacy')


## ----mult_muSig_plotSigSBS_advan1,eval=FALSE,echo=TRUE,tidy=FALSE-------------
## plotSignatures(nmfRes = laml.sig.ext,
##                title_size = 1.2,
##                contributions = FALSE,
##                show_title = TRUE,
##                sig_db = 'SBS')


## ----mult_muSig_plotSigSBS_advan2,eval=TRUE,echo=FALSE,tidy=FALSE,fig.align="center"----
plotSignatures(nmfRes = laml.sig.ext, 
               title_size = 1.2,
               contributions = FALSE,
               show_title = TRUE,
               sig_db = 'SBS')


## ----mult_muSig_plotSigSAM_advan1,eval=FALSE,echo=TRUE,tidy=FALSE-------------
## plotSignatures(nmfRes = laml.sig.ext,
##                title_size = 0.8,
##                contributions = TRUE,
##                show_title = TRUE)


## ----mult_muSig_plotSigSAM_advan2,eval=TRUE,echo=FALSE,tidy=FALSE,fig.align="center"----
plotSignatures(nmfRes = laml.sig.ext,
               title_size = 0.8,
               contributions = TRUE,
               show_title = TRUE)


## ----mult_muSig_enrich_advan1,eval=FALSE,echo=TRUE,tidy=FALSE-----------------
## laml.se = signatureEnrichment(maf = laml,
##                               sig_res = laml.sig.ext)


## ----mult_muSig_enrich_advan2,eval=TRUE,echo=FALSE,tidy=FALSE,fig.align="center",results=FALSE----
laml.se = signatureEnrichment(maf = laml, sig_res = laml.sig.ext)


## ----mult_muSig_enrichGene_advan----------------------------------------------
laml.se$groupwise_comparision[1:2,]


## ----mult_muSig_enrichGenePlot_advan1,fig.align="center",out.height="70%"-----
plotEnrichmentResults(enrich_res = laml.se, pVal = 0.05)

