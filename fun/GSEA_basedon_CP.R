#GSEA based on clusterProfiler
# reorganize data for dotplot visualization
# author: Chang Lu
# c.lu@maastrichtuniversity.nl

library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(stringr)
library("ReactomePA")
library(dplyr)
library(foreach)
library(doParallel)

wpgmtfile <- "./referencefile/wikipathways-20210110-gmt-Homo_sapiens.gmt"
wp2gene <- read.gmt(wpgmtfile)
wp2gene <- wp2gene %>% tidyr::separate(term, c("name","version","wpid","org"), "%")
wpid2gene <- wp2gene %>% dplyr::select(wpid, gene) #TERM2GENE
wpid2name <- wp2gene %>% dplyr::select(wpid, name) #TERM2NAME


# GSEA of GO, KEGG, REACTOME, WP
# GES_df: female or male GES
# ifsim: if you want to simplify the GO results
getdfforplot <- function(GES_df, pvalueCutoff = 0.05, ifsim = FALSE){
  # GES_df <- switch(gender,
  #                  male = GES_ma_df,
  #                  female = GES_fe_df)
  genelist <- GES_df$t
  names(genelist) <- rownames(GES_df)
  genelist <- sort(genelist,decreasing = TRUE)
  
  # GO GSEA
  gse = list(GOBP = NULL,
             GOCC = NULL,
             GOMF = NULL,
             KEGG = NULL,
             WP = NULL,
             REACTOME = NULL)
  gse_sim = list(GOBP = NULL,
                 GOCC = NULL,
                 GOMF = NULL,
                 KEGG = NULL,
                 WP = NULL,
                 REACTOME = NULL)
  
  for (i in names(gse)[1:3]) {
    gse[[i]] <- gseGO(geneList=genelist, 
                      ont = str_remove(i,'GO'), 
                      keyType ="ENTREZID", 
                      pvalueCutoff = pvalueCutoff, 
                      pAdjustMethod = "fdr",
                      verbose = TRUE, 
                      OrgDb = org.Hs.eg.db)
    if(ifsim == TRUE){
      gse_sim[[i]] <- clusterProfiler::simplify(gse[[i]], cutoff=0.5, by="p.adjust", select_fun=min)
    }else{
      gse_sim[[i]] <- gse[[i]]
    }
  }
  
    
  gse_sim[['KEGG']] <- gseKEGG(geneList=genelist, 
                           pvalueCutoff = pvalueCutoff, 
                           pAdjustMethod = "fdr",
                           verbose = TRUE)
  gse_sim[['WP']] <- GSEA(genelist, TERM2GENE = wpid2gene, TERM2NAME = wpid2name, verbose=FALSE)
  
  gse_sim[['REACTOME']] <- gsePathway(genelist, 
                                organism = "human",  pvalueCutoff = pvalueCutoff, pAdjustMethod = "fdr")
  
  dfforplot = list(GOBP = NULL,
                   GOCC = NULL,
                   GOMF = NULL,
                   KEGG = NULL,
                   WP = NULL,
                   REACTOME = NULL)
  for (i in names(gse_sim)) {
    if (nrow(gse_sim[[i]]@result!= 0)){
      dfforplot[[i]] <- enrichplot:::fortify.gseaResult(
        model= gse_sim[[i]],
        showCategory = 10,
        by = "geneRatio",
        split = ".sign",
        includeAll = TRUE)
    }
  }
  return(dfforplot)
}



# combine and reorganize the GSEA results of female and male cohort for dotplot visualization
# femalePathwayList, malePathwayList: female's and male's output of 'getdfforplot'
reorganClusterProfile <- function(femalePathwayList, malePathwayList){
  
  fedf <- do.call('rbind',lapply(1:length(femalePathwayList),function(x){
    if(!is.null(femalePathwayList[[x]]) ){
      Y = femalePathwayList[[x]]
      rownames(Y) = stringr::str_c(rownames(Y),'_fe')
       return(Y)}
  } ))
  
  classfe <-  do.call('rbind', lapply(1:length(femalePathwayList),function(x)
    if(!is.null(femalePathwayList[[x]]) ){ cbind(rep('female' , times = nrow(femalePathwayList[[x]])),
                                           rep(names(femalePathwayList)[x], times = nrow(femalePathwayList[[x]]))) }
    ))

  
  
  madf <- do.call('rbind',lapply(1:length(malePathwayList),function(x){
    if(!is.null(malePathwayList[[x]]) ){
    Y = malePathwayList[[x]]
    rownames(Y) = stringr::str_c(rownames(Y),'_fe')
    return(Y)}
  } ))
  
  classma <-  do.call('rbind', lapply(1:length(malePathwayList),function(x)
    if(!is.null(malePathwayList[[x]]) ){ cbind(rep('male' , times = nrow(malePathwayList[[x]])),
                                                 rep(names(malePathwayList)[x], times = nrow(malePathwayList[[x]]))) }
  ))
  
  
  df = rbind(fedf,madf) 
  cla = data.frame(rbind(classfe,classma),row.names = rownames(df))
  colnames(cla) <- c('gender','source')
  PathwayAll_df <- merge(df,cla,all=T,by='row.names')
  PathwayAll_df <- PathwayAll_df %>%
    dplyr::rename(pathway = Description, AdjPvalu = p.adjust) %>%
    column_to_rownames(var = "Row.names") %>%
    add_column( class = stringr::str_c(PathwayAll_df$gender, "(", PathwayAll_df$.sign, ")"))

  return(PathwayAll_df)
}



