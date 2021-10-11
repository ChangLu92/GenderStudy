library(viper)
library(ggplot2)

## organize regulon data for viper enrichment analysis
df2regulon <- function(df) {
  regulon_list = split(df, df$tf)
  viper_regulons = lapply(regulon_list, function(regulon) {
    tfmode = stats::setNames(regulon$mor, regulon$target)
    list(tfmode = tfmode, likelihood = rep(1, length(tfmode)))
  })
  return(viper_regulons)
}

# # gene permutation
# ViperPerm <- function(GES_mat, regulons, k = 1000, z_scores =  TRUE, minsize = 5, eset.filter = FALSE,  cores = 1, verbose = FALSE, nes = TRUE){
#   scores <- viper(GES_mat, regulons, minsize = 5, eset.filter = FALSE,  cores = 1, verbose = FALSE, nes = TRUE)
#   scores <- scores %>% as.data.frame()
#   
#   null_dist_t <- replicate(k, sample(GES_mat,length(GES_mat),  replace = FALSE))
#   rownames(null_dist_t) = rownames(GES_mat)
#   
#   null_dist_scores <- apply(null_dist_t, 2, function(x) viper(x, regulons, minsize = 5, eset.filter = FALSE,  cores = 1, verbose = FALSE, nes = TRUE))
#   
#   res = data.frame(row.names = rownames(scores))
#   
#   if(z_scores) {
#     scores$mean <- apply(null_dist_scores,1,mean)
#     scores$sd <- apply(null_dist_scores,1,sd)
#     resListCurrent <- (scores[,1]-scores[,2])/scores[,3]
#     names(resListCurrent) <- rownames(scores)
#     pvalue2sided=2*pnorm(-abs(resListCurrent)) 
#     res$NES = resListCurrent
#     res$PV = pvalue2sided
#     
#   } else {
#     for(j in seq(1, nrow(null_dist_scores))) {
#       ecdf_function <- ecdf(null_dist_scores[j,])
#       scores[j,1] <- ecdf_function(scores[j,1])
#     }
#     score_probas <- scores*2-1
#     resListCurrent <- score_probas[,1]
#     names(resListCurrent) <- rownames(scores)
#     pvalue2sided= ifelse(resListCurrent>0, (1-resListCurrent)/2, (1+resListCurrent)/2)
#     res$NES = resListCurrent
#     res$PV = pvalue2sided    
#   }
#   return(res)
# }


# gene permutation
viperNullgene <- function(expset, per= 10000){
  dnull <- replicate(per, sample(expset,length(expset),  replace = FALSE))
  rownames(dnull) = rownames(expset)
  return(dnull)
}



# sample permutation
viperNulllimma <- function(expset, clin, per= 1000, repos= TRUE){
  dnull <- sapply(1:per, function(i, expset, clin, repos){
    lengthcvd =  sum(clin$indication== 'CVD')
    lengthhealth =  sum(clin$indication== 'healthy control')
    repeat{
      sorder <- sample(ncol(expset), replace=repos)
      if (length(unique(sorder[1:lengthcvd]))>1 & length(unique(sorder[-(1:lengthcvd)]))>1) break
    }
    newexpset <- filterColMatrix(expset, sorder)
    batch <- as.factor(clin$batch)
    design <- model.matrix(~0+clin$indication + batch + clin$age)
    colnames(design)[1:2] = c('CVD','Health')
    colnames(design)[ncol(design)] = 'age'
    fit.limma <- lmFit(newexpset, design)
    cont.matrix <- makeContrasts(
      CVDvsHealthy = CVD-Health,
      levels = design)
    fit.limma2 <- contrasts.fit(fit.limma, cont.matrix)
    fit.limma2 <- eBayes(fit.limma2)
    DEG_limma = topTable(fit.limma2, coef="CVDvsHealthy",   number = 20000, p.value = 1, adjust.method="fdr")
    t <- DEG_limma$t
    return(t)
  }, expset=expset, clin=clin, repos=repos )
  rownames(dnull) = rownames(expset)
  return(dnull)
}

#------------------------------------
# volcano plot for visualizing  a TF regulon's t-values or logFCs of females  or  males 
# highlighted the significant genes' name
# inputs:
# GES_HGNC: GES data from limma (e.g. GES_ma_HGNC)
# tfname: TF name 
# filename : the path and filename where you want to save your pic
# xaxis: the column you want to use for visualization (t or logFC)
TFvolcanoplot <- function(GES_HGNC, tfname, filename, xaxis = 't'){
  
  targets <- regulons$target[regulons$tf == tfname]
  GES_HGNC$diffexpressed <- NA
  GES_HGNC$diffexpressed[GES_HGNC[,xaxis] > 0 & GES_HGNC$P.Value < 0.05] <- "UP"
  GES_HGNC$diffexpressed[GES_HGNC[,xaxis] < 0 & GES_HGNC$P.Value < 0.05] <- "DOWN"
  # GES_HGNC$delabel <- NA
  # GES_HGNC$delabel[!is.na(GES_HGNC$diffexpressed)] <- GES_HGNC$ID[!is.na(GES_HGNC$diffexpressed)]
  GES_HGNC$delabel <- GES_HGNC$ID
  
  # jpeg(paste0("./results/limma+batch+age/CVD_vs_Control/",filename," ",tfname," volcanoplot.jpg"), 
  #      width = 800, height = 1000, res = 160)
  if(is.na(GES_HGNC$diffexpressed)){
    p = ggplot(data=GES_HGNC[GES_HGNC$ID %in% targets,], aes_string(x=xaxis, y="-log(P.Value)", label="delabel"))
  }else{
    p = ggplot(data=GES_HGNC[GES_HGNC$ID %in% targets,], aes_string(x=xaxis, y="-log(P.Value)", col="diffexpressed", label="delabel"))
  }
  # print( 
  p <-  p +
    geom_point(aes(col=diffexpressed)) +
    scale_colour_manual(name = "Differential \n Expressed Genes", values = c("UP" = "red", "DOWN" = "blue"), na.translate = TRUE, na.value = "grey")+ theme_bw()+
    theme(
      axis.title = element_text(size = 14),
      axis.text = element_text(size = 12),
      legend.title = element_text(size = 14),
      legend.text = element_text(size = 14))+
    geom_text_repel(aes(color = diffexpressed),fontface = "bold" , fontsize = 4)+
    ggtitle(paste("Volcano Plot of",tfname, 'on', filename)) 
  # dev.off()
  return(p)
}


#------------------------------------
#scatter for comparing a TF regulon's t-value or logFC of females and males 
# highlighted the significant genes' name
# inputs:
# GES_HGNC_ALL: merged data of female and male's GES 
# tfname: TF name 
# label.by : the column you want to use for defining 'sex-specific' genes (e.g. 'adj.P.Val')
# cutoff : cutoff for defining 'sex-specific' genes (e.g. 'adj.P.Val')
# axis: the column you want to use for visualization (t or logFC)
TF2dplot <- function(GES_HGNC_ALL, tfname,  axis = 't', label.by='adj.P.Val', cutoff = 0.05){
  cyl = data.frame( rep('Others',nrow(GES_HGNC_ALL)) ,row.names = GES_HGNC_ALL$ID) 
  cyl[GES_HGNC_ALL[,paste0(label.by,'.fe')]<cutoff,] <-  'female significant' 
  cyl[GES_HGNC_ALL[,paste0(label.by,'.ma')]<cutoff,] <- 'male significant'
  cyl[GES_HGNC_ALL[,paste0(label.by,'.ma')]<cutoff & GES_HGNC_ALL[,paste0(label.by,'.fe')]<cutoff,] <-  'Male & Female significant' 
  GES_HGNC_ALL$labels <- factor(cyl[,1])
  
  targets <- regulons$target[regulons$tf == tfname]
  GES_target <- GES_HGNC_ALL[GES_HGNC_ALL$ID %in% targets,]
  p <- ggplot(data = GES_target, aes_string(x=paste0(axis,'.ma'), y=paste0(axis,'.fe'))) + 
    ggtitle(tfname) + 
    geom_point(aes(col=labels)) + 
    scale_colour_manual(name = "Differential Expressed Genes",
                        values = c('female significant' = "red",
                                   'male significant' = "blue",
                                   'Male & Female significant' = 'black',
                                   'Others' = 'grey')) +
    theme_bw()+
    xlab(paste('Male',axis))+ylab(paste('Female',axis))+
    geom_text_repel(aes(label = ID, color = labels),fontface = "bold", data = dplyr::filter(GES_target, !is.na(labels)))+
    theme(
      axis.title = element_text(size = 14),
      axis.text = element_text(size = 12),
      legend.title = element_text(size = 14),
      legend.text = element_text(size = 14),
      legend.position="right")
  
  return(list(p,GES_target))
}

#------------------------------------
# Over respresentation analysis for gender-specific TFs
# geneset
# regulons
# NESname
# PVname
# genename
# cases
# bg
# pvcutoff
TF_ORA <- function(geneset,
                  regulons,
                  NESname = NULL,
                  PVname = NULL,
                  genename = NULL,
                  cases, 
                  bg,
                  pvcutoff = 0.05){
  if(is.null(PVname)){
    selectedgeneset <- switch(  
      cases[1],  
      up= geneset[geneset[,NESname]>0, genename],
      dn= geneset[geneset[,NESname]<0, genename] )
  }else{
    selectedgeneset <- switch(  
      cases[1],  
      up= geneset[(geneset[,NESname]>0 & geneset[,PVname]< pvcutoff), genename],
      dn= geneset[(geneset[,NESname]<0 & geneset[,PVname]< pvcutoff), genename])
  }
  pathwayfilelist <- switch(  
    cases[2],  
    pathway= pathways,
    GO= GOFILE
  ) 
  
  candidate_genes <- unique(c(selectedgeneset,regulons$target[regulons$tf %in% selectedgeneset]))
  sig_pathways <- runGSAhyper(genes = candidate_genes, universe = bg, gsc = loadGSC(pathwayfilelist),pcutoff = 1)
  sig_pathways_df <- as.data.frame(sig_pathways$resTab)  %>% tibble::rownames_to_column(var = "pathway")
  sig_pathways_df <- data.frame(t(apply(sig_pathways_df, 1, function(r){
    aux = unlist(strsplit( sub("_",";", r["pathway"]), ";" ))
    r["pathway"] = gsub("_", " ", aux[2])
    return(c(r, "source" = aux[1]))})))
  colnames(sig_pathways_df)[length(sig_pathways_df)] = 'source'
  
  #data for plotting
  PathwaysSelect <- sig_pathways_df %>%
    dplyr::rename(pvalue = `p.value`, AdjPvalu = `Adjusted.p.value`) %>% 
    dplyr::mutate(pathway = as.factor(pathway))
  return(PathwaysSelect)
  
}


