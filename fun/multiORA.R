#ORA AND GSEA
# Author : Chang Lu
# c.lu@maastrichtuniversity.nl

library(piano)
library(ggnewscale)
library(dplyr)
library(stringr)

GOFILE = './referencefile/c5.go.v7.4.symbols.gmt'
pathways = "./referencefile/c2.cp.v7.2.symbols.gmt"


# ORA based on Piano
# geneset: GES table from limma
# NESname: the column of  geneset that you want to use for selecting the Differential expressed genes (e.g. t-value or logFC)
# genename: symbol or Gene ID
# bg: background gene list
multiORA = function(geneset,
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
      dn= geneset[geneset[,NESname]<0, genename],
      all = geneset[, genename])
  }else{
    selectedgeneset <- switch(  
      cases[1],  
      up= geneset[(geneset[,NESname]>0 & geneset[,PVname]< pvcutoff), genename],
      dn= geneset[(geneset[,NESname]<0 & geneset[,PVname]< pvcutoff), genename],
      all = geneset[, genename])
  }
  
  pathwayfilelist <- switch(  
    cases[2],  
    pathway= pathways,
    GO= GOFILE
  ) 
  
  candidate_genes <- unique(selectedgeneset)
  sig_pathways <- runGSAhyper(genes = candidate_genes, universe = unique(bg) , gsc = loadGSC(pathwayfilelist),pcutoff = 1)
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



# # GSEA based on Piano
# multiGSEA = function(geneset,
#                      NESname = NULL,
#                      genename,
#                      cases){
#   
#   pathwayfilelist <- switch(  
#     cases,  
#     pathway= pathways,
#     GO= GOFILE
#   ) 
#   selectedgeneset <- geneset[,NESname]
#   names(selectedgeneset) <- geneset[,genename]
#   
#   sig_pathways <- runGSA(geneLevelStats = selectedgeneset,  
#                          gsc = loadGSC(pathwayfilelist),
#                          adjMethod = 'fdr')
#   
#   sig_pathways_df <- as.data.frame(sig_pathways$resTab)  %>% tibble::rownames_to_column(var = "pathway")
#   sig_pathways_df <- data.frame(t(apply(sig_pathways_df, 1, function(r){
#     aux = unlist(strsplit( sub("_",";", r["pathway"]), ";" ))
#     r["pathway"] = gsub("_", " ", aux[2])
#     return(c(r, "source" = aux[1]))})))
#   colnames(sig_pathways_df)[length(sig_pathways_df)] = 'source'
#   
#   #data for plotting
#   PathwaysSelect <- sig_pathways_df %>%
#     dplyr::rename(pvalue = `p.value`, AdjPvalu = `Adjusted.p.value`) %>% 
#     dplyr::mutate(pathway = as.factor(pathway))
#   return(PathwaysSelect)
#   
# }




##---------------
## combine several enrichment result together for visualization
## PathwayLists: a list of enrichment results (from 'multiORA' or 'multiGSEA')
reorganisedf <- function(PathwayLists){
  
  for (i in 1:length(PathwayLists)) {
    pl = PathwayLists[[i]]
    modulename = names(PathwayLists)[[i]]
    
    df <- do.call('rbind',lapply(1:length(pl),function(x) pl[[x]]))
    modulenames <- rep(modulename, times = nrow(df))
    class <- do.call('c', lapply(1:length(pl),
                                 function(x){
                                   rep(paste0(modulename,'_',names(pl)[x]), 
                                       times = nrow(pl[[x]]))}))
    if(i == 1){
      dfall = df
      classall = class
      modules = modulenames
    }else{
      dfall = rbind(dfall,df) 
      classall = c(classall,class)
      modules = c(modules,modulenames)
    }
    
  }
  
  sign = vector(mode="character",length = length(classall))
  sign[str_detect(classall,'up' )] = 'activated'
  sign[str_detect(classall,'dn' )] = 'suppressed'
  sign[str_detect(classall,'all' )] = 'all'
  # unique(sign)
  
  PathwayAll_df<- dfall %>%
    tibble::add_column(class = factor(classall), modules = modules ,.sign = sign)%>%  
    dplyr::rename(Count = `Significant..in.gene.set.`, OtherintheGeneset = `Non.significant..in.gene.set.`) 
  PathwayAll_df[,2:7] <- sapply(dfall[,2:7], as.numeric)
  
  return(PathwayAll_df)
}




##---------------------------------------------------------------
# ENRICHMETN ANALYSIS VISUALIZATION ; DOT PlOTS
# PathwayAll_df: output of function'reorganisedf'
showCompEnrichDotplot <- function(PathwayAll_df,
                                  pathwayname = 'GO',
                                  split = FALSE,
                                  top = 10,
                                  pcutoff = 0.05,
                                  filename = NULL,
                                  wid = 1200,
                                  break_by = 2,
                                  height = 1000){
  
  if(pathwayname == 'GO'){
    PathwayAll_df_filter = PathwayAll_df %>% 
      dplyr::arrange(AdjPvalu) %>%
      dplyr::filter(AdjPvalu < pcutoff & source %in% c("GOBP" , "GOMF","GOCC") )
  }else if(pathwayname == 'pathway'){
    PathwayAll_df_filter = PathwayAll_df %>% 
      dplyr::arrange(AdjPvalu) %>%
      dplyr::filter(AdjPvalu < pcutoff & !(source %in% c("GOBP" , "GOMF","GOCC")))}
  
  # if(ifClusterProfile){
  #   label = c('female(suppressed)'='female\ndn',
  #                         'male(suppressed)'='male\ndn',
  #                         "female(activated)"= "female\nup",
  #                         "male(activated)" = "male\nup")
  # }
  
  # }else if(pathwayname == 'ClusterProfile'){
  #   PathwayAll_df_filter = PathwayAll_df %>% 
  #     dplyr::arrange(AdjPvalu) %>%
  #     dplyr::filter(AdjPvalu < pcutoff)
  #    label = c('female(suppressed)'='female\ndn',
  #             'male(suppressed)'='male\ndn',   
  #             "female(activated)"= "female\nup", 
  #             "male(activated)" = "male\nup")
  
  # }
  
  
  clProf.reshape.df <- do.call("rbind", lapply(unique(as.character(PathwayAll_df_filter$class)), function(x, top){
    a <- PathwayAll_df_filter %>%
      dplyr::filter(class %in% x) %>% 
      dplyr::arrange(AdjPvalu) %>%
      mutate(logpv = -log10(AdjPvalu),)
    # a$pathway <- factor(a$pathway, levels =a$pathway[order(a$class)])
    
    
    a <- a[1:min(c(top,nrow(a))),]
  },top = top))
  
  
  
  
  # updata <- clProf.reshape.df %>% dplyr::filter(.sign == 'activated')
  # updata$pathway <- factor(updata$pathway, levels = updata$pathway[order(updata$class)])
  
  # downdata <- clProf.reshape.df %>% dplyr::filter(.sign == 'suppressed') 
  # downdata$pathway <- factor(downdata$pathway, levels = downdata$pathway[order(downdata$AdjPvalu)])
  
  
  clProf.reshape.df$logadjpv[clProf.reshape.df$.sign=='suppressed'] = log(clProf.reshape.df$AdjPvalu[clProf.reshape.df$.sign=='suppressed'])
  clProf.reshape.df$logadjpv[clProf.reshape.df$.sign=='activated'] = -log(clProf.reshape.df$AdjPvalu[clProf.reshape.df$.sign=='activated'])
  
  clProf.reshape.df = clProf.reshape.df[clProf.reshape.df$.sign!='all',]
  
  clProf.reshape.df$pathway = sapply( strwrap(clProf.reshape.df$pathway, 50, simplify=FALSE), paste, collapse="\n" )
  
  if(floor(min(clProf.reshape.df$logadjpv))>0){
    brs = seq(0 , floor(max(clProf.reshape.df$logadjpv)) ,by=break_by)
    labels = seq(0 , floor(max(clProf.reshape.df$logadjpv)) ,by=break_by)
  }else if(floor(max(clProf.reshape.df$logadjpv))<0){
    brs = seq(floor(min(clProf.reshape.df$logadjpv)) ,0, by=break_by)
    labels = seq(floor(min(clProf.reshape.df$logadjpv)) ,0, by=break_by)
  }else{
    brs = seq( floor(min(clProf.reshape.df$logadjpv))  , floor(max(clProf.reshape.df$logadjpv)) ,by=break_by)
    labels = abs(seq( floor(min(clProf.reshape.df$logadjpv))  , floor(max(clProf.reshape.df$logadjpv)) ,by=break_by))
  }
  
  
  p<- ggplot()+
    geom_point(data = clProf.reshape.df , aes(x = class, y = pathway, size = Count,color =logadjpv)) +
    geom_point(data = clProf.reshape.df , aes(x = class, y = pathway , size = Count), shape = 1,colour = "black")+
    scale_color_gradient2("-log(Adj.pv)",
                          # colours = c("blue","White","red"),
                          # values = c(0, -min(clProf.reshape.df$logadjpv)/(max(clProf.reshape.df$logadjpv)-min(clProf.reshape.df$logadjpv))   ,1),
                          midpoint = 0,
                          low = "blue", mid = "White",high = "red", space = "Lab" ,
                          
                          breaks= brs ,
                          labels= labels,
                          
                          # labels = abs(seq( floor(min(clProf.reshape.df$logadjpv))  , floor(max(clProf.reshape.df$logadjpv)) ,by=2)),
                          limits=c(min(clProf.reshape.df$logadjpv),max(clProf.reshape.df$logadjpv)),
                          guide=guide_colorbar(nbin = 500, barheight = 15))+
    theme(legend.text=element_text(size=rel(0.9)), axis.text.y = element_text(size=rel(1)), axis.text.x = element_text(size=rel(1.2)) )+
    ylab('')+xlab('')
  
  
  
  
  
  # if(nrow(downdata)!=0){
  #   p <- ggplot()+
  #     geom_point(data = downdata , aes(x = gender, y = as.character(AdjPvalu), size = Count,color = -log10(AdjPvalu))) + 
  #     geom_point(data = downdata , aes(x = gender, y = as.character(AdjPvalu) , size = Count), shape = 1,colour = "black")+
  #     # scale_x_discrete(labels = pathway)+
  #     scale_color_gradient2("DOWN\n-log(Adj.pv)", 
  #                           midpoint = -log10(downdata$AdjPvalu[1])/2, 
  #                           low = "white", mid = "lightblue", high = "blue", space = "Lab" ,
  #                           breaks=seq(1,-log10(downdata$AdjPvalu[1]),by=2),
  #                           limits=c(1, -log10(downdata$AdjPvalu[1])),
  #                           guide=guide_colorbar(nbin = 500, barheight = 7))+ new_scale_color()
  # }else{p <- ggplot()}
  #  
  #    
  #    # ggplot()+
  #    #   +  geom_point(data = clProf.reshape.df , aes(x = gender, y = pathway, size = Count,color = -log10(AdjPvalu))) + 
  #    #   + geom_point(data = clProf.reshape.df , aes(x = gender, y = pathway , size = Count), shape = 1,colour = "black")  
  #    
  #    
  #      
  # if(nrow(updata)!=0){
  #   p <- p +
  #     geom_point(data = updata , aes(x = gender, y = as.character(AdjPvalu), size = Count, color = -log10(AdjPvalu)))+ 
  #     geom_point(data = updata , aes(x = gender, y = as.character(AdjPvalu), size = Count), shape = 1,colour = "black")+
  #     scale_color_gradient2("UP\n-log(Adj.pv)", 
  #                           midpoint = -log10(updata$AdjPvalu[1])/2, 
  #                           low = "White", mid = "red",high = "darkred", space = "Lab" ,
  #                           breaks=seq(1,-log10(updata$AdjPvalu[1]),by=1),
  #                           limits=c(1, -log10(updata$AdjPvalu[1])),
  #                           guide=guide_colorbar(nbin = 500, barheight = 7))+
  #     theme(legend.text=element_text(size=rel(0.7)))+
  #     ylab('')+xlab('')
  # }
  
  
  # p <- ggplot()+
  #   geom_point(data = updata, aes(x = gender, y = pathway, size = Count,color = -log10(AdjPvalu))) +
  #   geom_point(data = updata, aes(x = gender, y = pathway, size = Count), shape = 1,colour = "black")+
  #   scale_color_gradient2("UP\n-log(Adj.pv)",
  #                         midpoint = -log10(updata$AdjPvalu[1])/2,
  #                         low = "White", mid = "red",high = "darkred", space = "Lab" ,
  #                         breaks=seq(1,-log10(updata$AdjPvalu[1]),by=1),
  #                         limits=c(1, -log10(updata$AdjPvalu[1])),
  #                         guide=guide_colorbar(nbin = 500, barheight = 7))+
  #   new_scale_color() +
  #   geom_point(data = downdata , aes(x = gender, y = pathway, size = Count, color = -log10(AdjPvalu)))+
  #   geom_point(data = downdata , aes(x = gender, y = pathway, size = Count), shape = 1,colour = "black")+
  #   scale_color_gradient2("DOWN\n-log(Adj.pv)",
  #                         midpoint = -log10(downdata$AdjPvalu[1])/2,
  #                         low = "white", mid = "lightblue", high = "blue", space = "Lab" ,
  #                         breaks=seq(1,-log10(downdata$AdjPvalu[1]),by=1),
  #                         limits=c(1, -log10(downdata$AdjPvalu[1])),
  #                         guide=guide_colorbar(nbin = 500, barheight = 7, reverse = TRUE))+
  #   # guides(color = guide_colourbar(reverse = TRUE))+
  #   theme(legend.text=element_text(size=rel(0.7)), axis.text.y = element_text(size=rel(1.2)) )+
  #   ylab('')+xlab('')
  
  
  
  
  # p <- ggplot()+
  #     geom_point(data = updata, aes(x = gender, y = pathway, size = Count,color = -log10(AdjPvalu))) +
  #     geom_point(data = updata, aes(x = gender, y = pathway, size = Count), shape = 1,colour = "black")+
  #     scale_color_gradient2("UP\n-log(Adj.pv)",
  #                         midpoint = -log10(updata$AdjPvalu[1])/2,
  #                         low = "White", mid = "red",high = "darkred", space = "Lab" ,
  #                         breaks=seq(1,-log10(updata$AdjPvalu[1]),by=1),
  #                         limits=c(1, -log10(updata$AdjPvalu[1])),
  #                         guide=guide_colorbar(nbin = 500, barheight = 7))+
  #     new_scale_color() +
  #     geom_point(data = downdata , aes(x = gender, y = pathway, size = Count, color = -log10(AdjPvalu)))+
  #     geom_point(data = downdata , aes(x = gender, y = pathway, size = Count), shape = 1,colour = "black")+
  #     scale_color_gradient2("DOWN\n-log(Adj.pv)",
  #                         midpoint = -log10(downdata$AdjPvalu[1])/2,
  #                         low = "white", mid = "lightblue", high = "blue", space = "Lab" ,
  #                         breaks=seq(1,-log10(downdata$AdjPvalu[1]),by=1),
  #                         limits=c(1, -log10(downdata$AdjPvalu[1])),
  #                         guide=guide_colorbar(nbin = 500, barheight = 7, reverse = TRUE))+
  #     # guides(color = guide_colourbar(reverse = TRUE))+
  #     theme(legend.text=element_text(size=rel(0.7)), axis.text.y = element_text(size=rel(1.2)) )+
  #     ylab('')+xlab('')
  
  
  if(split){p = p+ facet_wrap(~source, ncol = length(unique(clProf.reshape.df$source)))+ theme(strip.text = element_text(size = 14))}
  if(is.null(filename)){
    p
  }else{
    jpeg(filename, width = wid, height = height, res = 150)
    print(p)
    dev.off()
  }
  
}










