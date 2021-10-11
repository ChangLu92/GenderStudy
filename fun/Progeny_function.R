# get progeny score by top n gene
getprogenyScore <- function(pvmat = prog_pv,
                            top = 100,
                            zscore = TRUE,
                            perm.by = 'sample'){
  prog_mat_p <- prog_w[2:15]
  prog_mat_p <- sapply(1:14,function(x) {
    a = prog_mat_p[,x]
    a[order(pvmat[,x+1])[top:nrow(pvmat)]]=0
    return(a)
  })
  prog_mat_p = data.frame(prog_mat_p, row.names = prog_w$gene)
  colnames(prog_mat_p) <-  colnames(prog_w)[2:15]
  prog_mat_p <- prog_mat_p %>%
    rownames_to_column(var = 'GeneID')
  prog_mat_p[is.na(prog_mat_p)] = 0
  
  if(perm.by == 'gene'){
    progenyScore_ma = progenyPerm(GES_ma_df, prog_mat_p, k = 10000, z_scores = zscore, get_nulldist = FALSE) %>% t() 
    progenyScore_fe = progenyPerm(GES_fe_df, prog_mat_p, k = 10000, z_scores = zscore, get_nulldist = FALSE) %>% t()
  }else if(perm.by == 'sample'){
    progenyScore_ma = progenyPermSample(GES_df=GES_ma_df, prog_mat_p=prog_mat_p, dnull=dnullma, z_scores = zscore) 
    progenyScore_fe = progenyPermSample(GES_df=GES_fe_df, prog_mat_p=prog_mat_p, dnull=dnullfe, z_scores = zscore) 
  }
  
  
  # progenyScore_ma <- progeny(GES_ma_mat, scale=TRUE, organism= "Human", top = top, perm = 10000) %>% t()
  # progenyScore_fe <- progeny(GES_fe_mat, scale=TRUE, organism= "Human", top = top, perm = 10000) %>% t()
  progenyScore = merge(progenyScore_ma,progenyScore_fe, by="row.names")
  colnames(progenyScore) <- c('Pathway','male','female')
  return(list(prog_mat_p,progenyScore))
}



# get progeny score by adjusted pv cut off
getprogenyScorePV <- function(pvmat = prog_pv,
                              logpv = 1,
                              zscore = FALSE,
                              perm.by = 'sample',
                              k = 1000,
                              dnullma = NULL,
                              dnullfe = NULL){
  prog_mat_p <- prog_w[2:15]
  prog_mat_p[pvmat[,2:15]> 10^(-logpv)] = 0
  rownames(prog_mat_p) <-  prog_w$gene
  prog_mat_p = prog_mat_p[rowSums(prog_mat_p, na.rm = TRUE)!=0,]
  prog_mat_p = cbind(rownames(prog_mat_p),prog_mat_p)
  colnames(prog_mat_p)[1] = 'GeneID'
  prog_mat_p[is.na(prog_mat_p)] = 0
  
  if(perm.by == 'gene'){
    print(perm.by)
    progenyScore_ma = progenyPerm(GES_ma_df, prog_mat_p, k = k, z_scores = zscore, get_nulldist = FALSE) %>% t()
    progenyScore_fe = progenyPerm(GES_fe_df, prog_mat_p, k = k, z_scores = zscore, get_nulldist = FALSE) %>% t()
  }else if(perm.by == 'sample'){
    print(perm.by)
    progenyScore_ma = progenyPermSample(GES_df=GES_ma_df, prog_mat_p=prog_mat_p, dnull=dnullma, z_scores = zscore)
    progenyScore_fe = progenyPermSample(GES_df=GES_fe_df, prog_mat_p=prog_mat_p, dnull=dnullfe, z_scores = zscore)
  }
  
  progenyScore = merge(progenyScore_ma,progenyScore_fe, by="row.names")
  colnames(progenyScore) <- c('Pathway','male','female')
  return(list(prog_mat_p,progenyScore))
}


# Sample permutaion 
progenyPermSample <- function(GES_df, prog_mat_p, dnull, z_scores = TRUE){
  current_weights <- prog_mat_p
  names(GES_df)[1] <- "ID"
  names(current_weights)[1] <- "ID"
  common_ids <- merge(GES_df, current_weights, by = "ID")
  common_ids <- as.character(common_ids$ID)
  
  row.names(GES_df) <- GES_df$ID 
  GES_df <- as.data.frame(GES_df[common_ids,-1])
  
  row.names(current_weights) <- current_weights$ID
  current_weights <- as.data.frame(current_weights[common_ids,-1])
  current_mat <- as.matrix(GES_df)
  current_weights <- t(current_weights)
  scores <- as.data.frame(current_weights %*% current_mat)
  
  currentdnull <- dnull[common_ids,]
  null_dist_scores =  current_weights %*% currentdnull
  
  
  if(z_scores) {
    scores$mean <- apply(null_dist_scores,1,mean)
    scores$sd <- apply(null_dist_scores,1,sd)
    resListCurrent <- (scores[,1]-scores[,2])/scores[,3]
    names(resListCurrent) <- names(prog_mat_p[,-1])
  } else {
    for(j in seq(1, length(prog_mat_p[,-1]))) {
      ecdf_function <- ecdf(null_dist_scores[j,])
      scores[j,1] <- ecdf_function(scores[j,1])
    }
    score_probas <- scores*2-1
    resListCurrent <- score_probas[,1]
    names(resListCurrent) <- names(prog_mat_p[,-1])
  }
  resDf <- as.data.frame(resListCurrent)
  return(resDf)
}


## visualize the stability of NES by line chart
linePlotVis <- function(progScorelist, filename, toplist, xlabel ){
  pathwayname = progScorelist[[1]]
  PS_ma = sapply(seq(from = 2, to = length(progScorelist)-1, by = 3), function(x) progScorelist[[x]])
  PS_fe = sapply(seq(from = 3, to = length(progScorelist), by = 3), function(x) progScorelist[[x]]) 
  for (i in 1:14) {
    aa = data.frame(toplist, PS = PS_ma[i,], gender = rep('male', times = length(toplist)))
    bb = data.frame(toplist, PS = PS_fe[i,], gender = rep('female', times = length(toplist)))
    df = rbind(aa,bb) %>%
      add_column(class = factor(rep(pathwayname[i], times = length(toplist)*2)))
    
    if(i==1){
      df_all = df
    }else{
      df_all = rbind(df_all,df)
    }
    
  }
  
  jpeg(filename, height = 220, width = 2500, res = 150)
  print(ggplot(df_all, aes(x=toplist, y = PS,colour = gender)) + 
    geom_line() + geom_point()+
    facet_wrap(~class,ncol=14)+
    scale_y_continuous(breaks= c(-4,-1.96,0,1.96,4))+
    geom_hline(aes(yintercept = 0),linetype="dashed")+
    geom_hline(aes(yintercept = 1.96),linetype="dashed")+
    geom_hline(aes(yintercept = -1.96),linetype="dashed")+
    # ylim(-5, +5)+
    xlab(xlabel)+
    # xlab("top")+
    ylab("NES"))
  dev.off()
  
}


# visualize the progeny matrix weight and GES by vocalno plot
ProgenyvocalnoPlot <- function(tmat,       
                               filename, 
                               xname = 'logFC',
                               yname = "adj.P.Val",
                               weightname = 'Progeny.weight',
                               showoppolabel = TRUE,
                               labelby = 'Progeny.weight',
                               top = 20){
  
  tmat <- tmat %>% column_to_rownames(var = 'ID')
  idxop = tmat[,weightname] * tmat[,xname] < 0
  idxsame = tmat[,weightname] * tmat[,xname] > 0
  
  if(labelby == yname){
    dfop = tmat[idxop,] %>% arrange_(labelby)
    dfsa = tmat[idxsame,] %>% arrange_(labelby)
  }else{
    dfop = tmat[idxop,] %>% arrange_(paste0("desc(abs(",labelby,"))"))
    dfsa = tmat[idxsame,] %>% arrange_(paste0("desc(abs(",labelby,"))"))
  }
  
  
  tmat$delabel <- NA
  
  if(showoppolabel){
    tmat[rownames(dfop)[1:min(c(nrow(dfop),top))], "delabel"] <- rownames(dfop)[1:min(c(nrow(dfop),top))]
  }else{
    tmat[rownames(dfsa)[1:min(c(nrow(dfsa),top))], "delabel"]<- rownames(dfsa)[1:min(c(nrow(dfsa),top))]
  }
  
  
  
  jpeg(filename, width = 1000, height = 1000, res = 180)
  p = ggplot(data=tmat, aes_string(x=xname, y=sprintf("-log10(%s)", yname)))
  print(p +
          # geom_point(aes_string(color = weightname,size = sprintf("abs(%s)",weightname)) ) +
          geom_point(aes_string(color = weightname) , size = 3) +
          scale_color_gradient2(weightname, 
                                midpoint = 0, 
                                low = "darkblue", mid = "white",high = "darkred", space = "Lab" ,
                                limits=c(min(tmat[,weightname]), max(tmat[,weightname])))+
          geom_label_repel(aes_string(label = "delabel", color = weightname),fontface = "bold" )+
          ggtitle(gsub(".*/ *(.*?) *.jpg.*", "\\1",filename)) )
  dev.off()
  
}

## compare the NES of 2 groups (female vs males) by scatter plot
ProgenyScatterPlot <- function(tmat,       
                               filename, 
                               xname = 'male',
                               yname = "female",
                               valuename = 'logFC',
                               weightname = 'Progeny.weight',
                               labelby = 'Progeny.weight',
                               top = 20){
  library(ggnewscale)
  if(weightname != labelby) {
    tmat <- tmat %>% tibble::column_to_rownames(var = 'ID') %>% select_(weightname, labelby, paste0(valuename,'.',xname),paste0(valuename,'.',yname))
  }else{
    tmat <- tmat %>% tibble::column_to_rownames(var = 'ID') %>% select_(weightname, paste0(valuename,'.',xname),paste0(valuename,'.',yname))
  }
  
  
  # idxop = tmat[,weightname] * tmat[,xname] < 0
  # idxsame = tmat[,weightname] * tmat[,xname] > 0
  # dfop = tmat[idxop,] %>% arrange_(paste0("desc(abs(",labelby,"))"))
  # dfsa = tmat[idxsame,] %>% arrange_(paste0("desc(abs(",labelby,"))"))
  
  
  df = tmat %>% arrange_(paste0("desc(abs(",labelby,"))"))
  
  
  
  tmat$delabel <- NA
  tmat[rownames(df)[1:top], "delabel"] <- rownames(df)[1:top]
  
  # if(showoppolabel){
  #   tmat[rownames(dfop)[1:min(c(nrow(dfop),top))], "delabel"] <- rownames(dfop)[1:min(c(nrow(dfop),top))]
  # }else{
  #   tmat[rownames(dfsa)[1:min(c(nrow(dfsa),top))], "delabel"]<- rownames(dfsa)[1:min(c(nrow(dfsa),top))]
  # }
  # 
  
  
  jpeg(filename, width = 1000, height = 1000, res = 150)
  p = ggplot(data=tmat, aes_string(x=paste0(valuename,'.',xname), y=paste0(valuename,'.',yname)))
  
  print(p +
          geom_point(aes_string(color = weightname),size = 2) + 
          scale_color_gradient2(weightname, 
                                midpoint = 0, 
                                low = "darkblue", mid = "white",high = "darkred", space = "Lab" ,
                                limits=c(min(tmat[,weightname]), max(tmat[,weightname])))+
          # new_scale_color() +
          geom_label_repel(aes_string(label = "delabel", color = weightname),fontface = "bold" )+
          theme(
            axis.title = element_text(size = 16),
            axis.text = element_text(size = 15),
            legend.title = element_text(size = 15),
            legend.text = element_text(size = 15),
            legend.position="top")+
        # scale_colour_steps2(n.breaks = 3,
        #                     low = "darkblue",
        #                     mid = "white",
        #                     high = "darkred")+
        ggtitle(gsub(".*/ *(.*?) *.jpg.*", "\\1",filename))
  )
  
  dev.off()
  
}





