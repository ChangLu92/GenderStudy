library(CARNIVAL)
library(OmnipathR)
library(visNetwork)
library(dplyr)
source("./fun/assignPROGENyScores.r")
source("./fun/generateTFList.r")
# source('support_networks.r')


# get sif dataframe for CARNICAL
getSIF <- function(omniR){
  # signed and directed
  omnipath_sd <- omniR %>% dplyr::filter(consensus_direction == 1 &
                                           (consensus_stimulation == 1 | 
                                              consensus_inhibition == 1))
  # changing 0/1 criteria in consensus_stimulation/inhibition to -1/1
  omnipath_sd$consensus_stimulation[which( omnipath_sd$consensus_stimulation == 0)] = -1
  omnipath_sd$consensus_inhibition[which( omnipath_sd$consensus_inhibition == 1)] = -1
  omnipath_sd$consensus_inhibition[which( omnipath_sd$consensus_inhibition == 0)] = 1
  # check consistency on consensus sign and select only those in a SIF format
  sif <- omnipath_sd[,c('source_genesymbol', 'consensus_stimulation', 'consensus_inhibition', 'target_genesymbol')] %>%
    dplyr::filter(consensus_stimulation==consensus_inhibition) %>%
    unique.data.frame()
  sif$consensus_stimulation <- NULL
  colnames(sif) <- c('source', 'interaction', 'target')
  # remove complexes
  sif$source <- gsub(":", "_", sif$source)
  sif$target <- gsub(":", "_", sif$target)
  return(sif)
}


# #save SIF
# write_tsv(sif, "./results/omnipath_carnival.tsv")

# # get initial nodes
# iniMTX = base::setdiff(sif$source, sif$target)
# iniciators = base::data.frame(base::matrix(data = NaN, nrow = 1, ncol = length(iniMTX)), stringsAsFactors = F)
# colnames(iniciators) = iniMTX

# get CARNIVAL result
# PathwayActivity_CARNIVALinput: progeny pathway activities
# tf_activities_carnival: tf activities
getCARNIVALresult <- function(PathwayActivity_CARNIVALinput, tf_activities_carnival,sif){
  progenylist = assignPROGENyScores(progeny = t(PathwayActivity_CARNIVALinput), 
                                    progenyMembers = progenyMembers, 
                                    id = "gene", 
                                    access_idx = 1)
  tfList = as.data.frame(t(tf_activities_carnival))
  # run carnival
  carnival_result = runCARNIVAL(inputObj= iniciators,
                                measObj = tfList, 
                                netObj = sif, 
                                weightObj = as.matrix(progenylist[[1]]), 
                                solverPath = "D:/Newfolder/IBM/LOG/CPLEX_Studio201/cplex/bin/x64_win64/cplex.exe",
                                solver = "cplex",
                                keepLPFiles = FALSE)
                                # dir_name = "./results")
  #transoform to data.frame
  carnival_result$weightedSIF <- data.frame(carnival_result$weightedSIF, stringsAsFactors = F)
  carnival_result$weightedSIF$Sign <- as.numeric(carnival_result$weightedSIF$Sign)
  carnival_result$weightedSIF$Weight <- as.numeric(carnival_result$weightedSIF$Weight)
  carnival_result$nodesAttributes <- data.frame(carnival_result$nodesAttributes, stringsAsFactors = F)
  carnival_result$nodesAttributes$ZeroAct <- as.numeric(carnival_result$nodesAttributes$ZeroAct)
  carnival_result$nodesAttributes$UpAct <- as.numeric(carnival_result$nodesAttributes$UpAct)
  carnival_result$nodesAttributes$DownAct <- as.numeric(carnival_result$nodesAttributes$DownAct)
  carnival_result$nodesAttributes$AvgAct <- as.numeric(carnival_result$nodesAttributes$AvgAct)
  return(carnival_result)
}


CARNIVAL_Pipeline = function(cases, idx, sif){
  TF_activity_all <- switch(
    cases[3],
    GenePerm= TFactivity_Carnival_genePerm,
    samplePerm= TFactivity_Carnival_sampPerm
  )
  
  progScorelist <- switch(
    cases[3],
    GenePerm= progScorelist_GenPerm,
    samplePerm= progScorelist_SamPerm
  )
  
  progenyScore = sapply(seq(3*idx-1,3*idx), function(x) as.numeric(progScorelist[[x]]))  %>% 
    as.data.frame() %>% 
    add_column(Pathway = progScorelist[[1]])
  colnames(progenyScore) <- c('male','female','Pathway')
  progenyScore$Pathway[progenyScore$Pathway %in%"JAK-STAT"] = "JAK.STAT"
  
  progenyScore_ma = sapply(seq(from = 2, to = length(progScorelist), by = 3), function(x) as.numeric(progScorelist[[x]]))
  rownames(progenyScore_ma) <- progScorelist[[1]]
  info_ma = apply(progenyScore_ma, 1, function(x) c(sd(x),min(x),max(x)))
  unstableidx_ma = info_ma[1,]>0.2 & info_ma[2,]*info_ma[3,]<0
  
  progenyScore_fe = sapply(seq(from = 3, to = length(progScorelist), by = 3), function(x) as.numeric(progScorelist[[x]]))
  rownames(progenyScore_fe) <- progScorelist[[1]]
  info_fe = apply(progenyScore_fe, 1, function(x) c(sd(x),min(x),max(x)))
  unstableidx_fe = info_fe[1,]>0.2 & info_fe[2,]*info_fe[3,]<0
  
  progenyScoreR = progenyScore
  progenyScoreR$male[unstableidx_ma] <- 0
  progenyScoreR$female[unstableidx_fe] <- 0
  
  progenyScoreA<- switch(
    cases[2],
    normal= progenyScore,
    remo_unstable= progenyScoreR)
  
  
  PathwayActivity_CARNIVALinput <- switch(
    cases[1],
    male= progenyScoreA %>% tibble::column_to_rownames(var = 'Pathway')%>% dplyr::select(male),
    female= progenyScoreA %>% tibble::column_to_rownames(var = 'Pathway')%>% dplyr::select(female)
  )
  
  tf_activities_carnival = switch(
    cases[1],
    male= TF_activity_all %>% tibble::column_to_rownames(var = 'TF')%>% dplyr::filter(labels == 'Male specific TF') %>% dplyr::select(NES.male),
    female= TF_activity_all %>% tibble::column_to_rownames(var = 'TF')%>% dplyr::filter(labels == 'Female specific TF') %>% dplyr::select(NES.female)
  )
  carnival_result <- getCARNIVALresult(PathwayActivity_CARNIVALinput, tf_activities_carnival,sif)
  return(carnival_result)
}


draw_CYTOSCAPE_CARNIVAL <- function(evis,nvis,title, style.name){
  library(RCy3)
  cytoscapePing ()
  cytoscapeVersionInfo ()
  ## edges
  colnames(evis) = c('source', "interaction", "target", "weight")
  evis$interaction[evis$interaction == 1] = "activation"
  evis$interaction[evis$interaction == -1] = "inhibition"
  
  ## nodes
  nvis = nvis[which(nvis$ZeroAct!=100),]
  nvis$ZeroAct = NULL
  nvis = nvis[which(nvis$Node%in%union(evis$source, evis$target)),]
  nvis$label = nvis$Node
  
  colnames(nvis) = c("id", "UpAct", "DownAct", "AvgAct", "group", "label")
  nvis$group = replace(nvis$group, nvis$group=='T', 'TFs')
  nvis$group = replace(nvis$group, nvis$group=='S', 'Perturbed')
  nvis$group = replace(nvis$group, nvis$group=='', 'Protein')
  nvis$AvgAct = as.numeric(nvis$AvgAct) 
  createNetworkFromDataFrames(nvis,evis, title=title)
  setVisualStyle(style.name)
}






getCARNIVALdegree <- function(sift){
  count_degree = sift %>% degree_count()
  p <- data.frame(Node = count_degree$node,
                  p = as.numeric(count_degree$total_count) / nrow(sift), 
                  total_count = count_degree$total_count, 
                  row.names =count_degree$node) 
  return(p)
}




# saveRDS(carnival_result,paste0("./results/",filename,"_carnival_result.rds") )















