###########################################
#  - TFMir  Project                       #
#  - TFMir  main function                 #
#  - Main function to be called           #
#  - 2014-10-1                            #
#  - Copyright: Mohamed Hamed             #
###########################################
start.time <- Sys.time()
#R programming environments:
#- R studio IDE
#- R version 2.12.0 (2010-10-15)
#-Platform: x86_64-apple-darwin9.8.0/x86_64 (64-bit)
#-locale:   [1] C/en_US.UTF-8/C/C/C/C


setwd("../backend")

## =============================================================================
## Initialize calling required packages and install them if they are not and log
## =============================================================================

source("loadpackages.R")


## ==================================================
## read the configuration file and get the parameters 
## ==================================================

#source("readconfig.R")
#config

## ===================================================
## Load the statistics ,color venn, and grapth scripts
## ===================================================
source("statistics.R")
source("ColorVenn.R")

#source("graph.R")
#source("initialize.R")

  

## =============================================================================================================================
## Intialize and Read all files in config file (tranmir db , mirna disease association file and mirna function association file)
## =============================================================================================================================

#source("initialize.R")
#source("initialize_tissue.R")
#getMIRnaCategory("hsa-mir-212","function")
#getMIRnaCategory("hsa-mir-212","disease")

## ===============================
## call the main function of TFMir
## ===============================

readInput=function(path)
{
  molecule.df = read.delim(path, header=FALSE)
  molecule.df=molecule.df[!duplicated (molecule.df),]
  molecule.input=unique(as.character(unlist(molecule.df[1]))) 
  return(molecule.input)
}
end.time <- Sys.time()
time.taken.packages <- end.time - start.time
time.taken.packages

start.time <- Sys.time()
#speciesp = "Human"
TFMir =function(tf.path, mirna.path,enrich.pval,pval.cutoff=0.05,pp.cutoff,evidence,species,disease="",processes="", tissue="", output.path)
{
  process= unlist(strsplit(processes, split=":"))
  speciesp <<- species;
  source("readconfig.R")
  config
  #if(pp.cutoff != 0.8){
    pp.cutoffp <<-pp.cutoff;
    source("Database_construct.R")
  #}
  #config
  source("graph.R")
  source("initialize.R")

  ## ========================
  ## log the input parameters 
  ## ========================  
  writeToLog(tf.path)
  writeToLog(mirna.path)
  writeToLog(pval.cutoff)
  writeToLog(disease)
  writeToLog(evidence)
  writeToLog(output.path)
  writeToLog("=====================================")
  
  if(tolower(evidence)=="both") { evidence=c("Experimental","Predicted") }

  ## ==================================================================
  ## read the input files and test automatically which scenario will be
  ## ==================================================================
  tfs.df=data.frame("gene"=character(),"gene.reg"=numeric())
  if(! is.na(tf.path) & tf.path !="")
  {
    tfs.df = read.delim(tf.path, header=FALSE)    
  }else
  {
    mirnas.input=unique(tolower(readInput(mirna.path)))
    config$pval.cutoffgene.targets_regulators.fromMiRNA.inputlist = enrich.pval
    tf.pval=as.double(config$pval.cutoffgene.targets_regulators.fromMiRNA.inputlist)
    tfs.regulatorsofMiRNA= getTFsRegulatorsofMiRNAs(mirnas.input,tf.pval,evidence)
    tfs.targetsofMiRNA= getTFsTargetsofMiRNAs(mirnas.input,tf.pval,evidence)
    tfs.list=unique(c(tfs.targetsofMiRNA,tfs.regulatorsofMiRNA))
    if(length(tfs.list) > 0){
      tfs.df=data.frame("gene"=tfs.list,"gene.reg"=0) }
  }
  names(tfs.df)=c("gene","gene.reg")
  tfs.df=tfs.df[!duplicated (tfs.df),]
  tfs.input=unique(as.character(unlist(tfs.df$gene)))
  #printGeneEntrezIDsMap(tfs.input,output.path)
  #print(tfs.input)

  
  
  
  
  mirnas.df=data.frame("mirna"=character(),"mirna.reg"=numeric())
  if(! is.na(mirna.path) & mirna.path !="")
  {
    mirnas.df = read.delim(mirna.path, header=FALSE)    
  }else
  {
    config$pval.cutoffmirna.regulators_targets.fromTFS.inputlist = enrich.pval
    mirna.pval=as.double(config$pval.cutoffmirna.regulators_targets.fromTFS.inputlist)
    tfstargets.mirna= getMiRNAsTargetsofTFs(tfs.input,mirna.pval,evidence)
    tfsregulators.mirna= getMiRNAsRegulatorsofTFs(tfs.input,mirna.pval,evidence)
    mirnas.list=unique(c(tfstargets.mirna,tfsregulators.mirna))
    if(length(mirnas.list) > 0){
    mirnas.df=data.frame("mirna"=mirnas.list,"mirna.reg"=0) }
  }
  names(mirnas.df)=c("mirna","mirna.reg")
  mirnas.df$mirna=tolower(mirnas.df$mirna)
  mirnas.df=mirnas.df[!duplicated (mirnas.df),]
  mirnas.input=unique(as.character(unlist(mirnas.df$mirna)))  
  #print(mirnas.input)
  
  ## ==================================
  ## get the five kinds of interactions 
  ## ==================================
  tf.mirna.res=getInteractions(category="tf-mirna",reg.input=tfs.input,target.input=mirnas.input,process=process,tissue=tissue,disease=disease,evidence=evidence,output.path=output.path,pval.cutoff=pval.cutoff)    
  tf.gene.res=getInteractions(category="tf-gene",reg.input=tfs.input,target.input=tfs.input,process=process,tissue=tissue,disease=disease,evidence=evidence,output.path=output.path,pval.cutoff=pval.cutoff)  
  #mirna.gene.res=tf.mirna.res[tf.mirna.res$category=="mirna-mirna",]  # return initial empty structure
  mirna.gene.res=getInteractions(category="mirna-gene",reg.input=mirnas.input,target.input=tfs.input,process=process,tissue=tissue,disease=disease,evidence=evidence,output.path=output.path,pval.cutoff=pval.cutoff)
  gene.gene.res=getInteractions(category="gene-gene",reg.input=tfs.input,target.input=tfs.input,process=process,tissue=tissue,disease=disease,evidence=evidence,output.path=output.path,pval.cutoff=pval.cutoff)
  mirna.mirna.res=tf.mirna.res[tf.mirna.res$category=="mirna-mirna",]  # return initial empty structure
  if("Predicted" %in% evidence)   ## cause mirna-mirna interacctions are only predictions
    mirna.mirna.res=getInteractions(category="mirna-mirna",reg.input=mirnas.input,target.input=mirnas.input,process=process,tissue=tissue,disease=disease,evidence=evidence,output.path=output.path,pval.cutoff=pval.cutoff)
  
  ## ======================================================================================================================
  ## Combine these interactions and get those related to disease only (disease speccific network) (if disease is specified)
  ## ======================================================================================================================
  input=list( tf.genes=names(tf.gene.res), mirna.genes=names(mirna.gene.res),tf.mirna=names(tf.mirna.res),mirna.mirna=names(mirna.mirna.res),gene.gene=names(gene.gene.res))
  columns=Reduce(intersect,input)
  all.res=rbind(tf.mirna.res[,columns],mirna.mirna.res[,columns],mirna.gene.res[,columns],tf.gene.res[,columns],gene.gene.res[,columns])

  names(mirnas.df)=c("node","regulation")
  names(tfs.df)=c("node","regulation")
  nodes.input=rbind(mirnas.df,tfs.df)
  names(nodes.input)=c("target","target.reg")
  all.res=merge(all.res,nodes.input,by="target")
  names(nodes.input)=c("regulator","regulator.reg")
  all.res=merge(all.res,nodes.input,by="regulator")
  
  all.res.disease=all.res[all.res$is_regulator_in_disease==TRUE | all.res$is_target_in_disease==TRUE,]
  #write("finished", file=file.path(output.path,"finished.txt"),append=F,sep="\n")
  
  ## =================================================================================================
  ## Combine the interactions and get those related to tissue only (tissue speccific network) 
  ## =================================================================================================
  set_mirna = all.res[((all.res$category %in% "tf-mirna" & all.res$is_regulator_in_tissue==TRUE) | (all.res$category %in% "mirna-gene" & all.res$is_target_in_tissue==TRUE)),]
  all.res.tissue=all.res[(all.res$is_regulator_in_tissue==TRUE & all.res$is_target_in_tissue==TRUE),] 
  all.res.tissue = rbind(all.res.tissue,set_mirna)
  #print(c(dim(set_mirna)[1], dim(all.res.tissue)[1]))
  #print(all.res.tissue)
  
 ## =================================================================================================
 ## Combine the interactions and get those related to tissue only (tissue speccific network) 
 ## =================================================================================================
 all.res.process=all.res[all.res$is_regulator_in_process==TRUE | all.res$is_target_in_process==TRUE,]
#print(all.res.process)

## =================================================================================================
## Combine the interactions and get those related to disease and process only  
## =================================================================================================
all.res.disease.process=all.res[(all.res$is_regulator_in_process==TRUE & all.res$is_target_in_disease==TRUE)
                                | (all.res$is_regulator_in_disease==TRUE & all.res$is_target_in_process==TRUE),] 
#print(all.res.disease.process)

## =================================================================================================
## Combine the interactions and get those related to tissue and process only  
## =================================================================================================
all.res.tissue.process=all.res[(all.res$is_regulator_in_process==TRUE & all.res$is_target_in_tissue==TRUE)
                                | (all.res$is_regulator_in_tissue==TRUE & all.res$is_target_in_process==TRUE),] 
#print(all.res.tissue.process)

## =================================================================================================
## Combine the interactions and get those related to tissue and process only  
## =================================================================================================
all.res.disease.tissue=all.res[(all.res$is_regulator_in_disease==TRUE & all.res$is_target_in_tissue==TRUE)
                               | (all.res$is_regulator_in_tissue==TRUE & all.res$is_target_in_disease==TRUE),] 
#print(all.res.tissue.process)

  ## ===================================================================
  ## get tissue specific genes for disease and tissue and full networks.
  ## ===================================================================
  #tissue = "Brain - Cerebellum"
  all.res.tissue.specific = list()
  all.res.disease.tissue.specific = list()
  all.res.normal.tissue.specific = list()
  all.res.process.tissue.specific = list()
  all.res.tissue.process.tissue.specific = list()
  all.res.disease.process.tissue.specific = list()
  tissue.specific.gene = list()
  tissue.specific.gene = getTissue_specificGenes(tissue)
  names(tissue.specific.gene) = "regulator"
  if(dim(tissue.specific.gene)[1] > 0){
    all.res.tissue.specific = union(intersect(all.res$regulator,tissue.specific.gene$regulator),intersect(all.res$target,tissue.specific.gene$regulator))
    all.res.disease.tissue.specific = union(intersect(all.res.disease$regulator,tissue.specific.gene$regulator),intersect(all.res.disease$target,tissue.specific.gene$regulator))
    all.res.normal.tissue.specific = union(intersect(all.res.tissue$regulator,tissue.specific.gene$regulator),intersect(all.res.tissue$target,tissue.specific.gene$regulator))
    all.res.process.tissue.specific = union(intersect(all.res.process$regulator,tissue.specific.gene$regulator),intersect(all.res.process$target,tissue.specific.gene$regulator))
    all.res.disease.process.tissue.specific = union(intersect(all.res.disease.process$regulator,tissue.specific.gene$regulator),intersect(all.res.disease.process$target,tissue.specific.gene$regulator))
    all.res.disease.tissue.tissue.specific = union(intersect(all.res.disease.tissue$regulator,tissue.specific.gene$regulator),intersect(all.res.disease.tissue$target,tissue.specific.gene$regulator))
  }
  if(length(all.res.tissue.specific) > 0){
    all.res.tissue.specific = toString(all.res.tissue.specific)
  }

  if(length(all.res.disease.tissue.specific) > 0){
    all.res.disease.tissue.specific = toString(all.res.disease.tissue.specific)
  }
  if(length(all.res.normal.tissue.specific) > 0){
    all.res.normal.tissue.specific = toString(all.res.normal.tissue.specific)
  }
  if(length(all.res.process.tissue.specific) > 0){
    all.res.process.tissue.specific = toString(all.res.process.tissue.specific)
  }
  if(length(all.res.disease.process.tissue.specific) > 0){
    all.res.disease.process.tissue.specific = toString(all.res.disease.process.tissue.specific)
  }
  if(length(all.res.disease.tissue.tissue.specific) > 0){
    all.res.disease.tissue.tissue.specific = toString(all.res.disease.tissue.tissue.specific)
  }

  if(dim(all.res)[1] > 0)
    exportNetworkProperties (all.res,file.path(output.path,"all"), disease,pval.cutoff,all.res.tissue.specific)
  
  if(dim(all.res.disease)[1] > 0)
    exportNetworkProperties (all.res.disease,file.path(output.path,"disease"),disease,pval.cutoff,all.res.disease.tissue.specific)

  if(dim(all.res.tissue)[1] > 0)
    exportNetworkProperties (all.res.tissue,file.path(output.path,"tissue"),disease,pval.cutoff,all.res.normal.tissue.specific)
  
  if(dim(all.res.process)[1] > 0)
    exportNetworkProperties (all.res.process,file.path(output.path,"process"),disease,pval.cutoff,all.res.process.tissue.specific)
  
  if(dim(all.res.disease.process)[1] > 0)
    exportNetworkProperties (all.res.disease.process,file.path(output.path,"disease_process"),disease,pval.cutoff,all.res.disease.process.tissue.specific)
  
  if(dim(all.res.tissue.process)[1] > 0)
    exportNetworkProperties (all.res.tissue.process,file.path(output.path,"tissue_process"),disease,pval.cutoff,all.res.tissue.process.tissue.specific)
  
  if(dim(all.res.disease.tissue)[1] > 0)
    exportNetworkProperties (all.res.disease.tissue,file.path(output.path,"disease_tissue"),disease,pval.cutoff,all.res.disease.tissue.tissue.specific)
write("finished", file=file.path(output.path,"finished.txt"),append=F,sep="\n")
}

ExportMotifs =function(net.path,output.path,evidence,species,random.method)
{ 
  speciesp <<- species;
  source("readconfig.R")
  config
  source("graph.R")
  source("initialize.R")
  motifs.FFL = c()
  print(net.path)
  print(getwd())
  net=read.delim(net.path, header=TRUE) 
  if(tolower(evidence)=="both") { evidence=c("Experimental","Predicted") }
  
  #### extract all putative Tf-mirna paris who share target genes
  tfmir.pairs=getPutativeTFmiRPairs(net)
  #print(dim(tfmir.pairs)[1])
  #### extract all significant TF - miRNA pairs
  #tfmir.pairs=getSignificantTFmiRpairs(tfmir.pairs,evidence)
  
  #if(dim(tfmir.pairs)[1] >0) 
  #{
    #### relax and message the TF mir pairs who have more than one target
    #tfmir.pairs=relaxMultipleTargetsForTFmiRPairs(tfmir.pairs)
    #print(dim(tfmir.pairs)[1])
    #### get motif type 1 : composite-FFL
    #motifs.composite= getMotifs.composite(tfmir.pairs,net)
    #print("End of motifs.composite:")
    #### get motif type 2 : TF-FFL
    #motifs.TF.FFL= getMotifs.TF.FFL(tfmir.pairs,net)
    #print("End of motifs.TF.FFL:")
    #### get motif type 3 : miRNA-FFL
    #motifs.miRNA.FFL= getMotifs.miRNA.FFL(tfmir.pairs,net)
    #print("End of motifs.miRNA.FFL:")
    #### get motif type 4 : Coregulation-FFL
    #motifs.coregulation= getMotifs.coregulation(tfmir.pairs,net)
    #print("End of motifs.coregulation:")
    #motifs.FFL=rbind(motifs.composite,motifs.TF.FFL,motifs.miRNA.FFL,motifs.coregulation)
    #if( dim(motifs)[1] > 0 )
    #{
      #motifs.ids=paste("motif",seq(1:dim(motifs)[1]), sep="" )
      #motifs=cbind(motifs.ids,motifs)
      #write.table(motifs,file=file.path(output.path,"motifs.txt"),quote=F,row.names=F,col.names=T,sep="\t")
    #}
  #}
  #print(motifs.FFL)
  if(dim(tfmir.pairs)[1] >0) 
  {
    print(dim(tfmir.pairs)[1])
    tfmir.pairs=relaxMultipleTargetsForTFmiRPairs(tfmir.pairs)
    print(dim(tfmir.pairs)[1])
    motifs.FFL = c()
    motifs.FFL= get3NodeMotifs.FFL(tfmir.pairs,net,random.method)
  }
    ###change the net.path from res.txt to res.4node.txt which is necessary for 4node-motifs 
    A = net.path
    A=strsplit(A,"/")
    AF = strsplit(A[[1]][length(A[[1]])],"\\.")
    AF[[1]][length(AF[[1]])]=paste(AF[[1]][length(AF[[1]])-1],".4node.txt",sep="")
    A[[1]][length(A[[1]])] = AF[[1]][length(AF[[1]])]
    B=c()
    B=A[[1]][1]
    for(i in 2:length(A[[1]])){
      B=paste(B,"/",A[[1]][i],sep="")
    }
    net.path.4node = B
    ####
    net=read.delim(net.path.4node, header=TRUE) 
    #### extract all putative Tf-mirna paris who share target genes
    tfmir.pairs.tf.mirna=getPutativeTFmiRPairs4Node_tf_mirna(net)
    tfmir.pairs.tf=getPutativeTFmiRPairs4Node_tf(net)
    tfmir.pairs.mirna=getPutativeTFmiRPairs4Node_mirna(net)
    tfmir.pairs.gene.gene=getPutativeTFmiRPairs4Node(net)
  
    tfmir.pairs = list(tfmir.pairs.tf.mirna,tfmir.pairs.tf,tfmir.pairs.mirna,tfmir.pairs.gene.gene)
    if(length(tfmir.pairs) >0) 
    {
      #print("Begin of motifs.composite:")
      motifs = c()
      motifs= get4NodeMotifs.FFL(tfmir.pairs,net,random.method)
      #print(dim(motifs)[1])
      #print("End of motifs.composite:")
      if( dim(motifs)[1] > 0 ){
        motifs.FFL = rbind(motifs.FFL,motifs)
      }
      #print(c("motifs.FFL", length(motifs.FFL)))
      if( length(motifs.FFL) > 0 )
      {
        motifs.ids=paste("motif",seq(1:dim(motifs.FFL)[1]), sep="")
        #print(motifs.ids)
        motifs.FFL=cbind(motifs.ids,motifs.FFL)
        write.table(motifs.FFL,file=file.path(output.path,"motifs.txt"),quote=F,row.names=F,col.names=T,sep="\t")
      }
      else{
        motifs.FFL = c()
        write.table(motifs.FFL,file=file.path(output.path,"motifs.txt"),quote=F,row.names=F,col.names=T,sep="\t")
      }
      write("finished motifs 4 nodes", file=file.path(output.path,"finished4Nodemotifs.txt"),append=F,sep="\n")  
    }  
}
Export4NodeMotifs =function(net.path,output.path,evidence,species)
{
  #species="Human"
  speciesp <<- species;
  source("readconfig.R")
  config
  #config
  source("graph.R")
  source("initialize.R")
  print(net.path)
  print(getwd())
  net=read.delim(net.path, header=TRUE) 
  if(tolower(evidence)=="both") { evidence=c("Experimental","Predicted") }
  
  #### extract all putative Tf-mirna paris who share target genes
  tfmir.pairs.tf.mirna=getPutativeTFmiRPairs4Node_tf_mirna(net)
  tfmir.pairs.tf=getPutativeTFmiRPairs4Node_tf(net)
  tfmir.pairs.mirna=getPutativeTFmiRPairs4Node_mirna(net)
  tfmir.pairs.gene.gene=getPutativeTFmiRPairs4Node(net)
  
  
  tfmir.pairs = list(tfmir.pairs.tf.mirna,tfmir.pairs.tf,tfmir.pairs.mirna,tfmir.pairs.gene.gene)
  if(length(tfmir.pairs) >0) 
  {
    motifs.FFL= get4NodeMotifs.FFL(tfmir.pairs,net)
    if( length(motifs.FFL) > 0 )
    {
      motifs.ids=paste("motif",seq(1:dim(motifs.FFL)[1]), sep="")
      motifs.FFL=cbind(motifs.ids,motifs.FFL)
      write.table(motifs.FFL,file=file.path(output.path,"motifs.txt"),quote=F,row.names=F,col.names=T,sep="\t")
    }
    else{
      #write.table("No motifs found",file=file.path(output.path,"motifs.txt"),quote=F,row.names=F,col.names=T,sep="\t")
      motifs.FFL = c()
    }
      write("finished motifs 4 nodes", file=file.path(output.path,"finished4Nodemotifs.txt"),append=F,sep="\n")  
  }
}
PlotFunctionalSimilarity=function(genes,output.path)
{
###### @christian : these commented lines for testing only. u can try them urself
#   genes="ESR1, TP53, GRIN2D, AGER, AKT1, TERT, NCOA2, BBC3"
#   genes="CREB1, LTC4S, TLR9, IL5RA, MCAM, RPL10, RPS3A, ME2, CXCR4, SLC6A4, ERF, ID1, FLII, TGFB1, FLI1, UBE2I, PPRC1, CDC37, LRRFIP1, TGIF1, JAG1, TP53BP2, MSH6, MSH2"
#   genes="CREB1, RPL10, CXCR4, ID1, TGFB1, UBE2I, LRRFIP1, TGIF1, JAG1, TP53BP2, MSH6, MSH2"
#   genes="ESR1, TP53, GRIN2D, AGER, AKT1, TERT, NCOA2, BBC3"
#   genes="SPI1, BACH1, GNA13, SACM1L, FLI1, RAB23, POLE4, MSH2, SERTAD2, SKI, PHC2, ATP6V1C1, MSH6, DHX40, DPP7, RCN2, CHAF1A, PKN2, MECP2, ARL5B, MYO1E, B2M, TYROBP, FLII, MSR1, P2RY10, WAS"
#   genes="SPI1, BACH1, GNA13, SACM1L, FLI1, RAB23, POLE4, MSH2, SERTAD2, SKI, PHC2, ATP6V1C1, MSH6, DHX40, DPP7, RCN2, CHAF1A, PKN2, MECP2, ARL5B, MYO1E, B2M, TYROBP, FLII, MSR1, P2RY10, WAS"
#   output.path="output/disease/funsim.png"
#   ############################################
  
  print(output.path)

  genes=as.vector(unlist(strsplit(genes,",")))
  dput(genes)
  genes.entrez=unique(as.vector(unlist(mget(as.character(genes), envir=org.Hs.egALIAS2EG, ifnotfound=NA))))
  gosem=mgeneSim(genes.entrez,organism="human",measure="Wang")#,ont="BP"
  gosem=gosem[upper.tri(gosem)]
  
  all.entrez.genes <- mappedkeys(org.Hs.egACCNUM)
  #pvals.ks=c()
  #pvals.t=c()
  #pvals.wc=c()
  gosem.random.vector=c()
  for(i in 1: as.integer(config$NO_OF_random_permutations_for_functional_similarity))
  {
    genes.random=sample(all.entrez.genes,length(genes.entrez),replace = FALSE)
    gosem.random=mgeneSim(genes.random,organism="human",measure="Wang")#,ont="BP"
    gosem.random=gosem.random[upper.tri(gosem.random)]
    gosem.random.vector=c(gosem.random.vector,gosem.random)
    #     if(length(gosem.random)>1)
    #     {
    #       pvals.ks=c(pvals.ks,ks.test(gosem,gosem.random,alternative="l")$p.value)
    #       pvals.wc=c(pvals.wc,wilcox.test(gosem,gosem.random,alternative="g")$p.value)
    #       pvals.t=c(pvals.t,t.test(gosem,gosem.random,alternative="g")$p.value)      
    #     }
  }  
  #   pval.t.final= (length(pvals.t[pvals.t > 0.05]) /  length(pvals.t))
  #   pval.ks.final= (length(pvals.ks[pvals.ks > 0.05]) /  length(pvals.ks))
  #   pval.wc.final= (length(pvals.wc[pvals.wc > 0.05]) /  length(pvals.wc))
  #   pval=min(median(pvals.t),median(pvals.wc),median(pvals.ks))
  gosem.random.forplot=sample(gosem.random.vector,length(gosem))
  pval=ks.test(gosem,gosem.random.forplot,alternative="l")$p.value
  
  
  CairoPNG(bg="transparent",output.path,width=as.integer(config$funsimilarity.diagram.width),height=as.integer(config$funsimilarity.diagram.height))
  plot(ecdf(gosem),col="red", xlim=range(c(gosem, gosem.random.forplot)) , main="",xlab="Pair-wise similarity score", ylab="Cumulative distribution")
  #lines(ecdf(gosem),col="red",type="l")
  grid()
  lines(ecdf(gosem.random.forplot))
  #text(0.9,0.05, col="blue", paste("P-value < ",round(pval,4) ,sep=""),cex=1, adj = c(0.5, 0.5))
  text(0.9,0.05, col="blue", paste("P-value < ",format(pval, scientific = TRUE,digits=2) ,sep=""),cex=0.8, adj = c(0.5, 0.5))
  #mtext(paste("P-value < ",round(pval,3) ,sep=""), adj = 1,col="blue")
  legend(bty="n","topleft",c("Motif genes CDF","Random genes CDF") ,pch=c(19,19), col=c("red","black") ,cex=1)
  dev.off()
  
}
