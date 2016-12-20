###########################################
#  - TFMir  Project                       #
#  - TFMir  Initialize and read files     #
#  - 2014-10-1                            #
#  - Copyright: Mohamed Hamed             #
###########################################

#R programming environments:
#- R studio IDE
#- R version 2.12.0 (2010-10-15)
#-Platform: x86_64-apple-darwin9.8.0/x86_64 (64-bit)
#-locale:   [1] C/en_US.UTF-8/C/C/C/C


## ========================================
## Load the databases of interactions (dbs)
## ========================================
if(identical(speciesp,"Human")){
#if(species=="Human"){
  env=org.Hs.egALIAS2EG
  symbol="hsa-"
  load("databases2_TFmiR_human.RDATA")
  #load("databases2.RDATA")
  #load("databases2_TFmiR_mentha.RDATA")
  #load("databases2_TFmiR_string.RDATA")
  #load("databases_human.RDATA")
}
#if(species=="Mouse"){
if(identical(speciesp,"Mouse")){
  env=org.Mm.egALIAS2EG
  symbol="mmu-"
  #load("databases_mouse.RDATA")
  #load("databases2_TFmiR_string_mouse.RDATA")
  #load("databases2_TFmiR_mentha_mouse.RDATA")
  load("databases2_TFmiR_mouse.RDATA")
}
#load("databases2.RDATA")

#dbs=dbs[! (dbs$category=="tf-gene" & dbs$evidence=="Predicted"),]   # remove all predicted links
##Maryam ==================================
## read the gene.process.association and miRNA.process.association
## ==================================
mirna.process.association=read.delim(config$mirna.process.association.file,header=T,sep="\t")
names(mirna.process.association)=c("mirna","mirna.process")
mirna.process.association=mirna.process.association[!duplicated(mirna.process.association),]

gene.process.association=read.delim(config$gene.process.association.file,header=T,sep="\t")
names(gene.process.association)=c("gene","gene.process")
gene.process.association=gene.process.association[!duplicated(gene.process.association),]
##Maryam ==================================
## read the gene.tissue.association
## ==================================
gene.tissue.association=read.delim(config$gene.tissue.network.association.file,header=T,sep="\t")
names(gene.tissue.association)=c("gene","gene.tissue")
gene.tissue.association=gene.tissue.association[!duplicated(gene.tissue.association),]
getTissue_specificGenes = function(tissue.selected)
{
  genes.tissue = list()
  gene.tissue.association = read.delim(config$gene.tissue.association.file)
  tissue.option <- unique(as.character(unlist(gene.tissue.association$specific.tissue)))
  genes.tissue = c()
  for(i in 1:length(tissue.option)){
    if(tissue.option[i]==tissue.selected){
      genes.tissue <- gene.tissue.association[gene.tissue.association$specific.tissue==tissue.option[i],][1]
    }
  }
  return(genes.tissue)
}
## read the gene.tissue..network.association
## =========================================
gene.tissue.network.association = read.delim(config$gene.tissue.network.association.file)
tissue.option <- unique(as.character(unlist(gene.tissue.network.association$specific.tissue)))
genes.tissue = list()

## ==================================
## read the gene.disease.association
## ==================================

gene.disease.association=read.delim(config$gene.disease.association.file)
names(gene.disease.association)=c("gene","gene.disease")
gene.disease.association=gene.disease.association[!duplicated(gene.disease.association),]

# gene.disease.association.genecards=read.delim(config$gene.disease.association.file.genecard)
# names(gene.disease.association.genecards)=c("gene","gene.disease")
# gene.disease.association.disgenet=read.delim(config$gene.disease.association.file.disgenet)
# gene.disease.association.disgenet=gene.disease.association.disgenet[,c("geneSymbol","diseaseName")]
# names(gene.disease.association.disgenet)=c("gene","gene.disease")
# gene.disease.association=rbind(gene.disease.association.genecards,gene.disease.association.disgenet)
# gene.disease.association=gene.disease.association[!duplicated(gene.disease.association),]
# gene.disease.association=gene.disease.association[order(gene.disease.association$gene),]

## ==================================
## read the mirna.disease.association
## ==================================

mirna.disease.association=read.delim(config$mirna.disease.association.file,header=F)
mirna.disease.association=mirna.disease.association[,c(2,3)]
names(mirna.disease.association)=c("mirna","mirna.disease")

## =====================================
## read the mirna.disease.spectrum.width
## =====================================

mirna.dsw.association=read.delim(config$mirna.dsw.file,header=F,skip=16)
names(mirna.dsw.association)=c("mirna","count","mirna.dsw")
mirna.dsw.association=mirna.dsw.association[! duplicated(mirna.dsw.association$mirna),]

## ===================================
## read the mirna.function.association
## ===================================

con <- file(config$mirna.function.association.file, "r", blocking = FALSE, )
res=readLines(con) # empty
close(con)
mirna.function.association=data.frame("mirna"=character(),"mirna.function"=character()) 
for(i in 11: length(res))
{
  s=res[i]
  strs=strsplit(s,"\t")
  df=data.frame(mirna= strs[[1]][3:length(strs[[1]])])
  df["mirna.function"]=strs[[1]][2]
  mirna.function.association=rbind(mirna.function.association,df)
}
mirna.function.association
names(mirna.function.association)

## ===================================================
## Export the combined list of gene and miRNA diseases
## ===================================================
dis.gene=unique(as.character(unlist(gene.disease.association$gene.disease)))
dis.mirna=unique(as.character(unlist(mirna.disease.association$mirna.disease)))
dis=unique(c(dis.gene,dis.mirna))
dis=dis[order(dis)]
write.table(dis,file="disease.txt",quote=F,row.names=F,col.name=F)


## =============================================================================
## design a function that return the disease or the function of a specific mirna
## =============================================================================

getMIRnaCategory=function(mirna,category)
{
  mirna.category=" "
  
  if(category=="function")
  {
    mirna.category= toString(unlist(mirna.function.association[mirna.function.association$mirna==mirna,]$mirna.function))
    
  }else if(category=="disease")
  {
    mirna.category= toString(unlist(mirna.disease.association[mirna.disease.association$mirna==mirna,]$mirna.disease))  
  }
  
  return (mirna.category)
}

getMIRna_DSW=function(mirna)
{
  mirna.dsw=0
  number=mirna.dsw.association[mirna.dsw.association$mirna==mirna,]$mirna.dsw
  if (length(number) >0){ 
    mirna.dsw=round(as.numeric((as.character(number))),digits=5)
  }  
  return(mirna.dsw)
}


getGeneDiseases=function(gene)
{
  diseases= toString(unlist(gene.disease.association[gene.disease.association$gene==gene,]$gene.disease))
  unique(diseases)
}


getmiRNAforDisease=function(disease)
{
  mirnas=mirna.disease.association[agrep(disease,mirna.disease.association$mirna.disease,ignore.case=T),]$mirna
  unique(as.character(mirnas))
}

getGenesforDisease=function(disease)
{
  genes=gene.disease.association[agrep(disease,gene.disease.association$gene.disease,ignore.case=T),]$gene
  unique(as.character(genes))
}

getEntrezIDs =function(genes)
{
  genes.str=unlist(strsplit(genes,","))
  entrezids=toString(as.vector(unlist(mget(as.character(genes.str), envir=env, ifnotfound=NA))))
  return(entrezids)
}

printGeneEntrezIDsMap =function(genes,output.path)
{
  entrezids=mget(as.character(genes), envir=env, ifnotfound=NA)
  genemap=data.frame(gene.symbol = rep(names(entrezids), lapply(entrezids, length)),gene.entrezid = unlist(entrezids))
  genemap=genemap[! duplicated(genemap$gene.symbol),]  
  write.table(genemap,file=file.path(output.path,"genemap.txt"),quote=F,row.names=F,col.names=T,sep="\t")
}
getGeneTissues=function(gene)#not exported yet
{
  tissues= toString(unlist(gene.tissue.association[gene.tissue.association$gene==gene,]$gene.tissue))
  unique(tissues)
}

getGenesforTissue=function(tissue)
{
  genes=gene.tissue.association[agrep(tissue,gene.tissue.association$gene.tissue,ignore.case=T),]$gene
  unique(as.character(genes))
}

getGenesforProcess=function(process)
{
  genes = as.character()
  for(p in process){
    genes=union(genes,gene.process.association[agrep(p,gene.process.association$gene.process,ignore.case=T),]$gene)
  }
    unique(as.character(genes))
}

getmiRNAforProcess=function(process)
{
  mirnas = c()
  for(p in process){
    mirnas = union(mirnas,mirna.process.association[agrep(p,mirna.process.association$mirna.process,ignore.case=T),]$mirna)
  }
    unique(as.character(mirnas))
}
getInteractions=function(category,reg.input,target.input,process,tissue,disease,evidence,output.path,pval.cutoff)
{
  
  db=dbs.all[dbs.all$category==category & dbs.all$evidence %in% evidence, ]
  db.res=db[db$regulator %in% reg.input & db$target %in% target.input,]
  reg.cat=strsplit(category,"-")[[1]][1]
  target.cat=strsplit(category,"-")[[1]][2]  
  
  no.of.interactions.associated.with.input.tissue=0
  no.of.interactions.associated.with.input.process=0
  if( dim(db.res)[1] > 0) 
  {  
   
    db.res["is_regulator_in_tissue"]="FALSE"
   db.res["is_target_in_tissue"]="FALSE"
   if(! is.na(tissue) & tissue !="")
    {
    genes.tissue=getGenesforTissue(tissue)
    if(dim(db.res[db.res$regulator %in% genes.tissue,])[1] > 0)
      db.res[db.res$regulator %in% genes.tissue,]$is_regulator_in_tissue="TRUE"
    if(dim(db.res[db.res$target %in% genes.tissue,])[1] > 0)
      db.res[db.res$target %in% genes.tissue,]$is_target_in_tissue="TRUE"
    no.of.interactions.associated.with.input.tissue=dim(db.res[db.res$is_regulator_in_tissue==TRUE | db.res$is_target_in_tissue==TRUE , ])[1]
    #print(no.of.interactions.associated.with.input.tissue)  
   }
   
     db.res["is_regulator_in_process"]="FALSE"
     db.res["is_target_in_process"]="FALSE"
     #if(! is.na(process) & process !="")
     #{
     if(length(process > 0)){
       mirna.process=getmiRNAforProcess(process)
       genes.process=getGenesforProcess(process)
       print(typeof(mirna.process))
       print(length(mirna.process))
       print(length(genes.process))
       if(dim(db.res[db.res$regulator %in% mirna.process | db.res$regulator %in% genes.process,])[1] > 0)
         db.res[db.res$regulator %in% mirna.process | db.res$regulator %in% genes.process,]$is_regulator_in_process="TRUE"
       if(dim(db.res[db.res$target %in% mirna.process | db.res$target %in% genes.process,])[1] > 0)
         db.res[db.res$target %in% mirna.process | db.res$target %in% genes.process,]$is_target_in_process="TRUE"
       no.of.interactions.associated.with.input.process=dim(db.res[db.res$is_regulator_in_process==TRUE | db.res$is_target_in_process==TRUE , ])[1]
      print(no.of.interactions.associated.with.input.process)
     }
    output.path=file.path(output.path,category)
    if  (! file.exists(output.path)){
      dir.create(file.path(output.path))
    }  
   ## ===========================
   ## get attributes of regulators
   ## ===========================
    if(reg.cat %in% c("tf","gene"))
    {
      for(i in 1: dim(db.res)[1])
        {
          reg=as.character(db.res[i,]$regulator)
          db.res[i,"regulator.association.disease"]=getGeneDiseases(gene = reg)
          #reg.entrez=as.vector(unlist(mget(as.character(reg), envir=org.Hs.egALIAS2EG, ifnotfound=NA)))
          #if(length(reg.entrez) >1) reg.entrez=reg.entrez[1]
          #db.res[i,"regulator.david.report"]=paste("http://david.abcc.ncifcrf.gov/api.jsp?type=ENTREZ_GENE_ID&ids=", reg.entrez, "&tool=geneReportFull",sep="")
          db.res[i,"regulator.genecard"]=paste("http://www.genecards.org/cgi-bin/carddisp.pl?gene=", reg,sep="")    
        }
    }else{
        for(i in 1: dim(db.res)[1])   ## regulator is mirna
        {
          reg=as.character(db.res[i,]$regulator)
          db.res[i,"regulator.association.disease"]=getMIRnaCategory(mirna = reg,category = "disease")
          db.res[i,"regulator.association.function"]=getMIRnaCategory(mirna = reg,category = "function")
          db.res[i,"regulator.dsw"]=getMIRna_DSW(mirna = reg)
        }
      }
      
      ## ===========================
      ## get attributes of targets
      ## ===========================
      if(target.cat %in% c("tf","gene"))    {
        for(i in 1: dim(db.res)[1])
        {
          target=as.character(db.res[i,]$target)
          db.res[i,"target.association.disease"]=getGeneDiseases(gene = target)
          #db.res[i,"target.david.report"]=paste("http://david.abcc.ncifcrf.gov/api.jsp?type=ENTREZ_GENE_ID&ids=", as.vector(unlist(mget(as.character(target), envir=org.Hs.egALIAS2EG, ifnotfound=NA))), "&tool=geneReportFull",sep="")
          db.res[i,"target.genecard"]=paste("http://www.genecards.org/cgi-bin/carddisp.pl?gene=", target,sep="")    
        }
      }else{   ##### target is mirna
        for(i in 1: dim(db.res)[1])
        {
          target=as.character(db.res[i,]$target)
          db.res[i,"target.association.function"]=getMIRnaCategory(mirna = target,category = "function")
          db.res[i,"target.association.disease"]=getMIRnaCategory(mirna = target,category = "disease")
          db.res[i,"target.dsw"]=getMIRna_DSW(mirna = target)
        }
      }
      
      
      no.of.interactions.associated.with.input.disease=0
      db.res["is_regulator_in_disease"]="FALSE"
      db.res["is_target_in_disease"]="FALSE"
      if(! is.na(disease) & disease !="")
      {
        mirna.dis=getmiRNAforDisease(disease)
        genes.dis=getGenesforDisease(disease)
        if(dim(db.res[db.res$regulator %in% mirna.dis | db.res$regulator %in% genes.dis,])[1] > 0)
          db.res[db.res$regulator %in% mirna.dis | db.res$regulator %in% genes.dis,]$is_regulator_in_disease="TRUE"
        if(dim(db.res[db.res$target %in% mirna.dis | db.res$target %in% genes.dis,])[1] > 0)
          db.res[db.res$target %in% mirna.dis | db.res$target %in% genes.dis,]$is_target_in_disease="TRUE"
        
        no.of.interactions.associated.with.input.disease=dim(db.res[db.res$is_target_in_disease==TRUE | db.res$is_regulator_in_disease==TRUE , ])[1]
      }
      
      
      write.table(db.res,file=file.path(output.path,"res.txt"), quote=F,row.names=F,col.names=T,sep="\t")    
      
      ## ==============================================================================
      ## Venn Diagram and HYper geometric test / simualtion test and / permutation test
      ## ==============================================================================
      
      alldys.targets.by.reg.inputs=as.character(unique(db[db$regulator %in% reg.input, ]$target))
      overlap=intersect(tolower(target.input),tolower(alldys.targets.by.reg.inputs))
      #   nodes=unique(c(as.vector(dbs$regulator),as.vector(dbs$targets)))
      #   total.mirna=length(grep("hsa-",nodes))
      #   total.gene=length(nodes)-total.mirna
      total.target=1
      if(target.cat %in% c("tf","gene")){ 
        total.target=as.integer(config$Total.No.of.genes)   ### 30000 as by the gene database
      } else{
        total.target=as.integer(config$Total.No.of.miRNA)   ### 1881 as by the mirbase database
      }
      venndiagram.overlap.pval.simulation=overlapSignificance_By_simulation(total=total.target,numgA = length(alldys.targets.by.reg.inputs),numgB = length(target.input), overlap = length(overlap),no_of_simulation = as.integer(config$NO.of.simulation.for.overlap.significance))
      #venndiagram.overlap.pval.python=overlapSignificance_By_Python(total=total.mirna,numgA = length(alldys.mirnas.by.tfs.inputs),numgB = length(mirnas.input), overlap = length(overlap))
      venndiagram.overlap.pval.hypergeom=overlapSignificance_By_HyperGEOM(total=total.target,numgA = length(alldys.targets.by.reg.inputs),numgB = length(target.input),overlap = length(overlap))
      
      list=list(Regulators_Targets=alldys.targets.by.reg.inputs,Dereg_Targets =target.input)
      plotColorVenn(list,file.path(output.path,"venn.png"),width=as.integer(config$venndiagram.width),hight=as.integer(config$venndiagram.hight))
      Venn.image.path=file.path(output.path,"venn.png")
      #   plotColorVennBar(list,file.path(output.path,"vennbar.png"),width=as.integer(config$venndiagram.width),hight=as.integer(config$venndiagram.hight))
      #   Venn.bar.image.path=file.path(output.path,"vennbar.png")   
      
      db.res.nodes=unique(c(as.vector(db.res$regulator),as.vector(db.res$target)))
      db.res.mirna=db.res.nodes[grep(symbol,db.res.nodes)]
      #db.res.gene=db.res.nodes[- grep("hsa-",db.res.nodes)]
      db.res.gene=setdiff(db.res.nodes,db.res.mirna)
      
      if((target.cat =="mirna" | reg.cat =="mirna"  ) & length(db.res.mirna) > 0)
      {
        ## =========================================
        ## ORA analysis of the miRNAs in the results 
        ## =========================================
        mirna.res.ora.function=DO_ORA_FOR_MIRNA(db.res.mirna,category="function",pval.cutoff=pval.cutoff)
        mirna.res.ora.disease=DO_ORA_FOR_MIRNA(db.res.mirna,category="disease",pval.cutoff = pval.cutoff)
        mirna.res.ora=rbind(mirna.res.ora.function, mirna.res.ora.disease)
        if(dim(mirna.res.ora)[2] == 1)
          mirna.res.ora = data.frame("Category"="NA","Term"="NA","Count"="NA","Mir"="NA","Percentage"="NA","Pval"="NA","Pval.BH"="NA","Pval.Bonf"="NA")
        #mirna.res.ora=c()
        write.table(mirna.res.ora,file=file.path(output.path,"mirna.ora.txt"), quote=F,row.names=F,col.names=T,sep="\t")    
      }
      ## =========================================
      ## ORA analysis of the genes in the results 
      ## =========================================
      db.res.gene.etrezID=toString(as.vector(unlist(mget(as.character(db.res.gene), envir=env, ifnotfound=NA))))
      david.BP=paste("http://david.abcc.ncifcrf.gov/api.jsp?type=ENTREZ_GENE_ID&ids=",db.res.gene.etrezID,",&tool=chartReport&annot=GOTERM_BP_ALL",sep="")
      david.KEGG=paste("http://david.abcc.ncifcrf.gov/api.jsp?type=ENTREZ_GENE_ID&ids=",db.res.gene.etrezID,",&tool=chartReport&annot=KEGG_PATHWAY",sep="")
      david.OMIM=paste("http://david.abcc.ncifcrf.gov/api.jsp?type=ENTREZ_GENE_ID&ids=",db.res.gene.etrezID,",&tool=chartReport&annot=OMIM_DISEASE",sep="")    
      david.functional.clust=paste("http://david.abcc.ncifcrf.gov/api.jsp?type=ENTREZ_GENE_ID&ids=",db.res.gene.etrezID,",&tool=term2term&annot=GOTERM_BP_ALL",sep="")    
      
      ## ======================================
      ## print the results to the output folder
      ## ======================================
      summary.file=file.path (output.path,"summary.txt")
      
      write(paste("Category.of.Interactions=",category,sep=""),      file=summary.file,append=F,sep="\n")
      write(paste("Evidence=",toString(evidence),sep=""),      file=summary.file,append=T,sep="\n")
      write(paste("Target.input.set=",length(target.input),sep=""),  file=summary.file,append=T,sep="\n")
      write(paste("Regulator.input.set=",length(reg.input),sep=""),      file=summary.file,append=T,sep="\n")
      write(paste("ORA.pval.cutoff=",pval.cutoff,sep=""),      file=summary.file,append=T,sep="\n")
      write(paste("Total.Targets.inDb=",length(unique(unlist(db$target))),sep=""),      file=summary.file,append=T,sep="\n")
      write(paste("Total.Regulators.inDb=",length(unique(unlist(db$regulator))),sep=""),      file=summary.file,append=T,sep="\n")
      write(paste("Deregulated.targets.by.inputRegulator=",length(alldys.targets.by.reg.inputs),sep=""),      file=summary.file,append=T,sep="\n")
      write(paste("Deregulated.targets.of.inputRegulator.and.in.input.targetlist.overlap=",length(overlap),sep=""),      file=summary.file,append=T,sep="\n")
      write(paste("no.of.resultant.interactions=",dim(db.res)[1],sep=""),      file=summary.file,append=T,sep="\n")    
      write(paste("no.of.interactions.associated.with.input.disease=",no.of.interactions.associated.with.input.disease,sep=""),      file=summary.file,append=T,sep="\n")  
      #write(paste("total.mirna.inDb=",length(unique(transmir$mirna)),sep=""),      file=summary.file,append=T,sep="\n")
      write(paste("venn.pval.hypergeom=",venndiagram.overlap.pval.hypergeom,sep=""),      file=summary.file,append=T,sep="\n")
      write(paste("venn.pval.simulation=",venndiagram.overlap.pval.simulation,sep=""),      file=summary.file,append=T,sep="\n")
      write(paste("Venn.image.path=",Venn.image.path,sep=""),      file=summary.file,append=T,sep="\n")
      #  write(paste("Venn.bar.image.path=",Venn.bar.image.path,sep=""),      file=summary.file,append=T,sep="\n")
      
      if(target.cat %in% c("tf","gene") | reg.cat %in% c("tf","gene")  )
      {
        summary.file=file.path (output.path,"genes.ora.txt")
        write(db.res.gene.etrezID,file=summary.file,append=F,sep="\n")
        
        #write(paste("david.BP.link=",david.BP,sep=""),      file=summary.file,append=F,sep="\n")
        #write(paste("david.KEGG.link=",david.KEGG,sep=""),      file=summary.file,append=T,sep="\n")
        #write(paste("david.OMIM.link=",david.OMIM,sep=""),      file=summary.file,append=T,sep="\n")
        #write(paste("david.functional.clust.link=",david.functional.clust,sep=""),  file=summary.file,append=T,sep="\n")
      }
      return (db.res)
  }  

  rest.columns=data.frame(regulator.association.disease=character(),target.association.disease=character(),
                          is_regulator_in_disease=numeric(),is_target_in_disease=numeric(),is_regulator_in_tissue=numeric(),is_target_in_tissue=numeric()
                          ,is_regulator_in_process=numeric(),is_target_in_process=numeric())
  return (cbind(db.res,rest.columns ))                
}        
"
getTissueNetworkGene = function(tissue.selected)
{
  genes.tissue = c()
  for(i in 1:length(tissue.option)){
    if(tissue.option[i] == tissue.selected){
      genes.tissue <- gene.tissue.network.association[gene.tissue.network.association$specific.tissue==tissue.option[i],][2]
    }
  }
  return(genes.tissue)
}
#all.res.tissue.reg = intersect(all.res$regulator,tissue.specific.gene$regulator)
    #all.res.tissue.target = intersect(all.res$target,tissue.specific.gene$regulator)
    #all.res.tissue.specific = union(as.list(all.res.tissue.reg),as.list(all.res.tissue.target))
    
    #all.res.disease.tissue.reg = intersect(all.res.disease$regulator,tissue.specific.gene$regulator)
    #all.res.disease.tissue.target = intersect(all.res.disease$target,tissue.specific.gene$regulator)
    #all.res.disease.tissue.specific = union(as.list(all.res.disease.tissue.reg),as.list(all.res.disease.tissue.target))
    
    #all.res.normal.tissue.reg = intersect(all.res.tissue$regulator,tissue.specific.gene$regulator)
    #all.res.normal.tissue.target = intersect(all.res.tissue$target,tissue.specific.gene$regulator)
    #all.res.normal.tissue.specific = union(as.list(all.res.normal.tissue.reg),as.list(all.res.normal.tissue.target))
    
    #all.res.process.tissue.reg = intersect(all.res.process$regulator,tissue.specific.gene$regulator)
    #all.res.process.tissue.target = intersect(all.res.process$target,tissue.specific.gene$regulator)
    #all.res.process.tissue.specific = union(as.list(all.res.process.tissue.reg),as.list(all.res.process.tissue.target))
"
