###########################################
#  - TFMir  Project                       #
#  - prepare DBs of gene-gene interacti   #
#  - 2016-2-25                            #
#  - Copyright: Maryam Nazarieh           #
###########################################
#cutoff_gg = 0.8
## ===============================
## read the Experimental databases
## ===============================
mentha = read.delim(config$mentha.file, header = FALSE)
mentha = mentha[which(mentha[3] > pp.cutoffp),]
mentha = mentha[,c(1,2)]
mentha["evidence"] = "Experimental"
mentha["source"] = "mentha"
names(mentha) = c("gene","gene","evidence","source")
mentha=mentha[!duplicated(mentha),]

## ===============================
## read the String databases
## ===============================
string = read.delim(config$string.file, header = FALSE)
string = string[which(string[3] > pp.cutoffp),]
string = string[,c(1,2)]
string["evidence"] = "Predicted"
string["source"] = "String"
names(string) = c("gene","gene","evidence","source")
string=string[!duplicated(string),]


gene.gene.db = rbind(mentha,string)
#gene.gene.db = mentha
#gene.gene.db = string
#gene.gene.db=gene.gene.db[! duplicated(mentha),]
gene.gene.db["category"]="gene-gene"
gene.gene.db=gene.gene.db[,c(5,1,2,3,4)]