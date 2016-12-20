cat="mirna-gene"

unique(dbs.all[dbs.all$category ==cat ,"source"])

df=dbs.all[dbs.all$category ==cat & dbs.all$source %in% c("Starbase") , ]

df=dbs.all
print( "## no of edges>>>"  )
dim(df)
print("==================")

nodes=unique( c(as.character(df$regulator), as.character(df$target)  ))
mirnas=nodes[grep("hsa-",nodes)]
print( " no of mirnas is >>>>>")
print(length(unique(mirnas)))
print("==================")

genes=setdiff(nodes,mirnas)
print( " no of genes is >>>>>")
length(unique(genes))



