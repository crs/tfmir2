###########################################
#  - TFMir  Project                       #
#  - TFMir  Demo test for the TFMIR script#
#  - 2014-10-1                            #
#  - Copyright: Mohamed Hamed             #
###########################################
#R programming environments:
#- R studio IDE
#- R version 2.12.0 (2010-10-15)
#-Platform: x86_64-apple-darwin9.8.0/x86_64 (64-bit)
#-locale:   [1] C/en_US.UTF-8/C/C/C/C
setwd("/Library/WebServer/Documents/TFmiR/backend")
options(warn=1)
tf.path="./Data/genes.bc.diff.txt"
#mirna.path="mirna.sample.txt"
#tf.path=""
#tf.path="genes.bc.diff.txt"
mirna.path="./Data/mirna.bc.diff.txt"
#tf.path="/Users/Maryam/Desktop/prostate_mouse.txt"
#tf.path="Mouse_diabetes.txt"
#tf.path="genes.bc.diff_mouse.txt"
#mirna.path="mirna.bc.diff_mouse.txt"
#mirna.path=""
tissue = "Breast - Mammary Tissue"
#tissue="Testis"
#process="Cell differentiation"
processes = "cell death:Cell differentiation"
#process = "cell death"
#tissue=""
#mirna.path=""
#tissue = "Testis"
#tissue=Thyroid
#tf.path="genes.ad.diff.txt"
#mirna.path="mirna.ad.diff.txt"



#mirna.path=""
pval.cutoff=0.05
#disease="Melanoma"
#disease="Neoplasm"
#disease="halabessahbelzabadymolokhia"
disease="Breast Neoplasms"
#disease="Prostate Cancer"
#disease="Diabetes Mellitus"
#disease="Polycystic Liver Disease"
#disease="alzheimer"
#disease=""

output.path="output"
evidence="Experimental"
#evidence="Predicted"
#evidence="both"
enrich.pval = 0.5

species="Human"
#species="Mouse"
pp.cutoff = 0.8

source("TFMIR.R")

####generate networks
start.time <- Sys.time()
TFMir(tf.path,mirna.path,enrich.pval,pval.cutoff,pp.cutoff,evidence,species,disease,processes,tissue,output.path)
end.time <- Sys.time()
time.taken.TFMir <- end.time - start.time
time.taken.TFMir

####find 3-node and 4-node motifs
start.time <- Sys.time()
net.path="./output/disease/res.txt"
output.path="./output/disease"
ExportMotifs(net.path,output.path,evidence,species)
end.time <- Sys.time()
time.taken.TFMir <- end.time - start.time
time.taken.TFMir

####run 4-node motifs separately from 3-node motifs
start.time <- Sys.time()
net.path="./output/disease.4node/res.txt"
output.path="./output/disease.4node"
Export4NodeMotifs(net.path,output.path,evidence,species)
end.time <- Sys.time()
time.taken.TFMir <- end.time - start.time
time.taken.TFMir

