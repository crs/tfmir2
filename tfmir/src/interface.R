
start <- function(tf, mirna, enrich_pvalue, pval, ppIcut, evidence, species, disease, process, tissue, output) {


#setwd('/Users/chris/Website/TF-site/backend');

zz <- file(file.path(output, "all.Rout"), open="wt")
sink(zz, type="message")

#setwd(normalizePath(file.path(getwd(),'../../../backend')));
setwd(normalizePath(file.path(getwd(),'../backend')));
#setwd("/Library/WebServer/Documents/TFmiR/backend");
	#print (c('parameters are ', tf, mirna))

warning(getwd());
print(version);
print(getwd());
print('Before normalizing');
print(mirna);
tf <- normalizePath(tf)
mirna <- normalizePath(mirna)
output <- normalizePath(output)
print('Source TF Mir');
print(source('TFMir.R'));
#setwd('../tf');
print(c('Used directories', tf,mirna));
print('Call TFMir');
if (!file.exists(mirna)) {
	mirna = "";
}
if (!file.exists(tf)) {
	tf = "";
}
print(species);
#tissue = "Breast - Mammary Tissue"
TFMir(tf, mirna, enrich_pvalue, pval, ppIcut, evidence, species, disease, process,tissue, output);
print('TFMIR finished');
}
