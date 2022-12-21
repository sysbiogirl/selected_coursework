library(edgeR) 
infile     <- 'w08-data.13.txt' 
group      <- factor(c(1,1,1,2,2,2))
outfile    <- 'w08-data.13_out_norm.txt' 
x     <- read.table(infile, sep='	', row.names=1) 
y     <- DGEList(counts=x,group=group) 
y     <- calcNormFactors(y) 
y     <- estimateDisp(y) 
et    <- exactTest(y) 
tab   <- topTags(et, nrow(x)) 

write.table(tab, file=outfile) 
