### DO NOT CHANGE THIS PART
d = readRDS(file="data.rds")
gs = c(
  "ATAD2", "SYCP3",                # Chromatin 
  "H19", "IGF2", "NNAT", "BLCAP",  # inprinted genes
  "BRD4", "BRDT", "NUTM1",         # Testis restricted
  "MAGEB6", "TUBA3C",              # Testis specific
  "SMYD3", "MAP3K2", "KDR",        # K methyl transferase
  "TP53", "KRAS", "BRAF",          # Tumor supressors
  "FASTKD1", "ND2", "ND3", "ND4"   # Mitochondiral activity
)

### PUT YOUR SCRIPT HERE

#Example 1
pred = as.matrix(d[,gs])
pred[] = sample(na.omit(pred), nrow(pred)*ncol(pred), replace=TRUE)
pred[!is.na(as.matrix(d[,gs]))] = NA # Sparse, only submit NA values!

#Example 2
#pred = sapply(gs, function(g, d){
#  m = lm(d[[g]] ~ d$tissue)
#  m = lm(d[[g]] ~ d$age+d$sex+d$tissue+d$tissue_status+d$project)
#  predict(m, d)
#}, d)
#pred[!is.na(as.matrix(d[,gs]))] = NA # Sparse, only submit NA values!


### DO NOT CHANGE THIS PART
saveRDS(pred, "results.rds")
zip_filename = paste(sep="",  "results_", format(Sys.time(), format="%m_%d_%Y_%s"), ".zip")
zip(zip_filename, "results.rds")
print(zip_filename)
