# submission_script.r 
# Authors: Magali Richard (CNRS), Florent Chuffart (INSERM)
#
#---------------------------------------------

### DO NOT CHANGE THIS PART
d = readRDS("data.rds")
gs = colnames(d)[11:31]

### PUT YOUR SCRIPT HERE
# Example 1
pred = as.matrix(d[,gs])
pred[] = sample(na.omit(pred), nrow(pred)*ncol(pred), replace=TRUE)

# # Example 2
# pred = sapply(gs, function(g, d){
#  m = lm(d[[g]] ~ d$tissue)
# #  m = lm(d[[g]] ~ d$age+d$sex+d$tissue+d$tissue_status+d$project)
#  predict(m, d)
# }, d)

### DO NOT CHANGE THIS PART
pred[!is.na(as.matrix(d[,gs]))] = NA # Sparse, only submit NA values!
saveRDS(pred, "results.rds")
zip_filename = paste(sep="",  "results_", format(Sys.time(), format="%m_%d_%Y_%s"), ".zip")
zip(zip_filename, "results.rds")
print(zip_filename)
