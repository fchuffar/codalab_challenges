# Authors: Raphael Bacher (UGA), Florent Chuffart (INSERM)
# raphael.bacher@univ-grenoble-alpes.fr
#
#---------------------------------------------

# SCORING
scoring_function = function(truth, pred, prefix=NULL){
  nbna = sum(is.na(pred))
  propNA = nbna/length(truth)

  lm = lm(truth~pred)
  res = residuals(lm)
  mse = mean(res^2)
  sc = list(MSE=mse, propNA=nbna/length(truth))
  if (!is.null(prefix)) {
    names(sc) = paste(names(sc), prefix, sep="_")
  }
  return(sc)
}







### SESSION
### DO NOT CHANGE THIS PART
# define input/output/ref/res from command line args (should in principle never be changed)
args = commandArgs(TRUE)
CHALLENGER_SESSION = length(args) == 0
print(paste("CHALLENGER_SESSION:", CHALLENGER_SESSION))
SERVER_SESSION = !CHALLENGER_SESSION
print(paste("SERVER_SESSION:", SERVER_SESSION)) 

# Environment variables
if (CHALLENGER_SESSION) {
  input = output = ref = res = "./"
} else {
  if (!exists("input"))  input = trimws(args[1]) #get args from command line and remove white spaces
  if (!exists("output")) output = trimws(args[2])
  if (!exists("ref"))    ref = "/ref/"
  if (!exists("res"))    res = "/res/"  
}

# read ref data (if on the server or admin)
data_full_filename = paste0(input, ref, "/data_full.rds")
ADMIN_SESSION = file.exists(data_full_filename)
print(paste("ADMIN_SESSION:", ADMIN_SESSION)) 

  
  
  
  
  
  
  
#  EVALUATION
if (ADMIN_SESSION) {
  # Ground truth
  d_full = readRDS(data_full_filename)
  gs = colnames(d_full)[11:31]
  dg_full = as.matrix(d_full[,gs])

  # get NA  
  nb_na = 20000
  gr = expand.grid(1:nrow(dg_full), 1:length(gs))
  set.seed(1) 
  idx = sample(1:nrow(gr), nb_na)
  na_coord = data.frame(sample=rownames(dg_full)[gr[idx,1]], gene=gs[gr[idx,2]], stringsAsFactors=FALSE)

  # generate seed with the secret (dg_full) 
  # Score compute on 10000 prediction vector, 
  # always the same but unkown.
  set.seed(round(sum(dg_full[,ncol(dg_full)])))
  idx_r2 = sample(1:nrow(na_coord), ceiling(nb_na/2))
  na_coord = na_coord[idx_r2,]
  rownames(na_coord) = NULL

  # Load submited results from participant
  predg_full = readRDS(paste0(input, res, "results.rds"))


  names(gs)=gs
  proj = unique(d_full$project)
  names(proj) = proj  
  proj = rev(sort(sapply(proj, function(pr) {nrow(na_coord[d_full[na_coord$sample,]$project %in% pr,])})))
  proj = proj[proj>0]
  proj = names(proj)
  names(proj) = substr(proj, 6, 1000)  

  sub_na_coord = list(
    list(GLOBAL=na_coord),
    list(NORMAL=na_coord[d_full[na_coord$sample,]$tissue_status=="normal",]),
    lapply(gs, function(g) {na_coord[na_coord$gene==g,]}),
    lapply(proj, function(pr) {na_coord[d_full[na_coord$sample,]$project %in% pr,]})
  )

  length(sub_na_coord)
  sub_na_coord = unlist(sub_na_coord, recursive=FALSE)
  length(sub_na_coord)
  sapply(sub_na_coord, nrow)

  
  scores = sapply(names(sub_na_coord), function(n){
    print(n)
    pred = apply(sub_na_coord[[n]], 1, function(l) {
      predg_full[l[["sample"]], l[["gene"]]]
    })
    truth = apply(sub_na_coord[[n]], 1, function(l) {
      dg_full[l[["sample"]], l[["gene"]]]
    })
    scoring_function(truth, pred)
  })

  write_scores = function(scores, output_file) {
    for (grp in colnames(scores)) {
      for (ind in rownames(scores)) {
        key = paste(grp, ind, sep="_")
        cat(sprintf(paste(key, ":%f\n", sep=""), scores[ind,grp]), file=output_file, append=!(grp==colnames(scores)[1]&ind==rownames(scores)[1]))              
      }
    }
    cat(readLines(output_file), sep = "\n")     
  }

  write_scores(scores, output_file=paste0(output,"/scores.txt"))
}







# EXPORT BUNDLE
# export codallab bundle (if local run with reference data, it means if admin run)
if (ADMIN_SESSION & CHALLENGER_SESSION) {
  print("Export codallab bundle, beause of (CHALLENGER_SESSION & ADMIN_SESSION).")

  # rmarkdown::render("overview.Rmd")
  # rmarkdown::render("evaluation.Rmd")
  # rmarkdown::render("submission_script.Rmd", output_file="data.html")


  write_board = function(scores) {
    output_file="board.yml"
        cat(paste("leaderboard:                                        \n" , sep=""), file=output_file, append=FALSE)              
        cat(paste("  columns:                                          \n" , sep=""), file=output_file, append=TRUE)              
    for (grp in colnames(scores)) {
      for (ind in rownames(scores)) {
        key = paste(grp, ind, sep="_")
        FIRST_ELEMENT = grp==colnames(scores)[1]&ind==rownames(scores)[1]
        cat(paste("    ", key, ":                                      \n" , sep=""), file=output_file, append=TRUE)              
        cat(paste("      label: ", key, "                              \n" , sep=""), file=output_file, append=TRUE)              
        cat(paste("      leaderboard: ", ifelse(FIRST_ELEMENT,"&","*"), "id001\n" , sep=""), file=output_file, append=TRUE)              
        if (FIRST_ELEMENT) {

        cat(paste("        label: Results                              \n" , sep=""), file=output_file, append=TRUE)              
        cat(paste("        rank: 1                                     \n" , sep=""), file=output_file, append=TRUE)              
        }        
        cat(paste("      rank: ",as.numeric(ind=="MSE"), "             \n" , sep=""), file=output_file, append=TRUE)              
        cat(paste("      sort: asc                                     \n" , sep=""), file=output_file, append=TRUE)              
      }
    }
        cat(paste("  leaderboards:           \n" , sep=""), file=output_file, append=TRUE)              
        cat(paste("    Results: *id001    \n" , sep=""), file=output_file, append=TRUE)              
        cat(readLines(output_file), sep = "\n")     
  }
  write_board(scores)

  # zip the bundle
  zip_filename = "reference_data.zip"
  zip(zip_filename, "data_full.rds")

  zip_filename = "scoring_program.zip"
  zip(zip_filename, "scoring_program/metadata")
  zip(zip_filename, "scoring_program/scoring.r")    

  zip_filename = "starting_kit.zip"
  zip(zip_filename, "data.rds")
  zip(zip_filename, "starting_kit.Rmd")
  zip(zip_filename, "overview.Rmd")
  zip(zip_filename, "evaluation.Rmd")
  zip(zip_filename, "submission_script.Rmd")
  zip(zip_filename, "scoring_program/metadata")
  zip(zip_filename, "scoring_program/scoring.r")
  zip(zip_filename, ".Rhistory")

  zip_filename = "./codalab_bundle.zip"
  zip(zip_filename, "competition.yaml")
  zip(zip_filename, "data.html")
  zip(zip_filename, "terms.html")
  zip(zip_filename, "evaluation.html")
  zip(zip_filename, "overview.html")
  zip(zip_filename, "get_starting_kit.html")
  zip(zip_filename, "logo.png")
  zip(zip_filename, "reference_data.zip")
  zip(zip_filename, "scoring_program.zip")
  zip(zip_filename, "starting_kit.zip")
  
  file.remove("reference_data.zip")
  file.remove("scoring_program.zip")
  file.remove("starting_kit.zip")    
}

# stop("EFN")

