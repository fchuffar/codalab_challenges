# Authors: Raphael Bacher (UGA), Florent Chuffart (INSERM)
# raphael.bacher@univ-grenoble-alpes.fr
#
#---------------------------------------------

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
rds_filename = paste0(input, ref, "/data_full.rds")
ADMIN_SESSION = file.exists(rds_filename)
print(paste("ADMIN_SESSION:", ADMIN_SESSION)) 

  
if (ADMIN_SESSION) {
  # Ground truth
  d_full = readRDS(rds_filename)
  gs = colnames(d_full)[11:31]
  d_full = as.matrix(d_full[,gs])

  # get NA  
  nb_na = 20000
  gr = expand.grid(1:nrow(d_full), 1:length(gs))
  set.seed(1) 
  idx = sample(1:nrow(gr), nb_na)
  na_coord = data.frame(sample=rownames(d_full)[gr[idx,1]], gene=gs[gr[idx,2]], stringsAsFatcors=FALSE)

  # generate seed with the secret (d_full) 
  set.seed(round(sum(d_full[,ncol(d_full)])))
  idx_r2 = sample(1:nrow(na_coord), ceiling(nb_na/2))
  na_coord = na_coord[idx_r2,]
  rownames(na_coord) = NULL

  # Load submited results from participant
  pred_full = readRDS(paste0(input, res, "results.rds"))

  # Define scoring function
  score = function(truth, pred, output_file, prefix){
    nbna = sum(is.na(pred))
    propNA = nbna/length(truth)
    lm = lm(truth~pred)
    res = residuals(lm)
    mse = mean(res^2)
    sc = list(MSE=mse, propNA=nbna/length(truth), nbNA=nbna)
    if (!missing(prefix)) {
      names(sc) = paste(names(sc), prefix, sep="_")
    }
      
    cat(sprintf("n: %f\n", length(truth)), file=output_file, append=FALSE)      
    foo = sapply(names(sc), function(key){
      cat(sprintf(paste(key, ":%f\n", sep=""), sc[[key]]), file=output_file, append=TRUE)      
    })
    cat(readLines(output_file), sep = "\n") 

    # cat(sprintf("MSE:   %f\n", mse  ), file=output_file, append=FALSE)
    # cat(sprintf("nbNA:  %f\n", nbna ), file=output_file, append=TRUE)
    # cat(readLines(output_file), sep = "\n")
    #
    #
    # cat(sprintf("MSE:   %f\n", mse  ), file=output_file, append=FALSE)
    # cat(sprintf("propNA:  %f\n", )), file=output_file, append=TRUE)
    # cat(sprintf("nbNA:  %f\n", nbna ), file=output_file, append=TRUE)
    return(sc)
  }
  
  # Evaluate prediction
  pred = apply(na_coord, 1, function(l) {
    pred_full[l[["sample"]], l[["gene"]]]
  })

  truth = apply(na_coord, 1, function(l) {
    d_full[l[["sample"]], l[["gene"]]]
  })
  foo = score(truth, pred, output_file=paste0(output,"/scores.txt"))
  
  # R2 R3 R4
  # un score global
  # Un score par gène
  # Nb NA
  # Un score sur normal
  # Un score sur tumoral
  # Un score sur homme
  # Un score sur femme
  # Un score par projet

  foo = sapply(gs, function(g){
    sub_na_coord = na_coord[na_coord$gene==g,] 
    pred = apply(sub_na_coord, 1, function(l) {
      pred_full[l[["sample"]], l[["gene"]]]
    })
    truth = apply(sub_na_coord, 1, function(l) {
      d_full[l[["sample"]], l[["gene"]]]
    })
    # score(truth, pred, output_file=paste0(output,"/scores.txt"), prefix=g)
  })
}


# export codallab bundle (if local run with reference data, it means if admin run)
if (ADMIN_SESSION & CHALLENGER_SESSION) {
  print("Export codallab bundle, beause of (CHALLENGER_SESSION & ADMIN_SESSION).")
  zip_filename = "reference_data.zip"
  zip(zip_filename, "data_full.rds")

  zip_filename = "scoring_program.zip"
  zip(zip_filename, "scoring_program/metadata")
  zip(zip_filename, "scoring_program/scoring.r")    

  zip_filename = "starting_kit.zip"
  zip(zip_filename, "starting_kit.Rmd")
  zip(zip_filename, "data.rds")
  zip(zip_filename, "scoring_program/metadata")
  zip(zip_filename, "scoring_program/scoring.r")
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

