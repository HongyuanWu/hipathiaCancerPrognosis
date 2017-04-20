# file:
#   runPredict.R
# description:
#   - runs double nested CV evaluation for prediction
#   - results will be saved to folder `Robj/`
#   - adapted to paralleling tasks on cluster run, see shell script`runPredict.sh` paired with parameters in `runPredict.param.txt` for an example
#   - requires many packages to run this script, check the headers of `../src/func.R`

# read cluster parameters from bash command line
ags <- commandArgs(trailingOnly = TRUE)
stopifnot(length(ags) == 10)
jobid <- as.integer(ags[1]) # job index or row index of param
xname <- as.character(ags[2]) # feature matrix or x
yname <- as.character(ags[3]) # label groups or y
prname <- as.character(ags[4]) # classifier or predictor
i.fold <- as.integer(ags[5]) # (outer) cv fold index ranging from 1 to (nfolds*nrepeats)
nfolds <- as.integer(ags[6]) # (outer) cv number of folds
nrepeats <- as.integer(ags[7]) # (outer) cv number of repeats
i.fold.inn <- as.integer(ags[8]) # inner cv fold index ranging from 0 to (nfolds.inn*nrepeats.inn)
nfolds.inn <- as.integer(ags[9]) # inner cv number of folds
nrepeats.inn <- as.integer(ags[10]) # inner cv number of repeats

datapath <- '../data/'

# start! ------------------------------------------------------------------

message("\nRunning job id : ", jobid, "\n")
message("\n\twith parameters ... \n", paste(ags[-1], collapse = "\t"), "\n")

if (i.fold.inn == 0) {
  # i.fold.inn equal to 0 indicates outer cv loop run
  objname <- paste0('ivres_', paste(ags[-1], collapse = '_'))
} else {
  # i.fold.inn not equal to 0 indicates corresponding inner cv loop run
  objname <- paste0('cvres_', paste(ags[-1], collapse = '_'))
}
objpath <- paste0('Robj/', objname, '.RData')

if (file.exists(objpath)) {
  message('job already done !!')
  quit(save = 'no')
}

message("loading functions ...")
source("../../src/func.R")

# load x
if (xname %in% c("eff.vals","path.vals","fun.vals","go.vals","path.genes.vals")) {
  xtr <- get(load(paste0(datapath,xname,".RData")))
  rm(list = xname)
} else if (xname == "other.genes.vals") {
  xtr <- cbind(
    get(load(paste0(datapath,"other.genes.vals.pt1.RData"))),
    get(load(paste0(datapath,"other.genes.vals.pt2.RData")))
  )
  rm(other.genes.vals.pt1, other.genes.vals.pt2)
} else if (xname == "genes.vals") {
  xtr <- cbind(
    get(load(paste0(datapath,"path.genes.vals.RData"))),
    get(load(paste0(datapath,"other.genes.vals.pt1.RData"))),
    get(load(paste0(datapath,"other.genes.vals.pt2.RData")))
  )
  rm(path.genes.vals, other.genes.vals.pt1, other.genes.vals.pt2)
} else if (xname == "eff.and.other.genes.vals") {
  xtr <- cbind(
    get(load(paste0(datapath,"eff.vals.RData"))),
    get(load(paste0(datapath,"other.genes.vals.pt1.RData"))),
    get(load(paste0(datapath,"other.genes.vals.pt2.RData")))
  )
  rm(eff.vals, other.genes.vals.pt1, other.genes.vals.pt2)
} else if (xname == "path.and.other.genes.vals") {
  xtr <- cbind(
    get(load(paste0(datapath,"path.vals.RData"))),
    get(load(paste0(datapath,"other.genes.vals.pt1.RData"))),
    get(load(paste0(datapath,"other.genes.vals.pt2.RData")))
  )
  rm(path.vals, other.genes.vals.pt1, other.genes.vals.pt2)
} else {
  stop("Indefinite xname specified!!")
}

# load y
ytr <- get(load(paste0(datapath,yname,".RData")))
rm(list = yname)

# keep only samples for which label info is available
samplelist <- intersect(rownames(xtr), names(ytr))
message('Sample size = ', length(samplelist))
xtr <- xtr[samplelist, ]
ytr <- ytr[samplelist]


message("creating data split for CV ...")
set.seed(94151402)
foldIndices <- caret::createMultiFolds(1:nrow(xtr), k = nfolds, times = nrepeats)
fold <- foldIndices[[i.fold]]
# create inner cv folds to select predictor
set.seed(19817919)
foldIndices.inn <- caret::createMultiFolds(fold, k = nfolds.inn, times = nrepeats.inn)
# set train-test split
if (i.fold.inn == 0) {
  train.fold <- fold
  test.fold <- setdiff(1:nrow(xtr), train.fold)
} else {
  train.fold <- fold[foldIndices.inn[[i.fold.inn]]]
  test.fold <- setdiff(fold, train.fold)
}


message('train and predict ... ')
res <- indepValidation(xtr = xtr[train.fold, , drop=F], ytr = ytr[train.fold], 
                       xtst = xtr[test.fold, , drop=F], ytst = ytr[test.fold], 
                       predictor = prname)
if (prname == "predictorRF") {
  message('and feature selection ... ')
  fsname <- gsub("^predictor", "featselect", prname)
  res$featlist.short <- names(get(fsname, mode = "function")(model = res$model))
}


assign(objname, res)
if (!dir.exists('Robj'))
  dir.create('Robj')
save(list = objname, file = objpath)
message('new job saved up !!')
quit(save = 'no')
