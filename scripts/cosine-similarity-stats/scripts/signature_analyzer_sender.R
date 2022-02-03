args = commandArgs(trailingOnly = TRUE)

# load functions. We modified one function and also removed some of the imputed variables
source("./utils_extraction_signature_analyzer.R")

#=======
# ARGS
#=======

mutation_file <- args[1]  # '/workspace/projects/mutfootprints/data/HMF/matrix_signatures/PanNotHypermutated_transcription.snvs.dlm
signatures <- strtoi(args[2])  # 72
outpath <- args[3]  # /home/opich/bg/temp/
pan <- args[4]  # PAN

# Number of BayesNMF runs
n.run <- 5

# tolerance based on whether we are analyzing Pan or not
tol <- 1.e-07

if (pan == "PAN") {
  tol <- 5.e-06
}

mutdf <- read.csv(mutation_file, sep = '\t')
SignatureAnalyzerExtraction(mutation_file, mutdf, outpath, n.run, signatures, tol)


# if it is a hypermutator, then run the nonhypermutator/hypermutator path
if(grepl("NotHypermutated", mutation_file, fixed=TRUE)){
    Hypermutators(mutation_file, outpath, n.run, signatures, tol)
}

