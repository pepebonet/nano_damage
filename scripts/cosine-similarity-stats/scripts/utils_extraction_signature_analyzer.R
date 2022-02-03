
#-------------------
# I am copying here the copyright from the Broad as stipulated. This script actually has a mixture of their code and
# my code. The core of the Signature Extraction is entirely Broad production
#------------------


# Copyright (c) 2017, Broad Institute
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are
# met:
#     Redistributions of source code must retain the above copyright
#     notice, this list of conditions and the following disclaimer.
#     Redistributions in binary form must reproduce the above copyright
#     notice, this list of conditions and the following disclaimer in
#     the documentation and/or other materials provided with the
#     distribution.
#     Neither the name of the Broad Institute nor the names of its
#     contributors may be used to endorse or promote products derived
#     from this software without specific prior written permission.
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
# A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
# HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
# SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
# LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
# DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
# THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.



# load functions. We modified one function and also removed some of the imputed variables
source("SignatureAnalyzer.PCAWG.function.R")

# Newly implemented function to retrieve the best run of the extraction
BestRun <- function(outpath, name_sample){

  summary_file <- paste(outpath, name_sample,  paste(name_sample,  "summaries", "tsv", sep = "."), sep ='/')
  print(summary_file)
  summary_df <- read.table(summary_file, sep ='\t')
  best_run <-summary_df[with(summary_df, order(K, -posterior)),]$Run[[1]]
  best_run_path <- paste(outpath, name_sample, "temp",
                         paste(name_sample, best_run, "RData", sep = "."), sep ='/')

  return(best_run_path)

}

# Newly implemented function to retrieve the results
ExtractedResults <- function(best_run_path, word){
  # WE LOAD THE EXTRACTION FROM THE NON-HYPERMUTATORS. It will be hardcoded given the best solution for the nonhyper.
  load(file=best_run_path)

  W <- res[[1]]
  H <- res[[2]]
  fix_cols <- gsub("\\.", '-', colnames(H))
  colnames(H)<- fix_cols
  second_fix <- gsub("(.*)-(.*)", "\\1\\.\\2", colnames(H))
  colnames(H)<- second_fix

  index <- colSums(W) > 1

  W <- W[,index]
  H <- H[index,]
  W1536.PRIMARY <- W
  H1536.PRIMARY <- H

  for (i in 1:ncol(W)) {
    H1536.PRIMARY[i,] <- H1536.PRIMARY[i,]*colSums(W)[i]
    W1536.PRIMARY[,i] <- W1536.PRIMARY[,i]*rowSums(H)[i]
  }

  colnames(W1536.PRIMARY) <- paste(word,seq(1:ncol(W1536.PRIMARY)),sep="")
  rownames(H1536.PRIMARY) <- colnames(W1536.PRIMARY)
  W1536.PRIMARY.norm <- apply(W1536.PRIMARY,2,function(x) x/sum(x))

  return(list(W=W1536.PRIMARY, H=H1536.PRIMARY, W_norm = W1536.PRIMARY.norm  ))
}

# Newly implemented function to retrieve cosinus similarity with PCAWG signatures
cosinus.correlation <- function(W1, mutation_file) {

  # this should be in the same order as in the mutation input. Extended signatures have different order than 96 channels
  if (grepl("full", mutation_file, fixed=TRUE)){
    W2 <- read.table('data/signatures/PCAWG/sigProfiler_SBS_signatures_2018_03_28.extended_indexed.csv',
                     sep ='\t', header = TRUE)
  }
  else{
    # 96 channels sorted in a different manner
    W2 <- read.table('data/signatures/PCAWG/sigProfiler_SBS_signatures_2018_03_28.indexed.csv',
                     sep ='\t', header = TRUE)
  }
  K1 <- ncol(W1)
  K2 <- ncol(W2)
  corr <- array(0,dim=c(K1,K2))
  for (i in 1:K1) {
    for (j in 1:K2) {
      if (sum(W1[,i])!=0 & sum(W2[,j])!=0) {
        corr[i,j] <- W1[,i]%*%W2[,j]/sqrt(sum(W1[,i]^2))/sqrt(sum(W2[,j]^2))
      }
    }
  }
  rownames(corr) <- colnames(W1)
  colnames(corr) <- colnames(W2)

  corr.max <- apply(corr,1,function(x) max(x))
  corr.id <- apply(corr,1,function(x) which.max(x))
  df.summary <- data.frame(colnames(W1), colnames(W2)[corr.id], corr.id, corr.max)
  colnames(df.summary)<- c("Discovered", "PCAWG", "PCAWG_ID", "MAX")
  return(df.summary)
}

# transform the 1536 channels to 96 (slightly modified from BROAD script)
get.lego96.from.lego1536 <- function(W0) {

  W <- W0[1:1536,]

  # sorted 1536 index
  index1536 <- read.table('data/signatures/PCAWG/lego1536.index.tsv',
                          header=F,sep='\t',as.is=T,skip=0, stringsAsFactors = FALSE)[,1]

  # sorted 96 context
  context96 <- read.delim('data/signatures/PCAWG/lego96.index.tsv',
                          header=F,sep='\t',as.is=T,skip=0)

  context96.label <- context96[,2]

  # 96 channels are triplet_Alt, 1536 are pentamer_alt
  contig.1536 <- paste(substring(index1536,2,4),substring(index1536,6,6),sep="")
  contig.96 <- contig.1536
  for (i in 1:96) {
    contig.96[contig.1536%in%context96.label[i]] <- i
  }
  mapping <- data.frame(index1536,contig.1536,contig.96)
  W96 <- array(0,dim=c(96,ncol(W)))
  for (i in 1:96) {
    W96[i,] <- colSums(W[mapping$contig.96==i,])
  }
  rownames(W96) <- context96.label
  colnames(W96) <- colnames(W0)
  norm_W = t(t(W96)/colSums(W96))
  return(norm_W)
}

#==
# Wrapper for SignatureAnalyzer. It will do at the same time the snvs, indels and dbs. This will work for both
# ttype-specific extraction and the Pan not hypermutant extraction. It contains code from the original script
#==
SignatureAnalyzerExtraction <- function(mutation_file, mutdf, outpath, n.run, Kcol, tol) {

  # get the outpath
  name_sample <- as.character(strsplit(basename(mutation_file), split = '.dlm', fixed = TRUE)[[1]][1])

  # create outpaths and define output
  dir.create(file.path(outpath), showWarnings = FALSE,  recursive = TRUE)
  dir.create(file.path(outpath, name_sample), showWarnings = FALSE)
  final_output <- paste(outpath, name_sample, "/", sep ="")
  final_output_temp <- paste(outpath, name_sample, "/", "temp/", sep ="")
  dir.create(file.path(final_output, "temp"), showWarnings = FALSE)

  OUTPUT <-final_output_temp

  # do iterations
  for (i in 1:n.run){

    # add seed so that the results are reproducible
    set.seed(12345+i)

    # first run of extraction
    res <- BayesNMF.L1W.L2H(as.matrix(mutdf), 2000000, 10, 5, tol, Kcol, Kcol, 1, final_output_temp)

    # save R object
    save(res, file = paste(OUTPUT, paste(name_sample,  i, "RData", sep = "."), sep = ""))

    W <- res[[1]]  # signature loading
    H <- res[[2]]  # activity loading
    index.W <- colSums(W)>1 ### only keep columns in W with non-zero contributions
    W <- W[,index.W]
    H <- H[index.W,]
    colsum <- colSums(W)
    rowsum <- rowSums(H)

    for (j in 1:ncol(W)) {
      W[,j] <- W[,j]*rowsum[j]
      H[j,] <- H[j,]*colsum[j]
    }

    write.table(W, file=paste(OUTPUT, paste(name_sample,  i, "processes.tsv", sep = "."), sep = ""),
                quote=FALSE, sep='\t', row.names =  FALSE)
    write.table(H, file=paste(OUTPUT, paste(name_sample,  i, "exposures.tsv", sep = "."), sep = ""),
                quote=FALSE, sep='\t', row.names =  FALSE)
  }

  # Create summary file
  summary.run <- data.frame()

  for (i in 1 : n.run) {
    load(file = paste(OUTPUT, paste(name_sample,  i, "RData", sep = "."), sep = ""))
    W <- res[[1]]
    index.W <- colSums(W) > 1
    W <- W[, index.W]
    K <- ncol(W)
    if (is.integer(K))
    {
      summary.run[i, 1] <- i
      summary.run[i, 2] <- K
      summary.run[i, 3] <- - res[[4]]
    }
  }

  summary.run <- summary.run[complete.cases(summary.run), ]

  colnames(summary.run) <- c("Run", "K", "posterior")
  summary.run <- data.frame(summary.run)
  summary.run <- summary.run[order(summary.run$posterior, decreasing = T),]

  # write summaries about runs
  write.table(summary.run, sep = "\t", file = paste(final_output, paste(name_sample, "summaries",
                                                                        "tsv", sep = "."), sep = ""))


  # extract the best results from the run and save them to files
  best_run <- BestRun(outpath, name_sample)
  run_signature_analyzer <- ExtractedResults(best_run, "W")

  W.norm <- run_signature_analyzer$W_norm
  W <- run_signature_analyzer$W
  H  <- run_signature_analyzer$H

  write.table(W.norm , file=paste(final_output, paste(name_sample, "processes.tsv", sep = "."), sep = ""),
              quote=FALSE, sep='\t', row.names =  FALSE)
  write.table(H, file=paste(final_output, paste(name_sample, "exposures.tsv", sep = "."), sep = ""),
              quote=FALSE, sep='\t', row.names =  FALSE)

  write.table(W, file=paste(final_output, paste(name_sample, "processes_full.tsv", sep = "."), sep = ""),
              quote=FALSE, sep='\t', row.names =  FALSE)

}

# if we are doing the full pipeline, this will extract also the hypermutants and do a final consensus merge
Hypermutators <- function(mutation_file, outpath, n.run, signatures, tol) {

  #===================
  # EXTRACTION OF HYPERMUTATORS
  #===================

  # type of mutation (snvs, dbs, indels)
  type_mutation <- as.character(strsplit(basename(mutation_file), split = '.', fixed = TRUE)[[1]][2])

  # hypermutant file path
  mutation_file_hypermutant <- gsub('NotHypermutated','Hypermutated',  mutation_file)

  # name of the cohort we are analyzing
  name_sample_not_hypermutant <- as.character(strsplit(basename(mutation_file), split = '.dlm', fixed = TRUE)[[1]][1])

  name_sample_hypermutant <- gsub('NotHypermutated', 'Hypermutated', name_sample_not_hypermutant)

  # get best run not hypermutant and load them
  best_run_not_hyper <- BestRun(outpath, name_sample_not_hypermutant)
  run_nothypermutant<- ExtractedResults(best_run_not_hyper, "Primary.W")

  print(paste("Doing ", name_sample_hypermutant, " with ", best_run_not_hyper, " as not hypermutant"))

  # read mutation file
  mutdf_hypermutants <- read.csv(mutation_file_hypermutant, sep ='\t')

  # merge DFs (not normalized)
  mutdf_hypermutants_nonhypermutants <- cbind(run_nothypermutant$W, mutdf_hypermutants)

  # signatures for DBS and INDELS
  signatures <- 40

  if (type_mutation == 'snvs'){
    signatures <- 96
  }

  SignatureAnalyzerExtraction(mutation_file_hypermutant, mutdf_hypermutants_nonhypermutants, outpath, n.run, signatures, tol)

  #==========================
  # MERGE AND GO SAMPLE BY SAMPLE
  #==========================

  name_sample_not_hypermutant <- as.character(strsplit(basename(mutation_file), split = '.dlm', fixed = TRUE)[[1]][1])
  name_sample_hypermutant <- gsub('NotHypermutated', 'Hypermutated', name_sample_not_hypermutant)

  full_mutation_file <-  gsub('NotHypermutated', '', mutation_file)
  name_sample_full <- as.character(strsplit(basename(full_mutation_file), split = '.dlm', fixed = TRUE)[[1]][1])

  lego1536.PAN<- read.csv(full_mutation_file, sep ="\t", check.names=FALSE)

    # create outpaths and define output
  dir.create(file.path(outpath), showWarnings = FALSE)
  dir.create(file.path(outpath, name_sample_full), showWarnings = FALSE)
  final_output <- paste(outpath, name_sample_full, "/", sep ="")
  final_output_temp <- paste(outpath, name_sample_full, "/", "temp/", sep ="")
  dir.create(file.path(final_output, "temp"), showWarnings = FALSE)

  OUTPUT <-final_output_temp

  best_run_hyper <- BestRun(outpath, name_sample_hypermutant)
  run_hypermutant <- ExtractedResults(best_run_hyper, "Secondary.W")

  W1536.PRIMARY.norm <- run_nothypermutant$W_norm
  W1536.PRIMARY <- run_nothypermutant$W
  H1536.PRIMARY  <- run_nothypermutant$H

  W1536.SECONDARY.norm <- run_hypermutant$W_norm
  W1536.SECONDARY <- run_hypermutant$W
  H1536.SECONDARY  <- run_hypermutant$H

  # cosine similairy between W1536.SECONDARY (signatures in the second step)  and W1536.PRIMARY (signatures in the first step)
  corr <- plot.W.correlation(W1536.SECONDARY.norm,W1536.PRIMARY.norm)

  # if it is a pentamer, select only the 1536 channels
  if (grepl("full", mutation_file, fixed=TRUE)){
    corr <- plot.W.correlation(W1536.SECONDARY.norm[1:1536,],W1536.PRIMARY.norm[1:1536,])
  }

  corr.max <- apply(corr,1,function(x) max(x)) ### maximum cosine similarity of secondary signatures to primary signatures
  corr.id <- apply(corr,1,function(x) which.max(x))  ### primary signature ID with corr.max

  ###### summary data-frame mapping the signatures in the second extration step to the signatures in the frist step
  df.summary <- data.frame(colnames(W1536.SECONDARY),colnames(W1536.PRIMARY)[corr.id],corr.id,corr.max) ### summary for the comparison of secondary signatures to primary ones
  colnames(df.summary) <- c("secondary","primary","primary.id","CS.max")
  df.summary[,"CS_0.99"] <- df.summary$CS.max > 0.99

  ###### identify primary signatures retained in the secondary signature extraction (W1536.1) and those attributions in 2624 samples (H1536.1)
  W1536.1 <- W1536.PRIMARY[,match(df.summary$primary[df.summary$CS_0.99],colnames(W1536.PRIMARY),nomatch=0)]
  H1536.1 <- H1536.PRIMARY[match(df.summary$primary[df.summary$CS_0.99],rownames(H1536.PRIMARY),nomatch=0),]
  colnames(W1536.1) <- df.summary$secondary[df.summary$CS_0.99]
  rownames(H1536.1) <- colnames(W1536.1)

  ###### identify signatures unique to the hyper-mutated samples (W1536.2) in the secondary signature extraction and those attributions in 156 samples (H1536.2.2)
  ###### H1536.2.1 is a signature attribution of the primary signatures in hyper-mutated samples.
  W1536.2 <- W1536.SECONDARY[,!df.summary$CS_0.99]

  H1536.2.1 <- H1536.SECONDARY[df.summary$CS_0.99,(ncol(W1536.PRIMARY)+1):ncol(H1536.SECONDARY)]
  H1536.2.2 <- H1536.SECONDARY[!df.summary$CS_0.99,(ncol(W1536.PRIMARY)+1):ncol(H1536.SECONDARY)]


  n.hyper <- ncol(W1536.2) ### # of signatures unique to hyper-mutated samples

  W.tmp1 <- W1536.1
  W.tmp2 <- cbind(W1536.1,W1536.2)

  #=======
  # ATRIBUTION TO NON HYPERMUTATORS
  #=======

  W1536.1.norm <- apply(W1536.1,2,function(x) x/sum(x))
  W1536.2.norm <- apply(W1536.2,2,function(x) x/sum(x))

  W0 <- W1536.1.norm #### fixed PRIMARY W
  H0 <- H1536.1      #### initialization for PRIMARY H
  V0 <- lego1536.PAN[,match(colnames(H0),colnames(lego1536.PAN),nomatch=0)] #### COMPOSITE lego matrix for 2624 samples
  K0 <- ncol(W0)     #### # of signatures


  # get the ttype names
  x1 <- sapply(colnames(V0),function(x) strsplit(x,"\\.")[[1]][1])
  x2 <- sapply(colnames(V0),function(x) strsplit(x,"\\.")[[1]][2])

  ttype <- x1
  ttype.unique <- unique(ttype)
  n.ttype.unique <- length(ttype.unique)

  Z0 <- array(1,dim=c(K0,ncol(V0)))

  colnames(Z0) <- colnames(H0)
  rownames(Z0) <- rownames(H0)

  a0 <- 10 # default parameter
  phi <- 1.5 # default parameter

  #---------------
  # this is a new implementation. For SNVs only allow signature 4 to selected tumor types (as they mention in the
  # algorithm, although it was harcoded for an already known signature number
  #--------------
  if (type_mutation =='snvs'){

    W0_to_test<-W0
    # if it is a pentamer, collapse it to 96 channels
    if (grepl("full", mutation_file, fixed=TRUE)){
      W0_to_test <-  get.lego96.from.lego1536(W0)
    }

    # if it is a 192 channels, collapse it to 96
    if ( grepl("replication", mutation_file, fixed=TRUE) || grepl("transcription", mutation_file, fixed=TRUE) || grepl("RT", mutation_file, fixed=TRUE) ){

      W1536.1.copy<-data.frame(W1536.1)
      W1536.1.copy<-cbind(data.frame(Class=((1:nrow(W1536.1.copy))+1) %/% 2), W1536.1.copy)
      W0_to_test_full<- aggregate(.~Class,data=W1536.1.copy, sum)
      W0_to_test_full <- subset(W0_to_test_full, select = -c(Class))
      W0_to_test <- apply(W0_to_test_full,2,function(x) x/sum(x))

    }
    cosinus_out <- cosinus.correlation(W0_to_test, mutation_file)
    similar_to_smoking <- cosinus_out[ (cosinus_out$MAX >0.9)&(cosinus_out$PCAWG =='SBS4'), ]

    for (i in similar_to_smoking$Discovered){
      similar_sig <- grep(paste(i, "$", sep =""), colnames(W0))
      Z0[similar_sig,!ttype%in%c("Head-and-neck","Lung","Unknown", "nan")] <- 0
    }
  }

  #=========================================
  # DO THE FITTING FOR THE NON-HYPERMUTATORS
  #=========================================
  for (i in 1:n.ttype.unique) {

    #### First step: determine a set of optimal signatures best explaining the observed mutations in each cohort.
    #### The automatic relevance determination technique was applied to the H matrix only, while keeping signatures (W0) frozen.
    cohort <- ttype.unique[i]
    print(cohort)
    W1 <- W0
    V1 <- as.matrix(V0[,ttype==cohort])
    Z1 <- Z0[,match(colnames(V1),colnames(Z0),nomatch=0)]
    H1 <- Z1*H0[,match(colnames(V1),colnames(Z0),nomatch=0)]

    # do it only if you have more than one sample
    if (length(colnames(H1))>1){
      lambda <- rep(1+(sqrt((a0-1)*(a0-2)*mean(V1,na.rm=T)/K0))/(nrow(V1) + ncol(V1) + a0 + 1)*2,K0)
      res0 <- BayesNMF.L1.KL.fixed_W.Z(as.matrix(V1),as.matrix(W1),as.matrix(H1),as.matrix(Z1),lambda,2000000,a0,1.e-07,1/phi)
      H2 <- res0[[2]]
      colnames(H2) <- colnames(V1)
      rownames(H2) <- colnames(W0)

      #### Second step: determine a sample-level signature attribution using selected sigantures in the first step.
      index.H2 <- rowSums(H2)>1 ### identify only active signatures in the cohort
      Z2 <- Z1
      Z2[!index.H2,] <- 0 ### only selected signatures in the first step are allowed + the original contraints on the signature availability from Z1.
      for (j in 1:ncol(H2)) {
        tmpH <- rep(0,ncol(W0))
        if (sum(V1[,j])>=5) {
          lambda <- 1+(sqrt((a0-1)*(a0-2)*mean(V1,na.rm=T)/K0))/(nrow(V1) + ncol(V1) + a0 + 1)*2
          res <- BayesNMF.L1.KL.fixed_W.Z.sample(as.matrix(V1[,j]),W0,as.matrix(H2[,j]),as.matrix(Z2[,j]),lambda,1000000,a0,1.e-07,1)
          tmpH <- res[[2]]
        }
        if (j==1) {
          H3 <- tmpH
        } else {
          H3 <- cbind(H3,tmpH)
          cat(j,'\n')
        }
      }
      colnames(H3) <- colnames(V1)
      rownames(H3) <- colnames(W0)
      if (i==1) {
        H2.all <- H2
        H3.all <- H3
      } else {
        H2.all <- cbind(H2.all,H2)
        H3.all <- cbind(H3.all,H3)
      }
    }
  }

  H.COMPOSITE.nohyper <- H3.all # attributions of COMPOSITE signatures in the PRIMARY sample set

  #=========================================
  # DO THE FITTING FOR THE HYPERMUTATORS
  #=========================================

  W0 <- cbind(W1536.1.norm, W1536.2.norm) ### fixed (PRIMARY + SECONDARY) W
  H0 <- rbind(H1536.2.1, H1536.2.2)       ### initialization for (PRIMARY + SECONDARY)

  V0 <- lego1536.PAN[,match(colnames(H0),colnames(lego1536.PAN),nomatch=0)] #### COMPOSITE lego matrix
  K0 <- ncol(W0)

  x0 <- colnames(H0)
  x1 <- sapply(colnames(V0),function(x) strsplit(x,"\\.")[[1]][1])
  ttype <- x1

  ttype.unique <- unique(ttype)
  n.ttype.unique <- length(ttype.unique)

  Z0 <- array(1,dim=c(K0,ncol(V0))) ### signature indicator matrix Z (0 = not allowed, 1 = allowed); all elements initialized by one.
  colnames(Z0) <- colnames(H0)
  rownames(Z0) <- rownames(H0)


  #---------------
  # this is a new implementation. For SNVs only allow signature 4 to selected tumor types (as they mention in the
  # algorithm, although it was harcoded for an already known signature number
  #--------------

  # for SNVs only allow signature 4 to selected tumor types
  if (type_mutation =='snvs'){

    W0_to_test<-W0
    # if it is a pentamer, collapse it to 96 channels
    if (grepl("full", mutation_file, fixed=TRUE)){
      W0_to_test <-  get.lego96.from.lego1536(W0)
    }

    # if it is a 192 channels, collapse it to 96
    if ( grepl("replication", mutation_file, fixed=TRUE) || grepl("transcription", mutation_file, fixed=TRUE) || grepl("RT", mutation_file, fixed=TRUE) ){

      W.copy<- cbind(W1536.1, W1536.2)
      W.copy<-cbind(data.frame(Class=((1:nrow(W.copy))+1) %/% 2),W.copy)
      W0_to_test_full<- aggregate(.~Class,data=W.copy, sum)
      W0_to_test_full <- subset(W0_to_test_full, select = -c(Class))
      W0_to_test <- apply(W0_to_test_full,2,function(x) x/sum(x))

    }
    cosinus_out <- cosinus.correlation(W0_to_test, mutation_file)
    similar_to_smoking <- cosinus_out[ (cosinus_out$MAX >0.925)&(cosinus_out$PCAWG =='SBS4'), ]

    for (i in similar_to_smoking$Discovered){
      similar_sig <- grep(paste(i, "$", sep =""), colnames(W0))
      Z0[similar_sig,!ttype%in%c("Head-and-neck","Lung","Unknown", "nan")] <- 0
    }
  }
  for (i in 1:n.ttype.unique) {
    #### First step: determine a set of optimal signatures best explaining the observed mutations in each cohort.
    #### The automatic relevance determination technique was applied to the H matrix only, while keeping signatures (W0) frozen.
    cohort <- ttype.unique[i]
    print(cohort)

    W1 <- W0
    V1 <- as.matrix(V0[,which(ttype==cohort)])

    Z1 <- as.matrix(Z0[,match(colnames(V1),colnames(Z0),nomatch=0)])
    H1 <- Z1*as.matrix(H0[,match(colnames(V1),colnames(Z0),nomatch=0)])
    if (length(colnames(H1))>1){

      lambda <- rep(1+(sqrt((a0-1)*(a0-2)*mean(V1,na.rm=T)/K0))/(nrow(V1) + ncol(V1) + a0 + 1)*2,K0)
      res0 <- BayesNMF.L1.KL.fixed_W.Z(as.matrix(V1),as.matrix(W1),as.matrix(H1),as.matrix(Z1),lambda,2000000,a0,1.e-07,1/phi)
      H2 <- res0[[2]]
      colnames(H2) <- colnames(V1)
      rownames(H2) <- colnames(W0)

      #### Second step: determine a sample-level signature attribution using selected sigantures in the first step.
      index.H2 <- rowSums(H2)>1 ### identify only active signatures in the cohort
      Z2 <- Z1
      Z2[!index.H2,] <- 0 ### only selected signatures in the first step are allowed + the original contraints on the signature availability from Z1.
      for (j in 1:ncol(H2)) {
        tmpH <- rep(0,ncol(W0))
        if (sum(V1[,j])>=5) {
          lambda <- 1+(sqrt((a0-1)*(a0-2)*mean(V1,na.rm=T)/K0))/(nrow(V1) + ncol(V1) + a0 + 1)*2
          res <- BayesNMF.L1.KL.fixed_W.Z.sample(as.matrix(V1[,j]),W0,as.matrix(H2[,j]),as.matrix(Z2[,j]),lambda,1000000,a0,1.e-07,1)
          tmpH <- res[[2]]
        }
        if (j==1) {
          H3 <- tmpH
        } else {
          H3 <- cbind(H3,tmpH)
          cat(j,'\n')
        }
      }
      colnames(H3) <- colnames(V1)
      rownames(H3) <- colnames(W0)
      if (i==1) {
        H2.all <- H2
        H3.all <- H3
      } else {
        H2.all <- cbind(H2.all,H2)
        H3.all <- cbind(H3.all,H3)
      }
    }
  }

  H.COMPOSITE.hyper <- H3.all # attributions of COMPOSITE signatures in the SECONDARY sample set

  tmp1 <- array(0,dim=c(n.hyper,ncol(H.COMPOSITE.nohyper)))
  rownames(tmp1) <- colnames(W1536.2)
  colnames(tmp1) <- colnames(H.COMPOSITE.nohyper)

  H.COMPOSITE <- cbind(rbind(H.COMPOSITE.nohyper,tmp1),H.COMPOSITE.hyper) # attributions of COMPOSITE signatures
  W.COMPOSITE <- W0

  dir.create(file.path(outpath, name_sample_full), showWarnings = FALSE)

  write.table(H.COMPOSITE, sep = "\t", file = paste(outpath, "/", name_sample_full, "/",  name_sample_full, ".exposures.tsv",sep=""))
  write.table(W.COMPOSITE, sep = "\t", file = paste(outpath, "/", name_sample_full, "/",  name_sample_full, ".processes.tsv",sep=""))
  W0_full <- cbind(W1536.1, W1536.2)
  write.table(W0_full, sep = "\t", file = paste(outpath, "/", name_sample_full, "/",  name_sample_full, ".processes_full.tsv",sep=""))
}

