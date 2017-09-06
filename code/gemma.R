# Fit the BVSR model (i.e., the BSLMM model with rmin=1) to the data
# using GEMMA, and output the marginal posterior inclusion
# probabilities (PIPs). Input argument ns specifies the length of the
# Markov chain for the MCMC algorithm, and input argument m specifies
# the interval for saving the state of the Markov chain. For example,
# if m = 100 and ns = 1e5, then the state will be saved 1,000 times.
mcmc.bvsr <- function (X, y, ns = 1e5, m = 100, seed = 1,
                       gemma.exec = "gemma", verbose = TRUE) {

  # Get the number of predictors/markers.
  p <- ncol(X)
    
  # Write the outcome (phenotype) data to file pheno.txt.
  write.table(data.frame(y),"pheno.txt",quote = FALSE,row.names = FALSE,
              col.names = FALSE)

  # Write the predictor (genotype) data to file geno.txt.
  ids <- colnames(X)
  X <- t(X)
  X <- as.data.frame(X)
  X <- round(X,digits = 6)
  X <- cbind(data.frame(id = ids,ref = "A",alt = "G"),X)
  write.table(X,"geno.txt",sep = " ",quote = FALSE,row.names = FALSE,
              col.names = FALSE)

  # Fit the BVSR model using GEMMA.
  system(sprintf(paste("%s -g geno.txt -p pheno.txt -bslmm 1 -rmin 1",
                       "-smax %d -w 0 -s %d -rpace %d -wpace 1 -seed %d",
                       "-mh 1"),
                 gemma.exec,p,ns,m,seed),
         ignore.stdout = !verbose)

  # Load the GEMMA results.
  out <- read.table("output/result.param.txt",sep = "\t",header = TRUE,
                    stringsAsFactors = FALSE)
  
  # Return the posterior inclusion probabilities (PIPs).
  pip        <- out$gamma
  names(pip) <- ids
  return(list(pip = pip))
}
