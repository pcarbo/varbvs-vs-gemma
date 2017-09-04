# TO DO: Explain here what this function does.
bslmm <- function (X, y, ns = 1e5, m = 100, seed = 1) {

  # Get the number of predictors / markers
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

  # Fit the model using GEMMA.
  system(sprintf(paste("gemma -g geno.txt -p pheno.txt -bslmm 1 -rmin 1",
                       "-smax %d -w 0 -s %d -rpace %d -wpace 1 -seed %d",
                       "-mh 1"),
                 p,ns,m,seed))

  # Load the GEMMA results.
  out <- read.table("output/result.param.txt",sep = "\t",header = TRUE,
                    stringsAsFactors = FALSE)
  
  # Return the posterior inclusion probabilities.
  pip        <- out$gamma
  names(pip) <- ids
  return(list(pip = pip))
}
