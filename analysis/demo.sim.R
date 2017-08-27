# TO DO: Explain what this script does, and how to use it.

library(varbvs)

# SCRIPT PARAMETERS
# -----------------
n  <- 800   # Number of samples.
p  <- 2000  # Number of variables (genetic markers).
na <- 20    # Number of quantitative trait loci (QTLs).
se <- 4     # Variance of residual.
r  <- 0.5   # Proportion of variance in trait explained by QTLs.

# Candidate values for the prior log-odds of inclusion.
logodds <- seq(-3,-1,0.1)

# Set the random number generator seed.
set.seed(1)

# GENERATE DATA SET
# -----------------
# Generate the minor allele frequencies so that they are uniform over range
# [0.05,0.5]. Then simulate genotypes assuming all markers are uncorrelated
# (i.e., unlinked), according to the specified minor allele frequencies.
cat("1. GENERATING DATA SET.\n")
maf <- 0.05 + 0.45*runif(p)
X   <- (runif(n*p) < maf) +
       (runif(n*p) < maf)
X   <- matrix(as.double(X),n,p,byrow = TRUE)

# Generate additive effects for the markers so that exactly na of them have
# a nonzero effect on the trait.
i       <- sample(p,na)
beta    <- rep(0,p)
beta[i] <- rnorm(na)

# Generate random labels for the markers.
colnames(X) <- paste0("rs",sample(1e6,p))

# Adjust the QTL effects so that we control for the proportion of variance
# explained (r). That is, we adjust beta so that r = a/(a+1), where I've
# defined a = beta'*cov(X)*beta. Here, sb is the variance of the (nonzero)
# QTL effects.
sb   <- r/(1-r)/var(c(X %*% beta))
beta <- sqrt(sb*se) * beta

# Generate a random intercept.
mu <- rnorm(1)

# Generate the quantitative trait measurements.
y <- c(mu + X %*% beta + sqrt(se)*rnorm(n))

