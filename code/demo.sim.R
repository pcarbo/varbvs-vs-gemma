# Empirical comparison of varbvs and GEMMA methods for genome-wide
# mapping of a quantitative trait in a simulated data set. The data
# are simulated assuming that all the genetic markers are uncorrelated
# with each other (i.e., they are "unlinked").

# SCRIPT SETTINGS
# ---------------
# Parameters controlling generation of data set.
n  <- 800   # Number of samples
p  <- 2000  # Number of variables or genetic markers.
na <- 20    # Number of quantitative trait loci (QTLs).
se <- 4     # Variance of the residual.
r  <- 0.5   # Proportion of variance in trait explained by QTLs.

# Length of Markov chain for simulating BVSR posterior (ns) and
# interval for saving state of Markov chain (m).
ns <- 1e5
m  <- 100

# Candidate values for the prior log-odds of inclusion in the varbvs
# model.
logodds <- seq(-3,-1,0.1)

# Initialize the pseudorandom number generator.
set.seed(1)

# Load the varbvs package, as well as some additional functions I
# implemented for this analysis.
library(varbvs)
source("../code/gemma.R")

# GENERATE DATA SET
# -----------------
# Generate the minor allele frequencies so that they are uniform over
# range [0.05,0.5]. Then simulate genotypes assuming all markers are
# uncorrelated (i.e., unlinked), according to the minor allele
# frequencies specified by vector maf.
cat("Generating data set.\n")
maf <- 0.05 + 0.45*runif(p)
X   <- (runif(n*p) < maf) +
       (runif(n*p) < maf)
X   <- matrix(as.double(X),n,p,byrow = TRUE)
colnames(X) <- paste0("rs",sample(1e6,p))

# Generate additive effects for the markers so that exactly "na" of them
# have a nonzero effect on the trait.
i       <- sample(p,na)
beta    <- rep(0,p)
beta[i] <- rnorm(na)

# Adjust the QTL effects so that we control for the proportion of
# variance explained (r). That is, we adjust beta so that r = a/(a+1),
# where I've defined a = beta'*cov(X)*beta. Here, sb is the variance
# of the (nonzero) QTL effects.
sb   <- r/(1-r)/var(c(X %*% beta))
beta <- sqrt(sb*se) * beta

# Generate a random intercept.
mu <- rnorm(1)

# Generate the quantitative trait measurements.
y <- c(mu + X %*% beta + sqrt(se)*rnorm(n))

# FIT VARBVS MODEL
# ----------------
# Fit the fully-factorized variational approximation to the posterior
# distribution of the coefficients for a linear regression model of a
# continuous outcome (quantitiative trait), with spike-and-slab priors on
# the coefficients.
cat("Fitting varbvs model.\n")
fit.varbvs <- varbvs(X,NULL,y,"gaussian",logodds = logodds,sa = sb,
                     verbose = FALSE)
					 
# FIT BVSR MODEL
# --------------
# Fit (nearly) the same model, but using an MCMC algorithm to simulate
# the posterior distribution instead of computing approximate
# posterior probabilities using a variational approximation. In the
# GEMMA software, this is called the BVSR model, short for "Bayesian
# variable selection in regression."
cat("Fitting BVSR model using GEMMA.\n")
fit.bvsr <- mcmc.bvsr(X,y,verbose = FALSE)

# SAVE RESULTS
# ------------
session.info <- sessionInfo()
save(list = c("ns","m","r","maf","mu","beta","X","y","se","fit.varbvs",
              "fit.bvsr","session.info"),file = "../output/demo.sim.out")


