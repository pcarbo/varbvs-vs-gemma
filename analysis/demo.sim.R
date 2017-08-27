# TO DO: Explain what this script does, and how to use it.
library(ggplot2)
library(cowplot)
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
# Generate the minor allele frequencies so that they are uniform over
# range [0.05,0.5]. Then simulate genotypes assuming all markers are
# uncorrelated (i.e., unlinked), according to the minor allele
# frequencies specified by vector "maf".
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

# FIT VARIATIONAL APPROXIMATION TO POSTERIOR
# ------------------------------------------
# Fit the fully-factorized variational approximation to the posterior
# distribution of the coefficients for a linear regression model of a
# continuous outcome (quantitiative trait), with spike and slab priors on
# the coefficients. 
cat("2. FITTING VARBVS MODEL TO DATA.\n")
fit.varbvs <- varbvs(X,NULL,y,"gaussian",logodds = logodds,n0 = 0)
cat("\n")

cat("3. FITTING BSLMM MODEL TO DATA.\n")
# TO DO.

# TO DO: Add comments here explaining what this code chunk does.
cat("4. SUMMARIZING VARBVS & BSLMM RESULTS.\n")

# Compute the false positive rate (FPR) and true positive rate (TPR)
# at each PIP threshold.
dat <- data.frame(pip.cutoff = seq(0,1,0.01),
                  tpr.varbvs = 0,
                  fpr.varbvs = 0,
                  tpr.bslmm  = 0,
                  fpr.bslmm  = 0)
for (i in 1:101) {

  # Get the PIP threshold.
  r <- dat$pip.cutoff[i]

  # Compute the varbvs TPR and FPR.
  np                <- sum(fit.varbvs$pip >= r)
  ntp               <- sum(fit.varbvs$pip >= r & beta != 0)
  nfp               <- sum(fit.varbvs$pip >= r & beta == 0)
  dat$tpr.varbvs[i] <- ntp / np
  dat$fpr.varbvs[i] <- nfp / np
}

p1 <- ggplot(dat,aes(x = pip.cutoff,y = tpr.varbvs)) +
  geom_line()

p2 <- ggplot(dat,aes(x = pip.cutoff,y = fpr.varbvs)) +
  geom_line()

print(plot_grid(p1,p2))
