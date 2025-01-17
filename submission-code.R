### Simulations to replicate Section 3.3 of the supplementary materials.
### Alex Stringer
### 2024/04

# Instructions: where you see the word CHANGE in a comment, change the corresponding file path
# to an appropriate one on your computer.
# This script does not require any external data to run.
# The packages required are listed in the pkgs variable below.
# The only one not on CRAN is the aghqmm package, available here: https://github.com/awstringer1/aghqmm
# Code to install it is below.

## Set paths ##
# CHANGE the base path to whatever you want on your machine
basepath <- "~/work/projects/num-quadpoints/approx-glmm-theory/code"
# CHANGE the name of the simulation to control how saved results are named
simname <- "sims-20241213-v1"
stopifnot(dir.exists(basepath))
resultspath <- file.path(basepath,'results')
if (!dir.exists(resultspath)) dir.create(resultspath)
simresultsname <- paste0(simname,".RData")
simsprocessedname <- paste0(simname,".csv")

# CHANGE the path to where you downloaded the aghqmm package repo
# from https://github.com/awstringer1/aghqmm
aghqmmpath <- "~/work/projects/mixedmodel-computation/aghqmm" # CHANGE this

## After this point, nothing should need to be changed. ##

## Load Packages ##

install.packages(aghqmmpath,repos=NULL,type="source")

pkgs <- c(
  'tidyverse',
  'lme4',
  'Rcpp',
  'RcppEigen',
  'parallel'
)
for (pkg in pkgs) {
  if (!require(pkg,character.only = TRUE,quietly = TRUE)) {
    cat(paste0("Could not find package ",pkg,", installing from CRAN.\n"))
    install.packages(pkg)
    require(pkg,character.only = TRUE,quietly = TRUE)
  }
}


## Set Parameters ##
numruns <- 1  # Number of times to execute the simulations
numsims <- 1000 # Number of simulations in each category PER RUN
m <- c(1000, 2000, 4000, 8000, 16000)
p <- 1
n <- seq(2, 10, by = 2)
k <- seq(1, 11, by = 2)
beta <- c(-4,2)
S <- 4
# Create simulation objects
simstodoframe <- expand.grid(
  n = n,
  p = p,
  m = m,
  k = k,
  idx = 1
)
simlist <- rep(split(simstodoframe, seq(nrow(simstodoframe))), numsims)
for (i in 1:length(simlist)) simlist[[i]]$idx <- i

options(mc.cores = parallel::detectCores())
RNGkind("L'Ecuyer-CMRG") # For reproducibility with parallel

### Function to execute simulation ###

dosim <- function(lst) {
  # lst: simulation parameters
  m <-            lst$m             # Number of groups
  n <-            lst$n             # Number of observations per group
  beta <-         beta          # Regression parameters
  S <-            S
  k <-            lst$k             # Number of quadrature points

  cat("Sim: ",lst$idx," of ",length(simlist),"|m=",m,"|n=",n,"|k=",k,"...",sep="")
  
  # Simulate data
  # simdata <- aghqmm::simulate_data(m,n,beta,S)
  simdata <- aghqmm:::simulate_data_nonconstant(m,n,beta,S)
  # Fit model
  # opt <- tryCatch(aghqmm::aghqmm(y ~ x + (1|id),simdata,k=k,method = "both",control = bfgscontrol),error = function(e) e)
  tm <- Sys.time()
  opt <- tryCatch(lme4::glmer(y ~ x + (1|id), data = simdata, nAGQ = k, family = binomial))
  comptime <- as.numeric(difftime(Sys.time(), tm, units = "secs"))
  if (inherits(opt,'condition')) return(opt)
  
  optsum <- summary(opt)

  # Return results
  simresults <- data.frame(
    m=m,n=n,k=k,p=lst$p,
    totaltime = comptime
  )
  # get the true beta
  paramresults <- data.frame(
    true = beta,
    est = opt@optinfo$val[2:3],
    lowerWald = optsum$coefficients[ , 1] - 1.96 * optsum$coefficients[ , 2],
    upperWald = optsum$coefficients[ , 1] + 1.96 * optsum$coefficients[ , 2]
  ) %>%
    dplyr::mutate(
      covrWald = lowerWald <= true & upperWald >= true,
      lengthWald = upperWald - lowerWald
    )
  cat(" completed.\n",sep="")
  list(simresults=simresults,paramresults=paramresults)
}

processsimulation <- function(sim) {
  if (inherits(sim,'condition')) return(NULL)
  # Change this if the parameters change
  out <- sim$simresults

  out$beta0true <- sim$paramresults[1,'true']
  out$beta0est <- sim$paramresults[1,'est']
  out$beta0bias <- out$beta0est-out$beta0true
  out$beta0covrWald <- sim$paramresults[1,'covrWald']
  out$beta0lengthWald <- sim$paramresults[1,'lengthWald']

  out$beta1true <- sim$paramresults[2,'true']
  out$beta1est <- sim$paramresults[2,'est']
  out$beta1bias <- out$beta1est-out$beta1true
  out$beta1covrWald <- sim$paramresults[2,'covrWald']
  out$beta1lengthWald <- sim$paramresults[2,'lengthWald']

  out
}

### Do Simulations ###
set.seed(4936298)
mc.reset.stream() # Reproducbility in parallel
# Do the simulations
cat("Doing",length(simlist),"simulations...\n")
tm <- Sys.time()
# execute the simulations numruns times
simruns <- list()
length(simruns) <- numruns
for (b in 1:numruns) {
  simruns[[b]] <- mclapply(simlist,dosim)
}
sims <- Reduce(c,simruns)
simtime <- as.numeric(difftime(Sys.time(),tm,units='secs'))
cat("Finished simulations, they took",simtime,"seconds.\n")
cat("Saving simulations...\n")
save(sims,file=file.path(resultspath,simresultsname))
cat("Saved simulations to file:",file.path(resultspath,simresultsname),"\n")
cat("Processing simulations...\n")
simsprocessed <- as_tibble(dplyr::bind_rows(lapply(sims,processsimulation)))
readr::write_csv(simsprocessed,file=file.path(resultspath,simsprocessedname))
cat("Finished processing simulations.\n")
cat("Wrote results to:",file.path(resultspath,simsprocessedname),"\n")


### Summarize Simulations ###

if (!exists("sims", .GlobalEnv)) {
  e <- new.env()
  load(file.path(resultspath,simresultsname), envir = e)
  sims <- e$sims
  rm(e)
  simsprocessed <- as_tibble(dplyr::bind_rows(lapply(sims,processsimulation)))
}

PLOTSIZE <- 16

simsummary <- simsprocessed %>%
  group_by(m, k, n) %>%
  summarize(covr0 = mean(beta0covrWald, na.rm = TRUE), 
            covr1 = mean(beta1covrWald, na.rm = TRUE),
            numsim = n(),
            numsuccess = sum(!is.na(beta0covrWald)),
            covr0SD = sqrt(mean(beta0covrWald, na.rm = TRUE) * (1 - mean(beta0covrWald, na.rm = TRUE)) / numsuccess),
            covr1SD = sqrt(mean(beta1covrWald, na.rm = TRUE) * (1 - mean(beta1covrWald, na.rm = TRUE)) / numsuccess),
            rmse0 = sqrt(mean(beta0bias^2)),
            rmse1 = sqrt(mean(beta0bias^2))
            ) %>%
  mutate(kf = factor(k), nf = factor(n))

# check the successes
xtabs( ~ numsuccess, data= simsummary)

k_labeller <- function(vl) paste0("Num. Points = ",vl)
n_labeller <- function(vl) paste0("Group Size = ",vl)

covrplt0 <- ggplot(simsummary, aes(x = m)) +
  theme_bw() +
  facet_grid(nf ~ kf, labeller = labeller(kf = k_labeller, nf = n_labeller)) +
  geom_hline(yintercept = .95, linetype = "dotted") +
  geom_line(aes(y = covr0)) +
  geom_line(aes(y = covr0 + 2 * covr0SD), linetype = "dashed") + 
  geom_line(aes(y = covr0 - 2 * covr0SD), linetype = "dashed") +
  geom_point(aes(y = covr0)) +
  theme(text = element_text(size = 20), axis.text.x = element_text(angle = 45)) +
  scale_x_continuous(labels = scales::comma_format(), breaks = m, trans = "log2") +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(title = expression("Empirical coverage of 95% Wald intervals for "~beta[0]),
       x = "Number of Groups", y = "Empirical Coverage") +
  geom_hline(yintercept = .95, linetype = "dotted")

ggsave(file = file.path(resultspath, "covrplot0.pdf"), plot = covrplt0, width = PLOTSIZE, height = PLOTSIZE)

covrplt1 <- ggplot(simsummary, aes(x = m)) +
  theme_bw() +
  facet_grid(nf ~ kf, labeller = labeller(kf = k_labeller, nf = n_labeller)) +
  geom_hline(yintercept = .95, linetype = "dotted") +
  geom_line(aes(y = covr1)) +
  geom_line(aes(y = covr1 + 2 * covr1SD), linetype = "dashed") + 
  geom_line(aes(y = covr1 - 2 * covr1SD), linetype = "dashed") +
  geom_point(aes(y = covr1)) +
  theme(text = element_text(size = 20), axis.text.x = element_text(angle = 45)) +
  scale_x_continuous(labels = scales::comma_format(), breaks = m, trans = "log2") +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(title = expression("Empirical coverage of 95% Wald intervals for "~beta[1]),
       x = "Number of Groups", y = "Empirical Coverage") +
  geom_hline(yintercept = .95, linetype = "dotted")

ggsave(file = file.path(resultspath, "covrplot1.pdf"), plot = covrplt1, width = PLOTSIZE, height = PLOTSIZE)

biasplt0 <- simsprocessed %>% 
  mutate(kf = factor(k), nf = factor(n), mf = factor(log(m, base = 2))) %>%
  ggplot(aes(x = mf)) +
  theme_bw() +
  facet_grid(nf ~ kf, labeller = labeller(kf = k_labeller, nf = n_labeller)) +
  geom_boxplot(aes(y = beta0bias)) + 
  scale_x_discrete(labels = scales::comma(m)) +
  theme(text = element_text(size = 20), axis.text.x = element_text(angle = 45)) +
  labs(title = expression("Empirical bias of approximate maximum likelihood estimates of "~beta[0]),
       x = "Number of Groups", y = "Empirical Bias") +
  geom_hline(yintercept = 0, linetype = "dotted")

ggsave(file = file.path(resultspath, "biasplt0.pdf"), plot = biasplt0, width = PLOTSIZE, height = PLOTSIZE)

biasplt1 <- simsprocessed %>% 
  mutate(kf = factor(k), nf = factor(n), mf = factor(log(m, base = 2))) %>%
  ggplot(aes(x = mf)) +
  theme_bw() +
  facet_grid(nf ~ kf, labeller = labeller(kf = k_labeller, nf = n_labeller)) +
  geom_boxplot(aes(y = beta1bias)) + 
  scale_x_discrete(labels = scales::comma(m)) +
  theme(text = element_text(size = 20), axis.text.x = element_text(angle = 45)) +
  labs(title = expression("Empirical bias of approximate maximum likelihood estimates of "~beta[1]),
       x = "Number of Groups", y = "Empirical Bias") +
  geom_hline(yintercept = 0, linetype = "dotted")

ggsave(file = file.path(resultspath, "biasplt1.pdf"), plot = biasplt1, width = PLOTSIZE, height = PLOTSIZE)

rmseplt0 <- ggplot(simsummary, aes(x = m)) +
  theme_bw() +
  facet_grid(nf ~ kf, labeller = labeller(kf = k_labeller, nf = n_labeller)) +
  geom_line(aes(y = rmse0)) +
  geom_point(aes(y = rmse0)) +
  theme(text = element_text(size = 20), axis.text.x = element_text(angle = 45)) +
  scale_x_continuous(labels = scales::comma_format(), breaks = m, trans = "log2") +
  labs(title = expression("Empirical RMSE of approximate maximum likelihood estimates of "~beta[0]),
       x = "Number of Groups", y = "Empirical RMSE")

ggsave(file = file.path(resultspath, "rmseplt0.pdf"), plot = rmseplt0, width = PLOTSIZE, height = PLOTSIZE)

rmseplt1 <- ggplot(simsummary, aes(x = m)) +
  theme_bw() +
  facet_grid(nf ~ kf, labeller = labeller(kf = k_labeller, nf = n_labeller)) +
  geom_line(aes(y = rmse1)) +
  geom_point(aes(y = rmse1)) +
  theme(text = element_text(size = 20), axis.text.x = element_text(angle = 45)) +
  scale_x_continuous(labels = scales::comma_format(), breaks = m, trans = "log2") +
  labs(title = expression("Empirical RMSE of approximate maximum likelihood estimates of "~beta[1]),
       x = "Number of Groups", y = "Empirical RMSE")

ggsave(file = file.path(resultspath, "rmseplt1.pdf"), plot = rmseplt1, width = PLOTSIZE, height = PLOTSIZE)

