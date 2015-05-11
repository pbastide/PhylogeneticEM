###############################################################################
## Bayou (see GitHub tutorial)
###############################################################################
# Data
library(bayou)
PATH <- "../Results/Chelonia/"
data(chelonia)
tree <- chelonia$phy
dat <- chelonia$dat
SE <- 0
# Prior
prior <- make.prior(tree, 
                    dists=list(dalpha="dlnorm",
                               dsig2="dlnorm",
                               dsb="dsb",
                               dk="cdpois",
                               dtheta="dnorm"),
                    param=list(dalpha=list(meanlog = -5, sdlog = 2.5),
                               dsig2=list(meanlog = 0, sdlog = 2),
                               dk=list(lambda=15, kmax=113),
                               dsb=list(bmax=1,prob=1),
                               dtheta=list(mean=3.5, sd=1.5)))
# MCMC
par(mfrow=c(2,3))
fit1 <- bayou.mcmc(tree,
                   dat,
                   SE=SE,
                   model="OU",
                   prior,
                   ngen=500000,
                   new.dir="../Results/Chelonia",
                   plot.freq=NULL,
                   ticker.freq=10000)
# chain
chain <- load.bayou(fit1, save.Rdata = TRUE, file = paste0(PATH, "bayou_chain_1.rds"), cleanup=TRUE)
chain <- set.burnin(chain, 0.3)
# sumary
bayou_sum <- summary(chain)
plot(chain)
# plot
par(mfrow=c(1,1))
plotSimmap.mcmc(tree, chain, burnin=0.3, circle=TRUE, fsize=0.4)
phenogram.density(tree, dat, chain=chain, burnin=0.3, pp.cutoff=0.3)
# Other chain for convergence
fit2 <- bayou.mcmc(tree,
                   dat, 
                   SE=SE,
                   model="OU", 
                   prior, 
                   ngen=500000, 
                   new.dir=TRUE, 
                   plot.freq=NULL, 
                   ticker.freq=10000)
chain2 <- load.bayou(fit2, save.Rdata = TRUE, file = paste0(PATH, "bayou_chain_2.rds"), cleanup=TRUE)
chain2 <- set.burnin(chain2, 0.3)
# Gelman's R statistic
RlnL <- gelman.R("lnL", chain1=chain, chain2=chain2, plot=TRUE, type="n", ylim=c(0.9, 2))
Ralpha <- gelman.R("alpha", chain1=chain, chain2=chain2, plot=TRUE, type="n", ylim=c(0.9, 2))
Rsig2 <- gelman.R("sig2", chain1=chain, chain2=chain2, plot=TRUE, type="n", ylim=c(0.9, 2))
# Position of shifts
L1 <- Lposterior(chain,tree, burnin=0.3)
L2 <- Lposterior(chain2,tree, burnin=0.3)
plot(L1$pp,L2$pp, xlim=c(0,1), ylim=c(0,1), xlab="Chain 1", ylab="Chain 2")
curve(1*x, add=TRUE, lty=2)
# Another chain for fun
fit3 <- bayou.mcmc(tree,
                   dat, 
                   SE=SE,
                   model="OU", 
                   prior, 
                   ngen=500000, 
                   new.dir=TRUE, 
                   plot.freq=NULL, 
                   ticker.freq=10000)
chain3 <- load.bayou(fit3, save.Rdata = TRUE, file = paste0(PATH, "bayou_chain_3.rds"), cleanup=TRUE)
chain3 <- set.burnin(chain3, 0.3)
# Combine chains
chains <- combine.chains(chain, chain2, burnin.prop=0.3)
chains <- set.burnin(chains, 0)
# plot
par(mfrow=c(1,1))
plotSimmap.mcmc(tree, chains, burnin=0.3, circle=TRUE, fsize=0.4)
phenogram.density(tree, dat, chain=chains, burnin=0.3, pp.cutoff=0.3)
# marginal likelihood
ss <- steppingstone(Bk=seq(0,1,length.out=5),
                    chains, 
                    tree, 
                    dat, 
                    SE=SE, 
                    prior=prior, 
                    new.dir=TRUE,
                    ngen=10000)
ss <- set.burnin(ss, 0.3)
ss
plot(ss)

# Summary results
OU_bayou <- list(Nbr_shifts = median(chains$k),
                 Nbr_regimes = median(chains$ntheta),
                 lnL = median(chains$lnL),
                 MlnL = ss$lnr,
                 alpha = median(chains$alpha),
                 half_life = median(log(2)/chains$alpha),
                 sigma = median(chains$sig2),
                 gamma = median(chains$sig2/(2*chains$alpha)))

save.image(paste0(PATH, "chelonia_bayou_res.RData"))
