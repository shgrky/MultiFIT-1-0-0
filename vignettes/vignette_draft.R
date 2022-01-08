set.seed(1)
# Generate data for two random vectors, each of dimension 2, 300 observations:
n = 300
x = matrix(0, ncol = 2, nrow = n)
y = matrix(0, ncol = 2, nrow = n)

# x1 and y1 are i.i.d Normal(0,1):
x[ , 1] = rnorm(n)
y[ , 1] = rnorm(n)

# x2 is a Uniform(0,1):  
x[,2]=runif(n)

# and y2 is depends on x2 as a noisy sine function:
y[,2]=sin(5*pi*x[,2]) + 0.6*rnorm(n)


# Full: -----
source("projects/MultiFit/MultiFIT-1-1-0/vignettes/helpers.R")
source("projects/MultiFit/MultiFIT-1-1-0/vignettes/multiple_testing.R")

library(MultiFIT)
fit = MultiFIT(x=x, y=y, verbose = TRUE, p.adjust.methods = c("H", "Hcorrected"),
               return.all.pvs = TRUE, compute.all.holm = TRUE)
fit$p.values.holistic
fit$p.values.resolution.specific
fit$res.by.res.pvs
fit$all.pvs

fit2 = MultiFit::multiFit(x=x, y=y, verbose = TRUE, p.adjust.methods = c("H", "Hcorrected"))

compute.all.holm = F
verbose = F
fit3 = MultiFitByResMH::multiFit_apply_stopping_MH(x=x, y=y, verbose = TRUE)
by_res_multiple_adjustment(fit3)

by_res_multiple_adjustment(fit2)

fit4 = MultiFitByResMHrknew::multiFit_apply_stopping_MH_rk_new(x=x, y=y, verbose = TRUE,
                                                               p.adjust.methods = c("H", "Hcorrected"))

fit$p.values.holistic
fit2$p.values
fit3$p.values
fit4$p.values

fit$p.values.resolution.specific
by_res_multiple_adjustment(fit2)
by_res_multiple_adjustment(fit3)
by_res_multiple_adjustment(fit4)


# Great: MultiFIT, multiFit and multiFit_apply_stopping_MH agree on full algorithm.


# Full, return pvs -----
library(MultiFIT)
fit = MultiFIT(x=x, y=y, verbose = TRUE, p.adjust.methods = c("H", "Hcorrected"),
               return.all.pvs = TRUE)

fit2 = MultiFit::multiFit(x=x, y=y, verbose = TRUE, p.adjust.methods = c("H", "Hcorrected"),
                          return.all.pvs = TRUE)

compute.all.holm = F
verbose = F
fit3 = MultiFitByResMH::multiFit_apply_stopping_MH(x=x, y=y, verbose = TRUE,
                                                   return.all.pvs = TRUE)

fit4 = MultiFitByResMHrknew::multiFit_apply_stopping_MH_rk_new(x=x, y=y, verbose = TRUE,
                                                               p.adjust.methods = c("H", "Hcorrected"),
                                                               return.all.pvs = TRUE)

fit$p.values.holistic
fit2$p.values
fit3$p.values
fit4$p.values

fit$p.values.resolution.specific
by_res_multiple_adjustment(fit2)
by_res_multiple_adjustment(fit3)
by_res_multiple_adjustment(fit4)

attributes(fit$all.pvs) <- NULL
attributes(fit2$pvs) <- NULL
attributes(fit3$pvs) <- NULL
attributes(fit4$pvs) <- NULL
all.equal(fit$all.pvs, fit2$pvs)
all.equal(fit$all.pvs, fit3$pvs) # this has MH
all.equal(fit$all.pvs, fit4$pvs)


# rank, M2: -----
fit = MultiFIT(x=x, y=y, verbose = TRUE, ranking.approximation = TRUE, M = 2,
               apply.stopping.rule = FALSE,
               p.adjust.methods = c("H", "Hcorrected"))
fit$p.values.holistic
fit$p.values.resolution.specific

fit2 = MultiFit::multiFit(x=x, y=y, verbose = TRUE, rk = TRUE, M = 2,
                          p.adjust.methods = c("H", "Hcorrected"))

compute.all.holm = FALSE
verbose = FALSE
fit3 = 
  MultiFitByResMHrknew::multiFit_apply_stopping_MH_rk_new(
    x=x, y=y, verbose = TRUE, rk = TRUE, apply_stopping_rule = FALSE, M = 2,
    p.adjust.methods = c("H", "Hcorrected"))

fit4 = 
  MultiFitByResMH::multiFit_apply_stopping_MH(
  x=x, y=y, verbose = TRUE, rk = TRUE, apply_stopping_rule = FALSE, M = 2,
  p.adjust.methods = c("H", "Hcorrected", "MH"))

fit$p.values.holistic
fit2$p.values
fit3$p.values
fit4$p.values

# Great:
# MultiFIT, multiFit_apply_stopping_MH_rk_new agree: they apply ranking AND p*
# They disagree with multiFit and multiFit_apply_stopping_MH because the latter
# pair only applies rank **without** cutting by p* as well.



# rank, M2, return pvs: -----
fit = MultiFIT(x=x, y=y, verbose = TRUE, ranking.approximation = TRUE, M = 2,
               apply.stopping.rule = FALSE,
               p.adjust.methods = c("H", "Hcorrected"),
               return.all.pvs = TRUE)
fit$p.values.holistic
fit$p.values.resolution.specific

fit2 = MultiFit::multiFit(x=x, y=y, verbose = TRUE, rk = TRUE, M = 2,
                          p.adjust.methods = c("H", "Hcorrected"),
                          return.all.pvs = TRUE)

compute.all.holm = FALSE
verbose = FALSE
fit3 = 
  MultiFitByResMHrknew::multiFit_apply_stopping_MH_rk_new(
    x=x, y=y, verbose = TRUE, rk = TRUE, apply_stopping_rule = FALSE, M = 2,
    p.adjust.methods = c("H", "Hcorrected"),
    return.all.pvs = TRUE)

fit4 = 
  MultiFitByResMH::multiFit_apply_stopping_MH(
    x=x, y=y, verbose = TRUE, rk = TRUE, apply_stopping_rule = FALSE, M = 2,
    p.adjust.methods = c("H", "Hcorrected", "MH"),
    return.all.pvs = TRUE)

fit$p.values
fit2$p.values
fit3$p.values
fit4$p.values

# Great:
# MultiFIT, multiFit_apply_stopping_MH_rk_new agree: they apply ranking AND p*
# They disagree with multiFit and multiFit_apply_stopping_MH because the latter
# pair only applies rank **without** cutting by p* as well.

# rank, M2:
fit = MultiFIT(x=x, y=y, verbose = TRUE, ranking.approximation = TRUE, M = 2,
               apply.stopping.rule = FALSE,
               p.adjust.methods = c("H", "Hcorrected"))
fit$p.values

fit2 = MultiFit::multiFit(x=x, y=y, verbose = TRUE, rk = TRUE, M = 2,
                          p.adjust.methods = c("H", "Hcorrected"))

compute.all.holm = FALSE
verbose = FALSE
fit3 = 
  MultiFitByResMHrknew::multiFit_apply_stopping_MH_rk_new(
    x=x, y=y, verbose = TRUE, rk = TRUE, apply_stopping_rule = FALSE, M = 2,
    p.adjust.methods = c("H", "Hcorrected"))

fit4 = 
  MultiFitByResMH::multiFit_apply_stopping_MH(
    x=x, y=y, verbose = TRUE, rk = TRUE, apply_stopping_rule = FALSE, M = 2,
    p.adjust.methods = c("H", "Hcorrected", "MH"))

fit$p.values
fit2$p.values
fit3$p.values
fit4$p.values

attributes(fit$pvs) <- NULL
attributes(fit2$pvs) <- NULL
attributes(fit3$pvs) <- NULL
attributes(fit4$pvs) <- NULL
# all.equal(fit$pvs, fit2$pvs)
all.equal(fit$all.pvs, fit3$pvs)
# all.equal(fit$pvs, fit4$pvs)

# Great:
# MultiFIT, multiFit_apply_stopping_MH_rk_new agree: they apply ranking AND p*
# They disagree with multiFit and multiFit_apply_stopping_MH because the latter
# pair only applies rank **without** cutting by p* as well.

# stopping, no ranking ----
fit = MultiFIT(x=x, y=y, verbose = TRUE, apply.stopping.rule = TRUE, return.all.pvs = FALSE,
               p.adjust.methods = c("H", "Hcorrected"))

fit$p.values.resolution.specific

fit2 = MultiFitByResMH::multiFit_apply_stopping_MH(
  x=x, y=y, verbose = TRUE, apply_stopping_rule = TRUE, return.all.pvs = FALSE,
  p.adjust.methods = c("H", "Hcorrected", "MH"))

fit3 = MultiFitByResMHrknew::multiFit_apply_stopping_MH_rk_new(
  x=x, y=y, verbose = TRUE, apply_stopping_rule = TRUE, return.all.pvs = FALSE,
  p.adjust.methods = c("H", "Hcorrected"))

fit$p.values.holistic
fit$p.values.resolution.specific
fit$res.by.res.pvs
fit2$p.values
fit3$p.values

# Great, they all agree for stopping rule without ranking!!!

# stopping, no ranking, return pvs ----
fit = MultiFIT(x=x, y=y, verbose = TRUE, apply.stopping.rule = TRUE,
               p.adjust.methods = c("H", "Hcorrected"),
               return.all.pvs = TRUE)

fit2 = MultiFitByResMH::multiFit_apply_stopping_MH(
  x=x, y=y, verbose = TRUE, apply_stopping_rule = TRUE,
  p.adjust.methods = c("H", "Hcorrected", "MH"),
  return.all.pvs = FALSE)

fit3 = MultiFitByResMHrknew::multiFit_apply_stopping_MH_rk_new(
  x=x, y=y, verbose = TRUE, apply_stopping_rule = TRUE,
  p.adjust.methods = c("H", "Hcorrected"),
  return.all.pvs = FALSE)

fit$p.values.resolution.specific
fit2$p.values
fit3$p.values

# stopping, ranking ----
fit = MultiFIT(x=x, y=y, verbose = TRUE, apply.stopping.rule = TRUE, return.all.pvs = FALSE,
               p.adjust.methods = c("H", "Hcorrected"),
               ranking.approximation = TRUE, M = 5)

fit2 = MultiFitByResMH::multiFit_apply_stopping_MH(
  x=x, y=y, verbose = TRUE, apply_stopping_rule = TRUE, return.all.pvs = FALSE,
  p.adjust.methods = c("H", "Hcorrected", "MH"),
  rk = TRUE, M = 5)

fit3 = MultiFitByResMHrknew::multiFit_apply_stopping_MH_rk_new(
  x=x, y=y, verbose = TRUE, apply_stopping_rule = TRUE, return.all.pvs = FALSE,
  p.adjust.methods = c("H", "Hcorrected"),
  rk = TRUE, M = 5)

fit$p.values.resolution.specific
fit2$p.values
fit3$p.values

# Great, stopping with ranking is the same.






# Data may also be transferred to the function as a single list:
library(MultiFIT)
xy = list(x=x,y=y)
fit = MultiFIT(xy, verbose=TRUE)
MultiFit::multiFit(xy, verbose=TRUE)

# And plot a DAG representation of the ranked tests:
library(png)
library(qgraph)
multiTree(xy=xy, fit=fit, filename="first_example")

fit1 = MultiFIT(xy, p_star = 0.1, verbose=TRUE)
multiSummary(xy=xy, fit=fit1, alpha=0.005, plot.tests=FALSE)

fit2 = MultiFit::multiFit(xy, p_star = 0.1, verbose=TRUE)
MultiFit::multiSummary(xy=xy, fit=fit2, alpha=0.005, plot.tests=FALSE)

# 1. set p_star=Inf, running through all tables up to the maximal resolution
# which by default is set to log2(n/100):
ex1 = MultiFIT(xy, p_star = 1, verbose = TRUE)
ex1 = MultiFit::multiFit(xy, p_star = 1, verbose = TRUE)

# 2. set both p_star=1 and the maximal resolution R_max=Inf.
# In this case, the algorithm will scan through higher and higher resolutions,
# until there are no more tables that satisfy the minimum requirements for 
# marginal totals: min.tbl.tot, min.row.tot and min.col.tot (whose default values 
# are presented below):
ex2 = MultiFIT(xy, p_star = 1, R_max=Inf,
               min.tbl.tot = 25L, min.row.tot = 10L, min.col.tot = 10L)

ex2$res.by.res.pvs
ex2$p.values.holistic
ex2$p.values.resolution.specific

ex2 = MultiFit::multiFit(xy, p_star = 1, R_max=Inf,
               min.tbl.tot = 25L, min.row.tot = 10L, min.col.tot = 10L)
ex2$p.values

# 3. set smaller minimal marginal totals, that will result in testing 
# even more tables in higher resolutions:
ex3 = MultiFIT(xy, p_star = 1, R_max=Inf,
               min.tbl.tot = 10L, min.row.tot = 4L, min.col.tot = 4L)
ex3$p.values.holistic

ex3 = MultiFit::multiFit(xy, p_star = 1, R_max=Inf,
               min.tbl.tot = 10L, min.row.tot = 4L, min.col.tot = 4L)

ex3$p.values

# Generate data for two random vectors, each of dimension 2, 800 observations:
n=800
x = matrix(0, ncol=2, nrow=n)
y = matrix(0, ncol=2, nrow=n)

# x1, x2 and y1 are i.i.d Normal(0,1):
x[,1]=rnorm(n)
x[,2]=rnorm(n)
y[,1]=rnorm(n)

# y2 is i.i.d Normal(0,1) on most of the space:
y[,2]=rnorm(n)
# But is linearly dependent on x2 in a small portion of the space:
w=rnorm(n)
portion.of.space = x[,2]>0 & x[,2]<0.7 & y[,2]>0 & y[,2]<0.7
y[portion.of.space,2] = x[portion.of.space,2]+(1/12)*w[portion.of.space]
xy.local = list(x=x, y=y)

fit.local = MultiFIT(xy=xy.local, R_star=4, verbose=TRUE)
multiSummary(xy=xy.local, fit=fit.local, plot.margin=TRUE, pch="`")

# Marginal signal:

# Generate data for two random vectors, each of dimension 3, 800 observations:
# Generate data for two random vectors, each of dimension 3, 800 observations:
set.seed(3212021)
n=800
x = matrix(0, ncol=3, nrow=n)
y = matrix(0, ncol=3, nrow=n)

# x1, x2, y1 and y2 are all i.i.d Normal(0,1)
x[,1]=rnorm(n)
x[,2]=rnorm(n)
y[,1]=rnorm(n)
y[,2]=rnorm(n)

# x3 and y3 form a noisy circle:
theta = runif(n,-pi,pi)
x[,3] = cos(theta) + 0.1*rnorm(n)
y[,3] = sin(theta) + 0.1*rnorm(n)

# par(mfrow=c(3,3))
# # par(mgp=c(0,0,0))
# par(mar=c(1.5,1.5,0,0))
# for (i in 1:3) {
#   for (j in 1:3) {
#     plot(x[,i],y[,j], col="black", pch=20, xlab=paste0("x",i), ylab=paste0("y",j),
#          xaxt="n", yaxt="n")
#   }
# }


# And now rotate the circle:
phi = pi/4
rot.mat = matrix(c(cos(phi), -sin(phi),  0,
                   sin(phi),  cos(phi),  0,
                   0,         0,         1), nrow=3, ncol=3)
xxy = t(rot.mat%*%t(cbind(x[,2],x[,3],y[,3])))

x.rtt = matrix(0, ncol=3, nrow=n)
y.rtt = matrix(0, ncol=3, nrow=n)

x.rtt[,1] = x[,1]
x.rtt[,2] = xxy[,1]
x.rtt[,3] = xxy[,2]
y.rtt[,1] = y[,1]
y.rtt[,2] = y[,2]
y.rtt[,3] = xxy[,3]

# par(mfrow = c(3,3))
# par(mgp = c(0,0,0))
# par(mar = c(1.5,1.5,0,0))
# for (i in 1:3) {
#   for (j in 1:3) {
#     plot(x.rtt[,i],y.rtt[,j], col = "black", pch = 20, xlab = paste0("x", i),
#          ylab = paste0("y", j), xaxt = "n", yaxt = "n")
#   }
# }

xy.rtt.circ = list(x = x.rtt, y = y.rtt)

fit.rtt.circ = 
  MultiFIT(xy = xy.rtt.circ, R_star = 2, verbose = TRUE, p.adjust.methods = "Hcorrected")

par(mfrow = c(1,2))
par(mar = c(4,4,3,3))
multiSummary(xy = xy.rtt.circ, fit = fit.rtt.circ, alpha = 0.001)
