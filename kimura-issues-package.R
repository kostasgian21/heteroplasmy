library(ggplot2)
library(gridExtra)
library(metR)
library(kimura)
library(heteroplasmy)

# problem 1

h = c(0,1)
hbar = mean(h)
s2 = var(h)
n = length(h)

nhat = hbar*(1-hbar)/s2
nhat

h0 = 0.1
s2 = sum((h-h0)**2)/(n-1)
s2

nhat = h0*(1-h0)/s2
nhat

# problem 2

h = c(0,1)
hbar = mean(h)
s2 = var(h)
n = length(h)

mu2 = sum((h-hbar)**2)/n
mu4 = sum((h-hbar)**4)/n
D4 = (n-1)/n**3 * ((n**2-3*n+3)*mu4 + 3*(2*n-3)*mu2**2)
se.s2 = sqrt(1/n*(D4 - (n-3)/(n-1)*s2**2))
se.s2

h0 = 0.1
s2 = sum((h-h0)**2)/(n-1)
mu2 = sum((h-h0)**2)/n
mu4 = sum((h-h0)**4)/n
D4 = (n-1)/n**3 * ((n**2-3*n+3)*mu4 + 3*(2*n-3)*mu2**2)
se.s2 = sqrt(1/n*(D4 - (n-3)/(n-1)*s2**2))
se.s2

# problem 3

mom.fit = estimate_parameters(h)
mom.fit
p = mom.fit[1]
b = mom.fit[2]
D4 = s2*((p-0.5)**2 * (3*(1-b-b**2)+b**3+b**4+b**5) + 0.25*(1- (b+b**2+b**3+b**4+b**5)/5))
se.s2 = sqrt(1/n*(D4 - (n-3)/(n-1)*s2**2))

# problem 4

h = rep(c(0, 0.9), 4)

test_kimura(h, num_MC = 10000)
ks.fit = estimate_parameters_ks(h)
test_kimura_param(h, ks.fit[1], ks.fit[2], num_MC = 10000)

# solution

h = c(0,1)
estimate_parameters_ml(h)

best = optim(c(0.5, 0.5), kimura_neg_loglik, h=h, h0=F, hessian=T)

best$b.hat = transfun(best$par[1])
best$h0.hat = transfun(best$par[2])


# LRT

h1 = list( c(0.2,0.4), c(0.6,0.7) )
h2 = list( c(0.1,0.7), c(0.3,1.0) )
kimura_lrt(h1, h2)
