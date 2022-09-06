library(ggplot2)
library(gridExtra)
library(metR)
library(kimura)
library(heteroplasmy)

grid = 20 # cell width for heatmaps, default 20
nbin = 20 # number of bins for for histograms, default 20

# the wrapper loop goes through different (synthetic and experimental) data sources
# the publication uses indices 1 (0,0.9) repeated 4 times; 4 (Freyer et al.); 9 (0.0.9) repeated 50 times
# others may be informative
# there are some index-specific tweaks to the plotting code below because different datasets have different scales etc
for(expt in c(1,4,9)) {
  if(expt == 1)        { h = rep(c(0, 0.9), 4)   }
  else if(expt == 2)   { h = rep(c(0.1, 0.5), 8) }
  else if(expt == 3)   { h = c(60, 68, 72, 78, 62, 52, 71, 59, 59, 64, 64, 65, 76, 56, 71, 74, 72, 63, 64, 64, 64, 65, 63, 42)/100 } # Zhang et al. P8A
  else if(expt == 4)   { h = c(0.61, 0.76, 0.47, 0.4, 0.54, 0.56, 0.41, 0.71, 0.66, 0.5, 0.4, 0.45, 0.69, 0.66, 0.64, 0.5, 0.66, 0.75, 0.77, 0.75, 0.67, 0.63, 0.74, 0.75, 0.78, 0.67, 0.65, 0.75, 0.64, 0.62, 0.37, 0.57, 0.69, 0.47, 0.51, 0.61, 0.77, 0.55, 0.53, 0.57, 0.51, 0.66, 0.67, 0.68, 0.66, 0.72, 0.79, 0.43, 0.75, 0.6, 0.72, 0.73, 0.73, 0.72, 0.69, 0.71, 0.75, 0.65, 0.56, 0.72, 0.72, 0.68, 0.64, 0.73, 0.69, 0.71, 0.78, 0.75, 0.77, 0.75, 0.75, 0.69, 0.65, 0.62, 0.6, 0.53, 0.66, 0.69, 0.42, 0.69, 0.58, 0.62, 0.65, 0.73, 0.77, 0.61, 0.73, 0.76, 0.69, 0.71, 0.64, 0.77, 0.72, 0.72, 0.64, 0.71, 0.71, 0.68, 0.72, 0.69, 0.74, 0.63, 0.68, 0.69, 0.53, 0.71, 0.75, 0.68, 0.56, 0.65, 0.59, 0.77, 0.71, 0.64, 0.54, 0.79, 0.68, 0.68, 0.73, 0.73, 0.74, 0.67, 0.73, 0.77, 0.6, 0.59, 0.74, 0.77, 0.66, 0.68, 0.69, 0.64, 0.66, 0.73, 0.7, 0.59, 0.62, 0.61, 0.64, 0.59, 0.73, 0.7, 0.46, 0.51, 0.59, 0.57, 0.65, 0.57, 0.68, 0.56, 0.56, 0.64, 0.49, 0.75, 0.32, 0.75, 0.74, 0.77, 0.78, 0.77, 0.76, 0.71, 0.73, 0.75, 0.68, 0.62, 0.71, 0.69, 0.71, 0.62, 0.67, 0.61, 0.72, 0.75, 0.62, 0.72, 0.74, 0.68, 0.65, 0.78, 0.76, 0.47, 0.74, 0.89, 0.6, 0.64, 0.69, 0.82, 0.78, 0.65, 0.69, 0.68, 0.7, 0.79, 0.67, 0.69, 0.71, 0.7, 0.53, 0.72, 0.62, 0.66, 0.64, 0.78, 0.74, 0.68, 0.47, 0.71, 0.65, 0.72, 0.8, 0.6, 0.79, 0.75, 0.64, 0.7, 0.72, 0.62, 0.52, 0.67, 0.58, 0.65, 0.71, 0.6, 0.59, 0.38, 0.58, 0.77, 0.72, 0.58, 0.65, 0.68, 0.45, 0.79, 0.84, 0.52, 0.38, 0.43, 0.71, 0.56, 0.6, 0.63, 0.73, 0.69, 0.52, 0.57, 0.62, 0.77, 0.78, 0.71, 0.8, 0.75, 0.76, 0.67, 0.72, 0.71, 0.6, 0.69, 0.77, 0.57, 0.69, 0.58, 0.61, 0.68, 0.72, 0.61, 0.48, 0.62, 0.67, 0.77, 0.77, 0.69, 0.74, 0.71, 0.77, 0.64, 0.6, 0.69, 0.57, 0.61, 0.76, 0.67, 0.84, 0.83, 0.81, 0.81, 0.81, 0.57, 0.76, 0.71, 0.64, 0.55) } # Freyer et al. > 70
  else if(expt == 5)   { h=c(rep(0,8),rep(5,35),rep(15,20),rep(25,10),rep(35,5),rep(45,4))/100 }
  else if(expt == 6)   { h=c(rep(0,8),runif(35, min=1,max=10),runif(20,min=10,max=20),runif(10,min=20,max=30),runif(5,min=30,max=40),runif(4,min=40,max=50))/100 }
  else if(expt == 7)   { h = c(rep(0, 15), ((1:40)/200)**2) }
  else if(expt == 8)   { h = round(c(rep(0, 15), ((1:40)/200)**2)*100)/100 }
  else if(expt == 9)   { h = rep(c(0, 0.9), 50)   }
  
  # first, estimate parameters using our various approaches
  mom.best = estimate_parameters(h)
  mom.p = mom.best[1]
  mom.b = mom.best[2]
  ml.best = estimate_parameters_ml(h)
  ml.p = ml.best[1]
  ml.b = ml.best[2]
  mks.best = estimate_parameters_ks(h)
  mks.p = mks.best[1]
  mks.b = mks.best[2]
  
  # use the Monte Carlo KS test with these different parameterisations
  t1 = test_kimura(h, num_MC = 10000)
  t1.a = test_kimura_par(h, mom.p, mom.b, num_MC = 10000)
  t2 = test_kimura_par(h, ml.p, ml.b, num_MC = 10000)
  t3 = test_kimura_par(h, mks.p, mks.b, num_MC = 10000)
  # first two vals should be identical, final value should be higher that all others
  c(t1$p.value, t1.a$p.value, t2$p.value, t3$p.value)
  
  ## PLOT 1 -- PDFs
  # produce the PDFs for these parameterisations
  hset = (0:nbin)/nbin
  pset1 = dkimura(hset, mom.p, mom.b)/nbin
  pset2 = dkimura(hset, ml.p, ml.b)/nbin
  pset3 = dkimura(hset, mks.p, mks.b)/nbin
  # the first and last elements correspond to point mass at 0 and 1 and need to be scaled to be compared with the other (integrated density) values
  pset1[1] = pset1[1]*nbin
  pset1[nbin+1] = pset1[nbin+1]*nbin
  pset2[1] = pset2[1]*nbin
  pset2[nbin+1] = pset2[nbin+1]*nbin
  pset3[1] = pset3[1]*nbin
  pset3[nbin+1] = pset3[nbin+1]*nbin
  
  # build data frame with these PDFs for plotting
  df = data.frame(h=hset,Fit="Moments",p=pset1)
  df = rbind(df,data.frame(h=hset,Fit="Max lik",p=pset2))
  df = rbind(df,data.frame(h=hset,Fit="Min KS",p=pset3))
  df$Fit = factor(df$Fit, levels=c("Moments", "Max lik", "Min KS"))
  
  # index-specific positioning of data scatter (below histograms)
  if(expt == 4) {yoff = -0.02}
  else {yoff = -0.05}
  
  # plot histograms and data scatter
  g1 = ggplot(df, aes(x=h, y=p, fill = Fit, color=Fit)) +
    geom_col(position="dodge") +
    scale_fill_manual(values = c("#FF8888", "#AA4444", "#8888FF")) +
    scale_color_manual(values = c("#FF8888", "#AA4444", "#8888FF")) +
    geom_jitter(data = data.frame(x = h, y = -0.05), color = "#000000", fill = "#000000", aes(x=x, y=y), width=0, height = 0.01) +
    theme_classic() +
    theme(legend.position="none", axis.text=element_text(size=12), axis.title=element_text(size=14))

  ## PLOT 2 -- likelihood surface
  # index-specific plotting of likelihood surface (only for experiment 1)
if(expt == 1) {
  # we're going to make a data frame storing likelihood values at each point in a parameter grid
  hobs = h
  ll.df = data.frame()
  minp = max(1/grid,        min(mom.p, ml.p, mks.p)-0.1)*grid
  maxp = min((grid-1)/grid, max(mom.p, ml.p, mks.p)+0.1)*grid
  minb = max(1/grid,	    min(mom.b, ml.b, mks.b)-0.1)*grid
  maxb = min((grid-1)/grid, max(mom.b, ml.b, mks.b)+0.1)*grid
  
  # lazily loop through the grid points and add likelihood values at each point
  for(h0 in (minp:maxp)/grid) {
    for(b in (minb:maxb)/grid) {
      ll = -kimura_neg_loglik(c(invtransfun(b), invtransfun(h0)), hobs)
      ll.df = rbind(ll.df, data.frame(p=h0, b=b, ll=ll))
    }
  }
  # some values diverge; replace with a cutoff
  ll.df$ll[abs(ll.df$ll) == Inf] = 1000
  
  # data frame for specific parameterisations
  models = data.frame(x=c(mom.p, ml.p, mks.p), y=c(mom.b, ml.b, mks.b), label=c("Moments", "Max lik", "Min KS"))
  
  # produce plot with contour map of surface and parameter points
  g2a = ggplot(ll.df, aes(x = p, y = b, z= ll)) +
    geom_contour(bins=30) +
    geom_text_contour(skip=0) +
    geom_point(data=models, aes(x=x, y=y)) +
    geom_text(data=models, aes(x=x+0.025,y=y,label=label)) +
    theme_classic()
    
  # compute specific parameterisation likelihoods for legend
  ml.ml = -kimura_neg_loglik(c(invtransfun(ml.b), invtransfun(ml.p)), hobs)
  mom.ml = -kimura_neg_loglik(c(invtransfun(mom.b), invtransfun(mom.p)), hobs)
  mks.ml = -kimura_neg_loglik(c(invtransfun(mks.b), invtransfun(mks.p)), hobs)
  
  # add legend giving specific values
  g2a.lab = g2a + annotate("label", x = minp/grid+0.1, y = minb/grid+0.1, label = paste(c("MoM MlL = ", round(mom.ml, digits=2), "\nML MlL = ", round(ml.ml, digits=2), "\nMin KS MlL = ", round(mks.ml, digits=2)), collapse=""))
  }
  
  ## PLOT 3 -- plot CDFs for data and distributions under different parameterisations
  # get CDFs for different parameterisations
  cdf_kimura.mom <- .pkimura_full(mom.p, mom.b)
  cdf_kimura.ml <- .pkimura_full(ml.p, ml.b)
  cdf_kimura.mks <- .pkimura_full(mks.p, mks.b)
  
  # empirical CDF
  ecdf_h <- (stats::ecdf(h))(seq(0, 1, 1e-04))
  
  # output KS distances
  c(max(abs(ecdf_h - cdf_kimura.mom)), max(abs(ecdf_h - cdf_kimura.ml)), max(abs(ecdf_h - cdf_kimura.mks)))
  
  # build data frame with all CDFs
  cdf.df = data.frame(h=rep(seq(0, 1, 1e-04),4), model=c(rep("Max lik", length(ecdf_h)), rep("Moments", length(cdf_kimura.mom)), rep("Min KS", length(cdf_kimura.mks)), rep("Data", length(cdf_kimura.ml))), CDF = c(cdf_kimura.ml, cdf_kimura.mom,  cdf_kimura.mks,  ecdf_h))
  cdf.df$model = factor(cdf.df$model, levels = c("Moments", "Max lik", "Min KS", "Data"))
  
  # plot CDFs
  g3 = ggplot(cdf.df, aes(x=h, y=CDF, colour=model)) +
    geom_point() +
    scale_color_manual(values = c("#FF8888", "#AA4444", "#8888FF", "#000000")) +
    theme_classic() +
    labs(colour="CDF") + 
    theme(legend.position="none", axis.text=element_text(size=12), axis.title=element_text(size=14))
  
  # create and add label giving statistics
  lab.3 = paste(c("MoM mKS = ", round(max(abs(ecdf_h - cdf_kimura.mom)), digits=2), ", p = ", round(t1$p.value, digits=2), "\nML mKS = ", round(max(abs(ecdf_h - cdf_kimura.ml)), digits=2), ", p = ", round(t2$p.value, digits=2), "\nMin KS mKS = ", round(max(abs(ecdf_h - cdf_kimura.mks)), digits=2), ", p = ", round(t3$p.value, digits=2)), collapse="")
  g3.lab = g3 + annotate("label", x = 0.25, y = 0.9, label = lab.3)
  
  ## PLOT 4 -- KS distance surface
  
  # build data frame to store KS distance values at different parameterisations
  ks.df = data.frame()
  # index-specific tweak of parameter range fo
  if(expt == 4) { h0range = 0.6 + 0.1*(1:(grid-1))/grid ; brange = 0.95 + 0.025*(1:(grid-1))/grid; }
  else { h0range = (1:(grid-1))/grid; brange = (1:(grid-1))/grid; }
  
  # lazily loop through parameter range and conpute KS distance for each
  for(h0 in h0range) {
    for(b in brange) {
      cdf_kimura <- .pkimura_full(h0, b)
      ks = max(abs(ecdf_h - cdf_kimura)) 
      ks.df = rbind(ks.df, data.frame(p=h0, b=b, ks=ks))
    }
  }
  
    # index-specific tweak for label offset
    if(expt == 4) {dy = -0.001}
    else { dy = 0.05 }
    
    # plot KS distance surface 
    g4 = ggplot(ks.df, aes(x = p, y = b, fill= ks)) + geom_tile() + scale_fill_gradient(low="#000088", high="#FFFFFF") + labs(fill="KS dist") + theme_classic()
    
    # add specific parameterisations and labels
    g4.lab = g4+ geom_point(data=models, aes(x=x, y=y), color="#FFFFFF") +
      geom_text(data=models, aes(x=x,y=y+dy,label=label), color="#FFFFFF") +
      theme(axis.text=element_text(size=12), axis.title=element_text(size=14))
    
    
    
  # now output selected plots to files
  # file with all labels for reference
  myres = 2
  png(paste(c("kimura-plots-red-", expt, ".png"), collapse=""), width=1200*myres, height=400*myres, res=72*myres)
  grid.arrange(g1, g3.lab, g4.lab, nrow=1)
  dev.off()
  # index-specific plot of likelihood surface
  if(expt == 1) {
  png(paste(c("kimura-plots-ml-", expt, ".png"), collapse=""), width=400*myres, height=400*myres, res=72*myres)
  grid.arrange(g2a.lab, nrow=1)
  dev.off()
  }
  # plot without several annotations for clarity
  png(paste(c("kimura-plots-pub-", expt, ".png"), collapse=""), width=800*myres, height=300*myres, res=72*myres)
  grid.arrange(g1, g3+theme(legend.position="none", axis.text.x=element_text(size=12), axis.text.y=element_text(size=12)), g4.lab, nrow=1)
  dev.off()

  # output stats for reference
  
  c(t1$p.value, t1.a$p.value, t2$p.value, t3$p.value)
}

