cl <- makeCluster(detectCores())

size = meanE = sigE = meanZ = sigZ =  biasE = biasZ = c()
EdistPars = ZdistPars = matrix(NA,ncol=3,nrow=length(setList))
colnames(EdistPars) = colnames(ZdistPars) = c("df","mu","sigma")

for(i in seq_along(setList)) {
  D2 = dataList[[paste0(setList[i],'_cal')]]
  E  = D2$E
  uE = D2$uE
  Z  = E / uE
  size[i] = M = length(uE)

  # E
  mu  = mean(E)
  sig = sd(E)/sqrt(length(E))
  meanE[i] = ErrViewLib::prettyUnc(mu,sig)

  bs = nptest::np.boot(
    x = E, statistic = sd,
    R = nBoot, level = 0.95, method = "bca",
    parallel = TRUE, cl = cl)
  mus  = bs$t0
  sigs = bs$se
  sigE[i] = ErrViewLib::prettyUnc(mus,sigs)
  biasE[i] = abs(mu/mus)

  ## Fit by Stutent's t to get effective df
  fit.t<-fitdistrplus::fitdist(
    E, "t_ls", start = list(df=20,mu=mean(E),sigma=sd(E)),
    keepdata = FALSE)
  pars = summary(fit.t)$estimate
  EdistPars[i,] = summary(fit.t)$estimate

  # Z
  mu  = mean(Z)
  sig = sd(Z)/sqrt(length(Z))
  meanZ[i] = ErrViewLib::prettyUnc(mu,sig)

  bs = nptest::np.boot(
    x = Z, statistic = sd,
    R = nBoot, level = 0.95, method = "bca",
    parallel = TRUE, cl = cl)
  mus  = bs$t0
  sigs = bs$se
  sigZ[i] = ErrViewLib::prettyUnc(mus,sigs)
  biasZ[i] = abs(mu/mus)

  ## Fit by Stutent's t to get effective df
  fit.t<-fitdistrplus::fitdist(
    Z, "t_ls", start = list(df=20,mu=mean(Z),sigma=sd(Z)),
    keepdata = FALSE)
  pars = summary(fit.t)$estimate
  ZdistPars[i,] = summary(fit.t)$estimate
}
stopCluster(cl)

# Presentation table
df = cbind(
  iD = 1:length(setList),
  data = setList,
  M = size,
  Reference = NA
)
sink(file =  file.path(tabDir,'tabDataPres.tex'))
print(knitr::kable(df, 'latex'))
sink()

df = cbind(
  iD = 1:length(setList),
  mean = meanE,
  sig  = sigE,
  bias = round(100*biasE),
  nu_eff = round(EdistPars[,1],1)
)
sink(file =  file.path(tabDir,'tabDataPropE.tex'))
print(knitr::kable(df, 'latex'))
sink()

df = cbind(
  iD = 1:length(setList),
  mean = meanZ,
  sig  = sigZ,
  bias = round(100*biasZ),
  nu_eff = round(ZdistPars[,1],1)
)
sink(file =  file.path(tabDir,'tabDataPropZ.tex'))
print(knitr::kable(df, 'latex'))
sink()


png(
  file = file.path(figDir, paste0('fig_dist_E.png')),
  width  = 2*gPars$reso,
  height = 2*gPars$reso
)
par(
  mfrow = c(3, 3),
  mar = c(3,2,2,0), #gPars$mar,
  mgp = gPars$mgp,
  pty = 's',
  tcl = gPars$tcl,
  cex = 1 * gPars$cex,
  cex.main = 1,
  lwd = gPars$lwd
)
for(i in seq_along(setList)) {
  D2 = dataList[[paste0(setList[i],'_cal')]]
  E  = D2$E
  uE = D2$uE
  xlim = sd(E)*c(-3,3)
  sel = abs(E) < xlim[2]

  h1 = hist(E[sel], nclass = 25, plot = FALSE)
  ylim = c(0,1.35*max(h1$density))
  hist(
    E[sel], freq = FALSE, col = gPars$cols_tr[1],
    xlim = xlim, main = paste0('Set ',i),
    xlab = paste0('Error ',D2$unit), ylim = ylim
  )
  abline(v=0, lty = 2, col = 'gray25', lwd = 2* gPars$lwd)
  curve(
    dnorm(x,mean(E),sd(E)),
    from = xlim[1], to = xlim[2], lwd = 3*gPars$lwd,
    n= 1000, col = gPars$cols[2], add=TRUE)

  curve(
    dt_ls(x,
          df = EdistPars[i,"df"],
          mu = EdistPars[i,"mu"],
          sigma = EdistPars[i,"sigma"]),
    from = xlim[1], to = xlim[2], lwd = 3*gPars$lwd,
    n= 1000, col = gPars$cols[6], add=TRUE)

  if(i==1)
    legend(
      'topright', bty = 'n', cex=0.8,
      legend = c('Data','Normal', 'Student'),
      col = c(gPars$cols[c(1,2,6)]),
      lty = c(0,1,1), lwd = 3 * gPars$lwd,
      pch = c(22,NA,NA), pt.bg = 'gray',
      pt.lwd = 2, pt.cex = 1.5
    )

}
dev.off()


png(
  file = file.path(figDir, paste0('fig_dist_Z.png')),
  width  = 2*gPars$reso,
  height = 2*gPars$reso
)
par(
  mfrow = c(3, 3),
  mar = c(3,2,2,0), #gPars$mar,
  mgp = gPars$mgp,
  pty = 's',
  tcl = gPars$tcl,
  cex = 1 * gPars$cex,
  cex.main = 1,
  lwd = gPars$lwd
)
for(i in seq_along(setList)) {
  D2 = dataList[[paste0(setList[i],'_cal')]]
  E  = D2$E
  uE = D2$uE
  Z  = E / uE
  xlim = sd(Z) * c(-3,3)
  sel = abs(Z)/sd(Z) < xlim[2]
  h = hist(Z[sel], plot = FALSE)
  hist(
    Z[sel], freq = FALSE, col = gPars$cols_tr[1],
    xlim = xlim, main = paste0('Set ',i), #nclass = 21,
    xlab = 'Z', ylim = 1.35 * c(0,max(h$density))
  )
  abline(v=0, lty = 2, col = 'gray25', lwd = 2* gPars$lwd)
  # Normal fit
  curve(
    dnorm(x,mean(Z),sd(Z)),
    from = xlim[1], to = xlim[2], lwd = 3*gPars$lwd,
    n= 1000, col = gPars$cols[2], add=TRUE)
  # t_ls fit

  curve(
    dt_ls(x,
          df = ZdistPars[i,"df"],
          mu = ZdistPars[i,"mu"],
          sigma = ZdistPars[i,"sigma"]),
    from = xlim[1], to = xlim[2], lwd = 3*gPars$lwd,
    n= 1000, col = gPars$cols[6], add=TRUE)

  if(i==1)
    legend(
      'topright', bty = 'n', cex=0.8,
      legend = c('Data','Normal', 'Student'),
      col = c(gPars$cols[c(1,2,6)]),
      lty = c(0,1,1), lwd = 3 * gPars$lwd,
      pch = c(22,NA,NA), pt.bg = 'gray',
      pt.lwd = 2, pt.cex = 1.5
    )

}
dev.off()



png(
  file = file.path(figDir, paste0('fig_dist_Esim.png')),
  width  = 2*gPars$reso,
  height = 2*gPars$reso
)
par(
  mfrow = c(3, 3),
  mar = c(3,2,2,0), #gPars$mar,
  mgp = gPars$mgp,
  pty = 's',
  tcl = gPars$tcl,
  cex = 1 * gPars$cex,
  cex.main = 1,
  lwd = gPars$lwd
)
for(i in seq_along(setList)) {
  D2 = dataList[[paste0(setList[i],'_cal')]]
  E  = D2$E; mu = mean(E)
  uE = D2$uE
  xlim = sd(E)*c(-3,3)
  sel = abs(E) < xlim[2]

  h1 = hist(E[sel], nclass = 25, plot = FALSE)
  ylim = c(0,1.35*max(h1$density))
  hist(
    E[sel], freq = FALSE, nclass = 25, col = gPars$cols_tr[1],
    xlim = xlim, main = paste0('Set ',i),
    xlab = paste0('Error ',D2$unit),
    ylim = ylim
  )
  abline(v=0, lty = 2, col = 'gray25', lwd = 2* gPars$lwd)

  df = ZdistPars[i,"df"]
  if(df > 2) {
    sample = as.vector(
      replicate(
        100,
        uE * rt_ls(length(uE), df = df, mu = 0, sigma = 1)/
          sqrt(df/(df-2))
      )
    )
    sel = abs(sample) < xlim[2]
    h2 = hist(sample[sel], nclass = 55, plot = FALSE)
    lines(h2$mids, h2$density, lwd = 3* gPars$lwd, col = gPars$cols[6])
  }

  sample = as.vector(replicate(100, uE * rnorm(length(uE))))
  sel = abs(sample) < xlim[2]
  h3 = hist(sample[sel], nclass = 55, plot = FALSE)
  lines(h3$mids, h3$density, lwd = 3* gPars$lwd, col = gPars$cols[2])

  if(i==1)
    legend(
      'topright', bty = 'n', cex=0.8,
      legend = c('Data','Normal', 'Student'),
      col = c(gPars$cols[c(1,2,6)]),
      lty = c(0,1,1), lwd = 3 * gPars$lwd,
      pch = c(22,NA,NA), pt.bg = 'gray',
      pt.lwd = 2, pt.cex = 1.5
    )

}
dev.off()


library(MCMCpack)
png(
  file = file.path(figDir, paste0('fig_dist_uE.png')),
  width  = 2*gPars$reso,
  height = 2*gPars$reso
)
par(
  mfrow = c(3, 3),
  mar = c(3,2,2,0), #gPars$mar,
  mgp = gPars$mgp,
  pty = 's',
  tcl = gPars$tcl,
  cex = 1 * gPars$cex,
  cex.main = 1,
  lwd = gPars$lwd
)

parstab = matrix(NA, nrow = length(setList), ncol = 2)
colnames(parstab) = c("shape","scale")

for(i in seq_along(setList)) {
  D2 = dataList[[paste0(setList[i],'_cal')]]
  uE = D2$uE

  fit.t<-fitdistrplus::fitdist(
    1/uE^2, "gamma", method = "mle",
    start = list(shape = 10, scale = 1),
    keepdata = FALSE)
  pars = summary(fit.t)$estimate
  shape = pars["shape"]
  scale = pars["scale"]

  sel = uE^2 < quantile(uE^2,probs = 0.99)
  X = uE[sel]^2
  h1 = hist(X, nclass = 33, plot = FALSE)
  ylim = c(0,1.35*max(h1$density))
  hist(
    X, freq = FALSE, col = gPars$cols_tr[1],
    main = paste0('Set ',i), nclass = 33,
    xlab = 'uE^2'
  )
  curve(
    MCMCpack::dinvgamma(x,shape = shape, scale = 1/scale),
    from = min(X), to = max(X), lwd = 3*gPars$lwd,
    n= 1000, col = gPars$cols[6], add=TRUE)

  parstab[i,]=c(shape, 1/scale)

}

dev.off()

sink(file =  file.path(tabDir,'tabParsFituE.tex'))
print(knitr::kable(signif(parstab,3), 'latex'))
sink()
