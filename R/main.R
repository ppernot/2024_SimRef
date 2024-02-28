setwd("~/Bureau/2024_SimRef/R")

figDir = '../Figs'
tabDir = '../Tabs'
tmpDir = '../Tmp'

library(fitdistrplus)
library(nptest)
library(parallel)
library(MCMCpack)
library(ErrViewLib)
gPars = ErrViewLib::setgPars(type = 'publish')

# Redo lengthy calculations ?
calcSens1  = FALSE
calcSens2  = FALSE
calcScores = FALSE
calcAppB   = FALSE

# Load functions ####
source('functions.R')

# Load datasets ####
source('getData.R')


# Results ####
nMC   = 10000
nBoot =  5000

## Data summary and properties ####
source("./datasetsProperties.R")

### Tables ####
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

### Fig. 2 ####

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

gi = bgm = c()

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

  # sel = uE^2 < quantile(uE^2,probs = 0.99)
  # X = uE[sel]^2
  X = log(uE^2)
  lMeaV = log(mean(uE^2))
  lMedV = log(median(uE^2))
  h1 = hist(X, nclass = 33, plot = FALSE)
  ylim = c(0,1.1*max(h1$density))
  hist(
    X, freq = FALSE, col = NULL,
    border = gPars$cols_tr2[1],
    main = paste0('Set ',i), nclass = 33,
    ylim = ylim, yaxs = 'i',
    xlab = 'log(uE^2)'
  )
  curve(
    MCMCpack::dinvgamma(exp(x),
                        shape = shape,
                        scale = 1/scale)*exp(x),
    from = min(X), to = max(X), lwd = 2*gPars$lwd,
    n= 1000, col = gPars$cols[6], add=TRUE)

  abline(v = lMeaV, lwd = 2*gPars$lwd,
         col = gPars$cols[2], lty = 1)
  abline(v = lMedV, lwd = 2*gPars$lwd,
         col = gPars$cols[3], lty = 1)
  if(i==3)
    legend(
      'topright', bty = 'n',
      legend = c('IG fit','MV', 'MedV'),
      pch = NA, lwd = 2*gPars$lwd,
      lty = c(1,1,1),
      col = gPars$cols[c(6,2,3)]
    )
  box()


  parstab[i,]=c(shape, 1/scale)

  gi[i] = ErrViewLib::gimc(uE)
  bgm[i] = ErrViewLib::skewgm(uE)
}

dev.off()

parstab = cbind(gi, bgm, parstab)
sink(file =  file.path(tabDir,'tabParsFituE.tex'))
print(knitr::kable(signif(parstab,3), 'latex'))
sink()

### Fig. 3 ####

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
    E[sel], freq = FALSE, col = NULL, nclass=25,
    border = gPars$cols_tr2[1], yaxs = 'i',
    xlim = xlim, main = paste0('Set ',i),
    xlab = paste0('Error ',D2$unit), ylim = ylim
  )
  abline(v=0, lty = 2, col = 'gray25', lwd = 2* gPars$lwd)
  curve(
    dnorm(x,mean(E),sd(E)),
    from = xlim[1], to = xlim[2], lwd = 2*gPars$lwd,
    n= 1000, col = gPars$cols[2], add=TRUE)

  curve(
    dt_ls(x,
          df = EdistPars[i,"df"],
          mu = EdistPars[i,"mu"],
          sigma = EdistPars[i,"sigma"]),
    from = xlim[1], to = xlim[2], lwd = 2*gPars$lwd,
    n= 1000, col = gPars$cols[6], add=TRUE)

  if(i==1)
    legend(
      'topright', bty = 'n', cex=0.75,
      legend = c('Data','Normal', 'Student'),
      col = gPars$cols[c(1,2,6)],
      lty = c(0,1,1), lwd = 2 * gPars$lwd,
      pch = c(22,NA,NA), pt.bg = 'white',
      pt.lwd = 2, pt.cex = 1.5
    )
  box()
}
dev.off()

### Fig. 4 ####

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
  h = hist(Z[sel], nclass = 25, plot = FALSE)
  hist(
    Z[sel], freq = FALSE, col = NULL, nclass = 25,
    border = gPars$cols_tr2[1],
    xlim = xlim, main = paste0('Set ',i), yaxs = 'i',
    xlab = 'Z', ylim = 1.35 * c(0,max(h$density))
  )
  abline(v=0, lty = 2, col = 'gray25', lwd = 2* gPars$lwd)
  # Normal fit
  curve(
    dnorm(x,mean(Z),sd(Z)),
    from = xlim[1], to = xlim[2], lwd = 2*gPars$lwd,
    n= 1000, col = gPars$cols[2], add=TRUE)
  # Student's fit
  curve(
    dt_ls(x,
          df = ZdistPars[i,"df"],
          mu = ZdistPars[i,"mu"],
          sigma = ZdistPars[i,"sigma"]),
    from = xlim[1], to = xlim[2], lwd = 2*gPars$lwd,
    n= 1000, col = gPars$cols[6], add=TRUE)

  if(i==1)
    legend(
      'topright', bty = 'n', cex=0.75,
      legend = c('Data','Normal', 'Student'),
      col = c(gPars$cols[c(1,2,6)]),
      lty = c(0,1,1), lwd = 3 * gPars$lwd,
      pch = c(22,NA,NA), pt.bg = 'white',
      pt.lwd = 2, pt.cex = 1.5
    )
  box()
}
dev.off()

### Fig. 5 ####

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
    E[sel], freq = FALSE, nclass = 25, col = NULL,
    border = gPars$cols_tr2[1],
    xlim = xlim, main = paste0('Set ',i),
    xlab = paste0('Error ',D2$unit),
    ylim = ylim, yaxs = 'i'
  )
  abline(v=0, lty = 2, col = 'gray25', lwd = 2* gPars$lwd)

  df = ZdistPars[i,"df"]
  if(df <= 2)
    df = 2.1
  sample = as.vector(
    replicate(
      100,
      uE * rt_ls(length(uE), df = df, mu = 0, sigma = 1)/
        sqrt(df/(df-2))
    )
  )
  sel = abs(sample) < xlim[2]
  h2 = hist(sample[sel], nclass = 55, plot = FALSE)
  lines(h2$mids, h2$density, lwd = 2* gPars$lwd,
        col = gPars$cols[6])

  sample = as.vector(replicate(100, uE * rnorm(length(uE))))
  sel = abs(sample) < xlim[2]
  h3 = hist(sample[sel], nclass = 55, plot = FALSE)
  lines(h3$mids, h3$density, lwd = 2* gPars$lwd,
        col = gPars$cols[2])

  if(i==1)
    legend(
      'topright', bty = 'n', cex=0.8,
      legend = c('Data','Normal', 'Student'),
      col = c(gPars$cols[c(1,2,6)]),
      lty = c(0,1,1), lwd = 2 * gPars$lwd,
      pch = c(22,NA,NA), pt.bg = 'white',
      pt.lwd = 2, pt.cex = 1.5
    )
  box()
}
dev.off()



## Sensitivity of theta_ref to uE ####
if(calcSens1) {
  source("sensitivityData.R")
} else {
  load(file = file.path(tmpDir,"sensData.Rda"))
}

### Fig. 6 ####
png(
  file = file.path(figDir, paste0('fig_sensData_s.png')),
  width  = 2*gPars$reso,
  height = 2.75*gPars$reso
)
par(
  mfrow = c(3,2),
  mar = c(4,4,1.5,1),
  mgp = gPars$mgp,
  pty = 'm',
  tcl = gPars$tcl,
  cex = 1.25*gPars$cex,
  cex.main = 1,
  lwd = gPars$lwd
)
nSets = length(setList)
label = 0
for(istat in seq_along(statst)) {
  sc  = scorest[,istat]
  usc = uscorest[,istat]
  ylim = range(c(sc-2*usc,sc+2*usc))
  plot( 1:nSets, sc, log = '', pch = 16,
        type = 'p', col = gPars$cols[5],
        xlab = "Set #", xaxt = 'n',
        ylim = ylim,
        ylab = expression(tilde(theta)[list(nu,ref)]),
        main = statst[istat],
        panel.first = grid()
  )
  axis(1, at = 1:nSets, labels = 1:nSets, gap.axis = 0)
  # grid()
  if(statst[istat]== "ZMS")
    abline(h = 1, lwd = gPars$lwd, lty = 2, col = gPars$cols[2] )
  for(is in 1:nSets)
    arrows(is,sc[is]-2*usc[is],is,sc[is]+2*usc[is],
           angle = 90, lwd = 2* gPars$lwd, col = gPars$cols[5],
           code = 3)
  box()
  label = label + 1
  mtext(
    text = paste0('(', letters[label], ')'),
    side = 3,
    adj = 1,
    cex = par()$cex,
    line = 0.3)

  if(statst[istat] == 'ENCE' | statst[istat] == 'ZMSE') {
    sc  = scorest[,istat] * sqrt(size)
    usc = uscorest[,istat] * sqrt(size)
    ylim = range(c(sc-2*usc,sc+2*usc))
    plot( 1:nSets, sc, log = '', pch = 16,
          type = 'p', col = gPars$cols[5],
          lwd = 2* gPars$lwd, lty = 1,
          xlab = "Set #", xaxt = 'n',
          ylim = ylim,
          ylab = expression(tilde(theta)[list(nu,ref)]),
          main = paste0(statst[istat],'* M^1/2'),
          panel.first = grid()
    )
    axis(1, at = 1:nSets, labels = 1:nSets, gap.axis = 0)
    for(is in 1:nSets)
      arrows(is,sc[is]-2*usc[is],is,sc[is]+2*usc[is],
             angle = 90, lwd = 2* gPars$lwd, col = gPars$cols[5],
             code = 3)
    box()
    label = label + 1
    mtext(
      text = paste0('(', letters[label], ')'),
      side = 3,
      adj = 1,
      cex = par()$cex,
      line = 0.3)

  }
}
dev.off()


## Sensitivity of theta_ref to D ####
if(calcSens2) {
  # Runtime about 1 hour
  source("sensitivity.R")
} else {
  load(file = file.path(tmpDir,"sensitivity.Rda"))
}

### Fig. 7 ####
png(
  file = file.path(figDir, paste0('fig_sensitivity78.png')),
  width  = 2*gPars$reso,
  height = 2*gPars$reso
)
par(
  mfrow = c(2,2),
  mar = c(4,4,1.5,1),
  mgp = gPars$mgp,
  pty = 'm',
  tcl = gPars$tcl,
  cex = 1.25*gPars$cex,
  cex.main = 1,
  lwd = gPars$lwd
)
label = 0
for(stat in stats) {
  sc  = t(t(scores[,,stat]))
  usc = t(t(uscores[,,stat]))
  ylim = range(c(sc-2*usc,sc+2*usc))
  matplot( nuSeq, sc, log = 'x', pch = 1,
           type = 'b', col = gPars$cols[c(5,2)],
           lwd = 2* gPars$lwd, lty = 1,
           xlab = expression(nu),
           ylim = ylim,
           ylab = expression(tilde(theta)[list(nu,ref)]),
           main = stat,
           panel.first = grid()
  )
  abline(h = scores[length(nuSeq),,stat],
         lwd = gPars$lwd, lty = 2, col = gPars$cols[c(5,2)] )
  for(inu in seq_along(nuSeq)) {
    nu = nuSeq[inu]
    arrows(nu,sc[inu,1]-2*usc[inu,1],nu,sc[inu,1]+2*usc[inu,1],
           angle = 90, lwd = 2* gPars$lwd, col = gPars$cols[5],
           code = 3)
    arrows(nu,sc[inu,2]-2*usc[inu,2],nu,sc[inu,2]+2*usc[inu,2],
           angle = 90, lwd = 2* gPars$lwd, col = gPars$cols[2],
           code = 3)
  }
  if(stat == "ZMS")
    legend(
      "topright", bty = "n",
      legend = c("Set 7", "Set 8"),
      col = gPars$cols[c(5,2)], lty = 1, pch = 1
    )
  box()
  label = label + 1
  mtext(
    text = paste0('(', letters[label], ')'),
    side = 3,
    adj = 1,
    cex = par()$cex,
    line = 0.3)
}
dev.off()

## Sensitivity of zeta-scores to D ####
if(calcScores) {
  source("calcZetaScores.R")
} else {
  load(file = file.path(tmpDir,"zetaScores.Rda"))
}

### Tables ####
for(i in seq_along(stats)) {
  df = data.frame(
    set      = 1:length(setList),
    stat     = signif(scores[,i],3),
    bias     = signif(bias[,i],2),
    CI       = ciScores[,i],
    zBS      = round(zmatBS[,i],2),

    muRefN   = signif(muSimN[,i],3),
    seRefN   = signif(seSimN[,i],2),
    zSimN    = round(zmatSimN[,i],2),
    ci2N     = ci2N[,i],
    zSim2N   = round(zmatSim2N[,i],2),

    muRefT   = signif(muSimT[,i],3),
    seRefT   = signif(seSimT[,i],2),
    zSimT    = round(zmatSimT[,i],2),
    ci2T     = ci2T[,i],
    zSim2T   = round(zmatSim2T[,i],2)
  )
  # print(knitr::kable(df))
  sink(file = file.path(tabDir,paste0('tab',stats[i],'.tex')))
  print(knitr::kable(df, 'latex'))
  sink()
}

### Fig. 8 ####
png(
  file = file.path(figDir, 'fig_scores_linear1.png'),
  width  = 1.5 * gPars$reso,
  height = 1.5 * gPars$reso
)
par(
  mfrow = c(2, 1),
  mar = gPars$mar,
  mgp = gPars$mgp,
  pty = 'm',
  tcl = gPars$tcl,
  cex = 0.85 * gPars$cex,
  lwd = gPars$lwd
)
nS = length(setList)
for(stat in 1:2) {
  df = c()
  for (meth in methods) {
    zm = get(paste0('zmat',meth))[,stat]
    if (meth == methods[1])
      df = zm
    else
      df = cbind(df,zm)
  }
  colnames(df) = methods
  rownames(df) = paste0('Set',1:nS)
  ylim = 1.2*range(df, na.rm = TRUE)
  barplot(t(df), beside=TRUE,
          ylim = ylim, ylab = expression(zeta - score),
          legend.text = stat == 1,
          args.legend = list(
            x = 'topright',
            ncol = 2,
            bty = 'n'
          ),
          xpd = FALSE,
          col = gPars$cols[1:length(methods)],
          main = stats[stat])
  abline(h= -15:15, col = 'gray50', lty = 3)
  abline(h=c(-1,0,1), col = gPars$cols[1], lty = c(2,1,2))
  barplot(t(df), beside=TRUE,
          xpd = FALSE,
          col = gPars$cols[1:length(methods)],
          add = TRUE)

  box()
}
dev.off()

### Fig. 9 ####
png(
  file = file.path(figDir, 'fig_scores_linear2.png'),
  width  = 1.5 * gPars$reso,
  height = 1.5 * gPars$reso
)
par(
  mfrow = c(2, 1),
  mar = gPars$mar,
  mgp = gPars$mgp,
  pty = 'm',
  tcl = gPars$tcl,
  cex = 0.85 * gPars$cex,
  lwd = gPars$lwd
)
nS = length(setList)
for(stat in 3:4) {
  df = c()
  for (meth in methods) {
    zm = get(paste0('zmat',meth))[,stat]
    if (meth == methods[1])
      df = zm
    else
      df = cbind(df,zm)
  }
  colnames(df) = methods
  rownames(df) = paste0('Set',1:nS)
  ylim = c(0,min(8,1.2*max(df, na.rm = TRUE)))
  barplot(t(df), beside=TRUE,
          ylim = ylim, ylab = expression(zeta - score),
          legend.text = stat == 1,
          args.legend = list(
            x = 'topright',
            ncol = 2,
            bty = 'n'
          ),
          xpd = FALSE,
          col = gPars$cols[1:length(methods)],
          main = stats[stat])
  abline(h= -15:15, col = 'gray50', lty = 3)
  abline(h=c(-1,0,1), col = gPars$cols[1], lty = c(2,1,2))
  barplot(t(df), beside=TRUE,
          xpd = FALSE,
          col = gPars$cols[1:length(methods)],
          add = TRUE)

  box()
}
dev.off()

# Appendices ####
## Appendix B ####
if(calcAppB)
  source('testENCE.R')

### Fig. 10 ####
png(
  file = file.path(figDir, paste0('fig_ENCE_ZMSE.png')),
  width  = 2*gPars$reso,
  height = 2*gPars$reso
)
par(
  mfrow = c(2,2),
  mar = c(4,4,1,1),
  mgp = gPars$mgp,
  pty = 's',
  tcl = gPars$tcl,
  cex = 1.0*gPars$cex,
  cex.main = 1,
  lwd = gPars$lwd
)
lwd = 1.5 *gPars$lwd
pch = c(0,1,2,3,6)

load(file = file.path(tmpDir,"testENCE.Rda"))
ylim = range(c(scorest[,,,1]-2*uscorest[,,,1],
               scorest[,,,1]+2*uscorest[,,,1]))
for(k in seq_along(nuSeq)) {
  sc  = scorest[,k,,1]
  usc = uscorest[,k,,1]
  if(k==1) {
    matplot(mSeq, sc, type = 'b', log = 'xy',
            pch = pch,
            lty = 1, col = gPars$cols[k],
            lwd = lwd,
            xlab = 'Set size (M)',
            ylim = ylim,
            ylab = expression(tilde(theta)[list(nu,ref)]),
            main = 'ENCE')
    grid(equilogs = FALSE)
  } else {
    matlines(mSeq, sc, type = 'b',
             lwd = lwd,
             pch = pch, lty = 1, col = gPars$cols[k])
  }
}
legend(
  'topright', bty = 'n', cex = 0.85, ncol = 2,
  title = expression(nu/2),
  legend = nuSeq,
  lty = 1, pch = 20, lwd = lwd,
  col = gPars$cols
)
legend(
  'bottomleft', bty = 'n', cex = 0.85,
  title = 'N',
  legend = nBinSeq,
  lty = 0, pch = pch, lwd = lwd,
  col = gPars$cols[1]
)
box()

xt = yt = c()
for(k in seq_along(nuSeq)) {
  for(j in seq_along(nBinSeq)) {
    nBint = nBinSeq[j]
    xt = c(xt, sqrt(nBint/mSeq))
    yt = c(yt, scorest[,k,j,1])
  }
}

xlim = c(0, 1.1*max(xt))
ylim = c(0, 1.1*max(yt))
for(k in seq_along(nuSeq)) {
  for(j in seq_along(nBinSeq)) {
    nBint = nBinSeq[j]
    x = sqrt(nBint/mSeq)
    y = scorest[,k,j,1]
    if(k == 1 & j == 1) {
      plot(x, y, type = 'b',
           pch = pch[j], lty = 1, col = gPars$cols[k],
           lwd = lwd,
           xlim = xlim, xaxs = 'i',
           xlab = '(N/M)^1/2',
           ylim = ylim, yaxs = 'i',
           ylab = expression(tilde(theta)[list(nu,ref)]),
           main = 'ENCE')
      grid()
    } else {
      lines(x, y, type = 'b',
            lwd = lwd,
            pch = pch[j], lty = 1, col = gPars$cols[k])
    }
  }
}
abline(reg = lm(yt~xt), lty = 2, col = 'gray25')
box()
print(coefficients(lm(yt~xt)))

load(file = file.path(tmpDir,"testZMSE.Rda"))
ylim = range(c(scorest[,,,1]-2*uscorest[,,,1],
               scorest[,,,1]+2*uscorest[,,,1]))
for(k in seq_along(nuSeq)) {
  sc  = scorest[,k,,1]
  usc = uscorest[,k,,1]
  if(k==1) {
    matplot(mSeq, sc, type = 'b', log = 'xy',
            pch = pch,
            lty = 1, col = gPars$cols[k],
            lwd = lwd,
            xlab = 'Set size (M)',
            ylim = ylim,
            ylab = expression(tilde(theta)[list(nu,ref)]),
            main = 'ZMSE')
    grid(equilogs = FALSE)
  } else {
    matlines(mSeq, sc, type = 'b',
             lwd = lwd,
             pch = pch, lty = 1, col = gPars$cols[k])
  }
}
legend(
  'topright', bty = 'n', cex = 0.85, ncol = 2,
  title = expression(nu/2),
  legend = nuSeq,
  lty = 1, pch = 20, lwd = lwd,
  col = gPars$cols
)
legend(
  'bottomleft', bty = 'n', cex = 0.85,
  title = 'N',
  legend = nBinSeq,
  lty = 0, pch = pch, lwd = lwd,
  col = gPars$cols[1]
)
box()

xt = yt = c()
for(k in seq_along(nuSeq)) {
  for(j in seq_along(nBinSeq)) {
    nBint = nBinSeq[j]
    xt = c(xt, sqrt(nBint/mSeq))
    yt = c(yt, scorest[,k,j,1])
  }
}

xlim = c(0, 1.1*max(xt))
ylim = c(0, 1.1*max(yt))
for(k in seq_along(nuSeq)) {
  for(j in seq_along(nBinSeq)) {
    nBint = nBinSeq[j]
    x = sqrt(nBint/mSeq)
    y = scorest[,k,j,1]
    if(k == 1 & j == 1) {
      plot(x, y, type = 'b',
           pch = pch[j], lty = 1, col = gPars$cols[k],
           lwd = lwd,
           xlim = xlim, xaxs = 'i',
           xlab = '(N/M)^1/2',
           ylim = ylim, yaxs = 'i',
           ylab = expression(tilde(theta)[list(nu,ref)]),
           main = 'ZMSE')
      grid()
    } else {
      lines(x, y, type = 'b',
            lwd = lwd,
            pch = pch[j], lty = 1, col = gPars$cols[k])
    }
  }
}
abline(reg = lm(yt~xt), lty = 2, col = 'gray25')
box()
print(coefficients(lm(yt~xt)))

dev.off()

### Fig. 11 ####

png(
  file = file.path(figDir, paste0('fig_ENCE_ZMSE_T6.png')),
  width  = 2*gPars$reso,
  height = 2*gPars$reso
)
par(
  mfrow = c(2,2),
  mar = c(4,4,1,1),
  mgp = gPars$mgp,
  pty = 's',
  tcl = gPars$tcl,
  cex = 1.0*gPars$cex,
  cex.main = 1,
  lwd = gPars$lwd
)
lwd = 1.5 *gPars$lwd
pch = c(0,1,2,3,6)

load(file = file.path(tmpDir,"testENCE_T6.Rda"))
ylim = range(c(scorest[,,,1]-2*uscorest[,,,1],
               scorest[,,,1]+2*uscorest[,,,1]))
for(k in seq_along(nuSeq)) {
  sc  = scorest[,k,,1]
  usc = uscorest[,k,,1]
  if(k==1) {
    matplot(mSeq, sc, type = 'b', log = 'xy',
            pch = pch,
            lty = 1, col = gPars$cols[k],
            lwd = lwd,
            xlab = 'Set size (M)',
            ylim = ylim,
            ylab = expression(tilde(theta)[list(nu,ref)]),
            main = 'ENCE')
    grid(equilogs = FALSE)
  } else {
    matlines(mSeq, sc, type = 'b',
             lwd = lwd,
             pch = pch, lty = 1, col = gPars$cols[k])
  }
}
legend(
  'topright', bty = 'n', cex = 0.85, ncol = 2,
  title = expression(nu/2),
  legend = nuSeq,
  lty = 1, pch = 20, lwd = lwd,
  col = gPars$cols
)
legend(
  'bottomleft', bty = 'n', cex = 0.85,
  title = 'N',
  legend = nBinSeq,
  lty = 0, pch = pch, lwd = lwd,
  col = gPars$cols[1]
)
box()

xt = yt = c()
for(k in seq_along(nuSeq)) {
  for(j in seq_along(nBinSeq)) {
    nBint = nBinSeq[j]
    xt = c(xt, sqrt(nBint/mSeq))
    yt = c(yt, scorest[,k,j,1])
  }
}

xlim = c(0, 1.1*max(xt))
ylim = c(0, 1.1*max(yt))
for(k in seq_along(nuSeq)) {
  for(j in seq_along(nBinSeq)) {
    nBint = nBinSeq[j]
    x = sqrt(nBint/mSeq)
    y = scorest[,k,j,1]
    if(k == 1 & j == 1) {
      plot(x, y, type = 'b',
           pch = pch[j], lty = 1, col = gPars$cols[k],
           lwd = lwd,
           xlim = xlim, xaxs = 'i',
           xlab = '(N/M)^1/2',
           ylim = ylim, yaxs = 'i',
           ylab = expression(tilde(theta)[list(nu,ref)]),
           main = 'ENCE')
      grid()
    } else {
      lines(x, y, type = 'b',
            lwd = lwd,
            pch = pch[j], lty = 1, col = gPars$cols[k])
    }
  }
}
abline(reg = lm(yt~xt), lty = 2, col = 'gray25')
box()
print(coefficients(lm(yt~xt)))

load(file = file.path(tmpDir,"testZMSE_T6.Rda"))
ylim = range(c(scorest[,,,1]-2*uscorest[,,,1],
               scorest[,,,1]+2*uscorest[,,,1]))
for(k in seq_along(nuSeq)) {
  sc  = scorest[,k,,1]
  usc = uscorest[,k,,1]
  if(k==1) {
    matplot(mSeq, sc, type = 'b', log = 'xy',
            pch = pch,
            lty = 1, col = gPars$cols[k],
            lwd = lwd,
            xlab = 'Set size (M)',
            ylim = ylim,
            ylab = expression(tilde(theta)[list(nu,ref)]),
            main = 'ZMSE')
    grid(equilogs = FALSE)
  } else {
    matlines(mSeq, sc, type = 'b',
             lwd = lwd,
             pch = pch, lty = 1, col = gPars$cols[k])
  }
}
legend(
  'topright', bty = 'n', cex = 0.85, ncol = 2,
  title = expression(nu/2),
  legend = nuSeq,
  lty = 1, pch = 20, lwd = lwd,
  col = gPars$cols
)
legend(
  'bottomleft', bty = 'n', cex = 0.85,
  title = 'N',
  legend = nBinSeq,
  lty = 0, pch = pch, lwd = lwd,
  col = gPars$cols[1]
)
box()

xt = yt = c()
for(k in seq_along(nuSeq)) {
  for(j in seq_along(nBinSeq)) {
    nBint = nBinSeq[j]
    xt = c(xt, sqrt(nBint/mSeq))
    yt = c(yt, scorest[,k,j,1])
  }
}

xlim = c(0, 1.1*max(xt))
ylim = c(0, 1.1*max(yt))
for(k in seq_along(nuSeq)) {
  for(j in seq_along(nBinSeq)) {
    nBint = nBinSeq[j]
    x = sqrt(nBint/mSeq)
    y = scorest[,k,j,1]
    if(k == 1 & j == 1) {
      plot(x, y, type = 'b',
           pch = pch[j], lty = 1, col = gPars$cols[k],
           lwd = lwd,
           xlim = xlim, xaxs = 'i',
           xlab = '(N/M)^1/2',
           ylim = ylim, yaxs = 'i',
           ylab = expression(tilde(theta)[list(nu,ref)]),
           main = 'ZMSE')
      grid()
    } else {
      lines(x, y, type = 'b',
            lwd = lwd,
            pch = pch[j], lty = 1, col = gPars$cols[k])
    }
  }
}
abline(reg = lm(yt~xt), lty = 2, col = 'gray25')
box()
print(coefficients(lm(yt~xt)))

dev.off()

## Appendix C ####
### Fig. 12 ####
png(
  file = file.path(figDir, paste0('fig_ZMSEval.png')),
  width  = 2.5*gPars$reso,
  height = 2.5*gPars$reso
)
par(
  mfrow = c(3,3),
  mar = c(4,4,1,1),
  mgp = gPars$mgp,
  pty = 's',
  tcl = gPars$tcl,
  cex = 1.25*gPars$cex,
  cex.main = 1,
  lwd = 2*gPars$lwd
)
lwd = 1.5 *gPars$lwd
for(i in seq_along(setList)) {
  D2 = dataList[[paste0(setList[i],'_cal')]]; print(setList[i])
  io = order(D2$uE)
  uEt = D2$uE[io]
  Et  = D2$E[io]
  res = varBinSize(Et, uEt)
  Mt = length(uEt)

  x = sqrt(res$nBins / Mt)
  y = res$ZMSE
  plot(x, y,
       col = gPars$cols[5], pch = 16, cex = 1,
       xlab = '(N/M)^1/2', xlim = c(0,1.1*max(x)), xaxs = 'i',
       ylab = 'ZMSE', ylim = c(0,1.1*max(y)), yaxs = 'i',
       main = paste0('Set ',i))
  grid()

  sel = res$nBins > 20
  x   = x[sel]
  y   = y[sel]
  reg = lm(y~x)
  int = summary(reg)$coefficients[1,1]
  Uint = 2*summary(reg)$coefficients[1,2]
  abline(reg = reg, lty = 1, lwd = lwd, col = gPars$cols[2])
  segments(
    0.0005,int-Uint,
    0.0005,int+Uint,
    lwd = 4*lwd,
    lty = 1,
    lend = 2,
    col = gPars$cols[2])
  abline(a=0, b=1.14, lty = 2, lwd = lwd,
         col = gPars$cols[1])
  abline(a=0.006, b=1.577, lty = 3, lwd = lwd,
         col = gPars$cols[1])
  box()
  if(i==1)
    legend(
      'topleft', bty = 'n', cex = 0.65,
      legend = c('Data','Linear fit','NIG ref.','T6IG ref.'),
      lty = c(0,1,2,3),
      pch = c(16,NULL,NULL,NULL),
      col = gPars$cols[c(5,2,1,1)]
    )
}
dev.off()



