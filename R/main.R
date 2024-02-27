setwd("~/Bureau/2024_SimRef/R")

figDir = '../Figs'
tabDir = '../Tabs'

library(fitdistrplus)
library(nptest)
library(parallel)
library(MCMCpack)
library(ErrViewLib)
gPars = ErrViewLib::setgPars(type = 'publish')

# Load functions ####
source('functions.R')

# Load datasets ####
source('getData.R')


# Results ####
nMC   = 10000
nBoot =  5000

## Data summary and properties ####
source("./datasetsProperties.R")

## Sensitivity ####

### Sensitivity of ref to dataset

source("sensitivityData.R")

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


### Sensitivity of ref to generative dist. D

calcSens = FALSE
if(calcSens) {
  # Runtime about 1 hour
  source("sensitivity.R")
} else {
  load("./sensitivity.Rda")
}

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

## Scores ####
calcScores = FALSE
if(calcScores) {
  source("calcZetaScores.R")
} else {
  load("./zetaScores.Rda")
}

### Results table and Fig ####
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
if(FALSE)
  source("testCLT.R")






