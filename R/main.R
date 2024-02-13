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
calcSens = FALSE
if(calcScores) {
  source("sensitivity.R")
} else {
  load("./sensitivity.Rda")
}


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
  file = file.path(figDir, 'fig_scores_linear.png'),
  width  = 1.5 * gPars$reso,
  height = 3 * gPars$reso
)
par(
  mfrow = c(4, 1),
  mar = gPars$mar,
  mgp = gPars$mgp,
  pty = 'm',
  tcl = gPars$tcl,
  cex = 0.85 * gPars$cex,
  lwd = gPars$lwd
)
nS = length(setList)
for(stat in 1:length(stats)) {
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

# Appendices ####

## Appendix B ####
if(FALSE)
  source("testCLT.R")






