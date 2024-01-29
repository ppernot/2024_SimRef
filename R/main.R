figDir = '../Figs'
tabDir = '../Tabs'
library(ErrViewLib)
library(CHNOSZ)
library(rcdk); sp <- get.smiles.parser()
library(fitdistrplus)
library(nptest)
library(parallel)
gPars = ErrViewLib::setgPars(type = 'publish')
scalePoints = 0.2

# Load functions
source('functions.R')

# Load datasets
newData = FALSE
if(newData) {
  source('getData.R')
} else {
  load(file = 'data.Rda')
}


# Results ####
nMC   = 10000
nBoot =  5000

## Data summary and properties ####
source("./datasetsProperties.R")

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
  height = 2 * gPars$reso
)
par(
  mfrow = c(3, 1),
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
            bty = 'n'
          ),
          xpd = FALSE,
          col = gPars$cols[1:length(methods)],
          main = stats[stat])
  abline(h= -5:5, col = 'gray50', lty = 3)
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

## Appendix C ####
doCalc = FALSE
if(doCalc) {
  # Takes 1/2 day to run...
  source("testValidRCE.R")
} else {
  load('testValidRCE.Rda')
}

pro = cilo = ciup = matrix(NA,nrow=length(dfList),ncol = 2)
for(j in seq_along(dfList)) {
  nu = dfList[j]
  zm = zList[[paste0(nu)]]
  tm = tList[[paste0(nu)]]
  print(colSums(zm == tm))
  pro[j,]  = colMeans(tm)
  success  = colMeans(tm) * nTry
  trials   = c(nTry, nTry)
  ci       = DescTools::BinomCI(success, trials, method = "wilsoncc")
  cilo[j,] = ci[,2]
  ciup[j,] = ci[,3]
}

png(
  file = file.path(figDir, paste0('fig_validRCE.png')),
  width  = gPars$reso,
  height = gPars$reso
)
par(
  mfrow = c(1, 1),
  mar = gPars$mar,
  mgp = gPars$mgp,
  pty = 'm',
  tcl = gPars$tcl,
  cex = gPars$cex,
  cex.main = 1,
  lwd = gPars$lwd
)
matplot(dfList, pro, type = 'b',
        pch=16:17, lty =1, col=gPars$cols[1:2],
        log = 'x',
        xlab = expression(nu),
        ylab = expression(p[val]), ylim = c(0.70,1.0))
grid(equilogs = FALSE)
abline(h=0.95,lty=2)
for(i in 1:2)
  segments(dfList,cilo[,i],dfList,ciup[,i],col=i)
box()
legend(
  'bottomright', bty = 'n',
  legend = stats,
  lty = 1,
  pch = 16:17,
  col = gPars$cols[1:2]
)
dev.off()




