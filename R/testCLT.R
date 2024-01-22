nBoot = 5000
nMC = 1e4
Mseq = c(
  seq(  30,   100, by=  10), 150,
  seq( 200,  900, by= 100), 1500,
  seq(2000,10000, by=1000)
)

# ZMS ####
tZSU95N = utZSU95N = 
  tZSU95T = utZSU95T = 
  tZSU95T3 = utZSU95T3 = c()
for (j in seq_along(Mseq)) {
  M = Mseq[j]; cat(M,'/ ')
  ZSU95N = muN = ZSU95T = muT = ZSU95T3 = muT3 = c()
  for(i in 1:nMC) {
    
    Z  = rnorm(M)
    res = ErrViewLib::ZMSCI(Z,'stud')
    muN[i] = mean(Z^2)
    u95 = diff(res$ci)/2
    ZSU95N[i] = (muN[i] - 1) / u95
    
    Z  = rT4(M, df=6)
    res = ErrViewLib::ZMSCI(Z,'stud')
    muT[i] = mean(Z^2)
    u95 = diff(res$ci)/2
    ZSU95T[i] = (muT[i] - 1) / u95
    
    Z  = rT4(M, df=4)
    res = ErrViewLib::ZMSCI(Z,'stud')
    muT3[i] = mean(Z^2)
    u95 = diff(res$ci)/2
    ZSU95T3[i] = (muT3[i] - 1) / u95
    
  }
  tZSU95N[j] = mean(abs(ZSU95N) <= 1)
  utZSU95N[j]= sd(abs(ZSU95N)   <= 1)/sqrt(nMC)

  tZSU95T[j] = mean(abs(ZSU95T) <= 1)
  utZSU95T[j]= sd(abs(ZSU95T)   <= 1)/sqrt(nMC)

  tZSU95T3[j] = mean(abs(ZSU95T3) <= 1)
  utZSU95T3[j]= sd(abs(ZSU95T3)   <= 1)/sqrt(nMC)
}

png(
  file = file.path(figDir, 'figLUP.png'),
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

x = Mseq
y = tZSU95N
uy = utZSU95N
plot(
  x, y, type = 'b',
  pch = 16, col = gPars$cols[1],
  log = 'x', xlab = 'M',
  ylim = c(0.80, 0.96), ylab = 'Coverage probability',
  main = ''
)
segments(x, y-2*uy, x, y+2*uy, col = gPars$cols[1])
grid(equilogs = FALSE)
abline(h = 0.95, lty = 2)

y = tZSU95T
uy = utZSU95T
points(x, y, type = 'b', pch = 17, col = gPars$cols[2])
segments(x, y-2*uy, x, y+2*uy, col = gPars$cols[2])

y = tZSU95T3
uy = utZSU95T3
points(x, y, type = 'b', pch = 18, col = gPars$cols[3])
segments(x, y-2*uy, x, y+2*uy, col = gPars$cols[3])

# y = tZSU95T3BS
# uy = utZSU95T3BS
# points(x, y, type = 'b', pch = 19, col = gPars$cols[4])
# segments(x, y-2*uy, x, y+2*uy, col = gPars$cols[4])

box()
legend(
  'bottomright', bty = 'n',
  legend = c('Normal','Student: nu = 6','Student: nu = 4'),
  col = gPars$cols[1:3],
  pch = 16:18,
  lty = 1
)
dev.off()
