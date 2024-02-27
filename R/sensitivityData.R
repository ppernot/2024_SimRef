set.seed(123)

nMC  = 10000
nBin = 50

# Get stats names
intrvt = ErrViewLib::genIntervals(1:10, nBin)
statst = names(calScoresBS2(1:10,cbind(1:10,1:10),intrvt))

# Real datasets ####

## Sensitivity to D ####
smct    = matrix(NA, nrow = nMC, ncol = length(statst))
scorest = uscorest =
  matrix(NA, nrow = length(setList), ncol = length(statst))

for(i in seq_along(setList)) {
  D2 = dataList[[paste0(setList[i],'_cal')]]
  print(setList[i])
  uEt = D2$uE
  Mt  = length(uEt)

  intrvt = ErrViewLib::genIntervals(1:Mt, nBin)

  for(kt in 1:nMC) {
    Ept = uEt * rnorm(Mt)
    smct[kt,] = calScoresBS2(1:Mt, cbind(Ept,uEt), intrvt)
  }
  scorest[i,]  = apply(smct, 2, mean, na.rm = TRUE)
  uscorest[i,] = apply(smct, 2, sd, na.rm = TRUE) / sqrt(nMC)
}


save(nBin, setList, statst, scorest, uscorest,
     file = 'sensData_s.Rda')
