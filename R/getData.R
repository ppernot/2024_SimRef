# Get data ####
dataList = list()
setList = c()
isl = 0


## PAL2022 ####

unit = '[eV]'
sets = c('Diffusion', 'Perovskite')[1:2]
methods = c('RF', 'LR', 'GPR_Bayesian')[1:3]
cases = c('Test_uncal', 'Test_cal')

for (method in methods) {
  for (set in sets) {
    fRoot = paste0(set, '_', method)
    isl = isl + 1
    setList[isl] = paste0(fRoot,'_Test')
    for (case in cases) {
      dataName = paste0(fRoot, '_', case)
      fName = file.path('..', 'Data', 'PAL2022',
                        paste0(fRoot, '_', case, '.csv'))
      if (!file.exists(fName))
        next
      data = read.csv(fName, header = TRUE)
      sel = data$uE > 1e-6 * sd(data$E)
      if (sum(!sel) > 0) {
        print(c(case, ' non-pos unc:', sum(!sel)))
        data = data[sel,]
      }
      dataList[[dataName]] =
        list(
          unit = unit,
          E    = data$E,
          uE   = data$uE
        )
    }
  }
}
## BUS2022 ####

unit = '[eV]'
sets = c('qm9_qm9_E_uncalibrated_test',
         'qm9_qm9_E_calibrated_isotonic_test')
aliases = c('QM9_E_Test_uncal', 'QM9_E_Test_cal')

isl = isl + 1
setList[isl] = 'QM9_E_Test'
for (i in seq_along(sets)) {
  fName = file.path('..', 'Data', 'BUS2022', paste0(sets[i],'.csv'))
  D = read.table(fName,
                 sep = ',',
                 header = TRUE,
                 check.names = FALSE,
                 stringsAsFactors = FALSE
  )

  S  = D[, "formula"]
  R  = D[, "E_"]
  C  = D[, "prediction"]
  uC = D[, "uncertainty_total"] ^ 0.5
  E  = R - C
  uE = uC

  dataList[[aliases[i]]] =
    list(
      unit = unit,
      E    = E,
      uE   = uE
    )
}

## RAS2023 ####
unit = ''
sets = c(
  'logP_10k_a_LS-GCN_test',
  'logP_150k_LS-GCN_test'
)

for (i in seq_along(sets)) {
  isl = isl + 1
  setList[isl] = sets[i]
  fName = file.path('..', 'Data', 'RAS2023', paste0(sets[i],'.csv'))
  D = read.table(fName,
                 sep = ',',
                 colClasses = c('integer','character','numeric',
                                'numeric', 'numeric', 'numeric'),
                 header = TRUE,
                 comment.char = "",
                 check.names = FALSE,
                 stringsAsFactors = FALSE
  )

  R  = D[, "logP"]
  uC = D[, "uq"]
  C  = D[, "y_pred"]
  E  = R - C
  uE = uC

  dataList[[paste0(setList[isl],'_cal')]] =
    list(
      unit = unit,
      E    = E,
      uE   = uE
    )
}

save(setList, dataList, file = 'data.Rda')
