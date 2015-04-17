##additional calls in RelocaTEi compared to RelocaTE
relocate <- read.table("RIL275_RelocaTE.sofia.sorted.table", skip=1)
relocatei <- read.table("RIL275_RelocaTEi.summary.table", skip=1)
temp <- read.table("RIL275_TEMP.summary.table", skip=1)

more <- relocatei[,8] > relocate[,3]
addition <- relocatei[,8] - relocate[,3]
confident <- summary(addition[more])
write('Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA', file='AdditionRelocaTEi.summary')
write('>Confident', file='AdditionRelocaTEi.summary', append=TRUE)
write(confident, file='AdditionRelocaTEi.summary', append=TRUE)

more <- relocatei[,8] > relocate[,3]
addition <- relocatei[,7] - relocate[,3]
all <- summary(addition[more])
write('>ALL', file='AdditionRelocaTEi.summary', append=TRUE)
write(all, file='AdditionRelocaTEi.summary', append=TRUE)
