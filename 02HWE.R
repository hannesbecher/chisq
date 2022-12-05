# chi-suared and HWE testing

# simulate genotypes

simGts <- function(p=0.1, nInd=100, ploidy=2){
  matrix(rbinom(n=nInd*ploidy, size=1, prob = p), nInd, ploidy)
}
simGts(.1)

# assuming biallelic of 0 and 1
getGtFreqs <- function(gtMat){
  rs <- rowSums(gtMat)
  c(aa=sum(rs==0),
    Aa=sum(rs==1),
    AA=sum(rs==2))
}
getGtFreqs(simGts())

# make contingency tab
makeContTab <- function(gtMat){
  rs <- rowSums(gtMat)
  matrix(c(sum(rs==0),
    sum(rs==1)/2,
    sum(rs==1)/2,
    sum(rs==2)),
    2,2)
}
ct01 <- makeContTab(gts01)



makeContTab(gts02)

getChisqVal <- function(cTab){
  expect <- (rowSums(cTab)/sum(cTab)) %o% colSums(cTab)
  sum(
    (expect - cTab)^2 / expect
    )
}

getChisqVal(makeContTab(gts01))
getChisqVal(makeContTab(gts02))



# same value as returned by test function
getChisqVal(makeContTab(gts01))
chisq.test(makeContTab(gts01), correct = F, simulate.p.value = F)
chisq.test(makeContTab(gts01), correct = T, simulate.p.value = F) # the default is simulate.p.value=F
#chisq.test(makeContTab(gts01), simulate.p.value = F)
chisq.test(makeContTab(gts01), simulate.p.value = T)


# run many simulations and record values
cVals <- sapply(1:10000, function(x) getChisqVal(makeContTab(simGts(nInd=100))))
#cVals <- sapply(1:10000, function(x) getChisqVal(makeContTab(simGts(nInd=1000))))
hist(cVals, freq=F)
curve(dchisq(x, df=1), 0, 15, add=T)
abline(v=qchisq(0.95, df=1), lty=2)
table(cVals > qchisq(0.95, df=1))
# with nInd=100, there are usually less than 5% false positives
# with nInd=1000, there are usually closer to 5% false positives and the curve matches the histogram better


# run many simulations and record values, different p
#cVals <- sapply(1:10000, function(x) getChisqVal(makeContTab(simGts(p=0.4, nInd=100))))
cVals <- sapply(1:10000, function(x) getChisqVal(makeContTab(simGts(p=0.4, nInd=1000))))
hist(cVals, freq=F)
curve(dchisq(x, df=1), 0, 15, add=T)
abline(v=qchisq(0.95, df=1), lty=2)
table(cVals > qchisq(0.95, df=1))
# with p=0.4, the false positives look more like 5% for nInd=100 and 1000




