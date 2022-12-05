# Testing for LD

# simulate genotypes
simGts <- function(p=0.1, nInd=100, ploidy=2){
  matrix(rbinom(n=nInd*ploidy, size=1, prob = p), nInd, ploidy)
}

simGts2loci <- function(p1=0.1, p2=0.4, nInd=100, ploidy=2){
  cbind(simGts(p1, nInd, ploidy),
        simGts(p2, nInd, ploidy)
        )
}

gt2l01 <- simGts2loci()

makeContTab2l <- function(gts2l){
  sum1 <- rowSums(gts2l[,1:2])
  sum2 <- rowSums(gts2l[,3:4])
  cMat <- matrix(
    c(
      sum(sum1==0 & sum2==0), sum(sum1==0 & sum2==1), sum(sum1==0 & sum2==2),
      sum(sum1==1 & sum2==0), sum(sum1==1 & sum2==1), sum(sum1==1 & sum2==2),
      sum(sum1==2 & sum2==0), sum(sum1==2 & sum2==1), sum(sum1==2 & sum2==2)
    ),
    3,3
  )
  if(0 %in% colSums(cMat)) cMat <- cMat[, colSums(cMat)!=0]
  if(0 %in% rowSums(cMat)) cMat <- cMat[, rowSums(cMat)!=0]
  cMat
}

makeContTab2l(gt2l01)


getChisqVal <- function(cTab){
  expect <- (rowSums(cTab)/sum(cTab)) %o% colSums(cTab)
  sum(
    (expect - cTab)^2 / expect
  )
}

getChisqValAndDf <- function(cTab){
  expect <- (rowSums(cTab)/sum(cTab)) %o% colSums(cTab)
  csq <- sum(
    (expect - cTab)^2 / expect
  )
  c(csq=csq,df=prod(dim(cTab)-1))
}

# takes as input the result of getChisqValAndDf
getChisqP <- function(x){
  y <- unname(x)
  c(p=1-pchisq(q=y[1], df=y[2]),
    chisq=y[1],
    df=y[2]
  )
}

getChisqVal(makeContTab2l(gt2l01))
getChisqValAndDf(makeContTab2l(simGts2loci()))
getChisqP(getChisqValAndDf(makeContTab2l(simGts2loci())))
vals <- t(sapply(1:10000, function(x) getChisqP(getChisqValAndDf(makeContTab2l(simGts2loci(nInd=100))))))
boxplot(p ~ df, data=vals)
summary(vals)
hist(vals[,1])
table(vals[,1] < 0.05)
plot(vals[,2:1], log="")
#cVals <- sapply(1:10000, function(x) getChisqVal(makeContTab2l(simGts2loci(nInd=100))))
cValsDfs <- t(sapply(1:10000, function(x) getChisqValAndDf(makeContTab2l(simGts2loci(nInd=100)))))
dfEmpirical <- mean(cValsDfs[,2])
table(cValsDfs[,2])
hist(cValsDfs[,1], freq=F, breaks="FD")
curve(dchisq(x, df=dfEmpirical), 0, 15, add=T)
curve((0.6331*dchisq(x, df=4) + 0.3669*dchisq(x, df=2)), 0, 15, add=T)

abline(v=qchisq(0.95, df=1), lty=2)
table(cVals > qchisq(0.95, df=1))



curve((dchisq(x, df=2) + dchisq(x, df=3))/2, 0, 15)
curve(dchisq(x, df=2.5), 0, 15, add=T, col=2, lty=2)

curve((dchisq(x, df=0) + dchisq(x, df=1))/2, 0, 15)
curve(dchisq(x, df=0.5), 0, 15, add=T, col=2, lty=2)
abline(v=1)
abline(v=1/2)
abline(v=1/3)
abline(v=1/4)
