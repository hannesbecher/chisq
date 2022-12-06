

# generate data with LD ---------------------------------------------------
library(genetics)



# D = p(AB) - p(A)*p(B)


# AB   Ab   p(A)
# aB   ab   p(a)
# p(B) p(b)

# p(ab) = 1 - pAB - (pA-pAB) - (pB-pAB)
# p(ab) = 1  -pAB - pA+pAB -pB+pAB
# p(ab) = 1 + pAB- pA -pB

# if
# 0 > 1 + pAB- pA -pB
# then
# pA+pB-1 > pAB


# generate a 2x2 table of gemete probabilities given pAB
makeContTabABAB <- function(pA, pB, pAB){
  if(pAB > pA) stop("pAB is greater than pA")
  if(pAB > pB) stop("pAB is greater than pB")
  if(pA+pB - 1 > pAB) stop("pA+pB - 1 is greater then pAB")
  matrix(c(pAB, pB-pAB,  pA-pAB, (1-pA)-(pB-pAB)),
         2,2)
}

# pAB = pA*pB + D
# pAb = pA- pAB
# if
# pAb < 0
# then 
# pA-pAB < 0
# pA-(pA*pB + D) < 0
# pA- pA*pB - D < 0
# D > pA-pA*pB

# generate a 2x2 table of gemete probabilities given D
makeContTabABD <- function(pA, pB, D){
  if(D < -pA*pB) stop("D is less than - pA * pB")
  #if(D > 1 - pA*pB) stop("D is greater  than 1 - pA * pB")
  if(D > pA - pA*pB) stop("D is greater  than pA - pA * pB")
  if(D > pB - pA*pB) stop("D is greater  than pB - pA * pB")
  pAB <- pA*pB +D
  matrix(c(pAB, pB-pAB,  pA-pAB, (1-pA)-(pB-pAB)),
         2,2)
}
makeContTabABD(0.1, .2, .07)



makeContTabABAB(.1, .2, .05)
as.vector(makeContTabABAB(.1, .2, .05)) # does it by column
makeContTabABAB(.1, .2, .11)
makeContTabABAB(.5, .5, .01)
makeContTabABAB(.8, .9, .71)
makeContTabABAB(.8, .9, .70)
makeContTabABAB(.8, .9, .69)
makeContTabABAB(.2, .1, .29)


# converts a 2x2 probability matrix to a vector naming the items
gameteProbsFromMat <- function(contTab){
  if(!all(dim(contTab) == c(2,2))) stop("`contTab` must be a 2x2 matrix.")
  prbs <- as.vector(contTab)
  names(prbs) <- c("p(AB)", "p(aB)", "p(Ab)", "p(ab)")
  if(sum(prbs) != 1) stop("Gamete probs do not sum to 1!")
  prbs
}


gameteProbsFromMat(makeContTabABAB(.1, .2, .09))

# generate diploid genotypes with LD
# returns 0/1/2 GTs for 2 loci
# either D of pAB must be specified
makeLdData <- function(N=100, pA=0.1, pB=0.2, pAB, D){
  pABexp <- pA*pB
  if(!missing(pAB) && missing(D)){
    prbs <- gameteProbsFromMat(makeContTabABAB(pA=pA, pB=pB,pAB=pAB))
    return(t(matrix(c(1, 1, 0, 1,  1, 0, 0, 0), 2, 4) %*% rmultinom(N, 2, prbs)))
    
  } else if(missing(pAB) && !missing(D)) {
    prbs <- gameteProbsFromMat(makeContTabABD(pA=pA, pB=pB,D=D))
    return(t(matrix(c(1, 1, 0, 1,  1, 0, 0, 0), 2, 4) %*% rmultinom(N, 2, prbs)))
  } else stop("Either D or pAB must be specified")
}




makeLdData(N=10, pA=0.7, pB=.8, pAB=.6)
makeLdData(N=10, pA=0.7, pB=.8, D=.11)
makeLdData(N=10, pA=0.7, pB=.8, D=0)


gts <- makeLdData(N=100, pA=0.7, pB=.8, pAB=.6)
gt1 <- as.genotype.allele.count(gts[,1])
gt2 <- as.genotype.allele.count(gts[,2])
LD(gt1, gt2)

gts02 <- makeLdData(N=100, pA=0.7, pB=.8, D=.14)
gt21 <- as.genotype.allele.count(gts02[,1])
gt22 <- as.genotype.allele.count(gts02[,2])
LD(gt21, gt22)


# No LD but deviation from HWE -------------------------------------------

# Haldane 1954, p 633, dropping n to get probabilities:
# pAA = p(p+fq) = p(p+f(1-p)) = p(p + f-fp) = p^2+pf-fp^2
# pAa = 2(1-f)pq
# paa = q(q+fp)


makeGtsDevHwe <- function(N, p, f){
  q <- 1-p
  pAA <- p*(p+f*q)
  pAa <- 2*(1-f)*p*q
  paa <- q*(q+f*p)
  as.vector(c(2,1,0) %*% rmultinom(n=N, size=1, prob = c(pAA, pAa, paa)) )
}


makeGtsDevHweTwoCol <- function(N, p, f){
  q <- 1-p
  pAA <- p*(p+f*q)
  pAa <- 2*(1-f)*p*q
  paa <- q*(q+f*p)
  gts <- as.vector(c(2,1,0) %*% rmultinom(n=N, size=1, prob = c(pAA, pAa, paa)) )
  matrix(c(0,0,1,0,1,1),ncol=2)[gts+1,]
}
makeGtsDevHweTwoCol(10, .1, .1)
makeGtsDevHwe(100, 0.3, .9)
makeGtsDevHwe(100, 0.3, 0)
gt00 <- as.genotype.allele.count(makeGtsDevHwe(100, 0.3, -.3))
gt00 <- as.genotype.allele.count(makeGtsDevHwe(100, 0.3, 0.8))
HWE.chisq(gt00)
HWE.exact(gt00)

pvs <- sapply(1:1000, function(f){
gt001 <- as.genotype.allele.count(makeGtsDevHwe(100, 0.1, 0.8))
gt002 <- as.genotype.allele.count(makeGtsDevHwe(100, 0.1, 0.8))
LD(gt001, gt002)$`P-value`
}
)
hist(pvs)
quantile(pvs, 0.05)
table(pvs < 0.05)/1000
# with deviation from HWE in both loci, there is nflation of significant chi-sq p-vals




pvs <- sapply(1:1000, function(f){
  gt001 <- as.genotype.allele.count(makeGtsDevHwe(1000, 0.1, 0))
  gt002 <- as.genotype.allele.count(makeGtsDevHwe(1000, 0.2, 0))
  c(pLD=LD(gt001, gt002)$`P-value`,
    HWE.chisq(gt001)$p.value,
    HWE.chisq(gt002)$p.value
  )
}
)

hist(pvs[1,], breaks=50)
hist(pvs[2,])
hist(pvs[3,])
library(rgl)

plot3d(pvs[1,], pvs[2,], pvs[3,])
cor(pvs[1,], pvs[2,])
cor(pvs[1,], pvs[3,])
cor(pvs[2,], pvs[3,])
table(pvs < 0.05)/1000

library(pegas)
?LD2
?loci

pegas::alleles2loci(gt001)
?pegas::loci2alleles()
loci2alleles(matrix(c(1, 0, 1, 0, 0, 0, 1, 1), ncol=2))
x <- matrix(c("A", "A", "A", "a"), 2)
colnames(x) <- c("Loc1", NA)
y <- alleles2loci(x)
print(y, details = TRUE)
loci2alleles(y)
LD2(y, y)
?LD2
data(jaguar)
str(jaguar)
?jaguar

pvls <- sapply(1:1000, function(x){
gts4col <- cbind(makeGtsDevHweTwoCol(1000, .1, .2),
      makeGtsDevHweTwoCol(1000, .2, .2)
)
LD2(alleles2loci(gts4col))$T2[3]
}
)
hist(pvls, breaks=20)
LD2


#ldd <- makeLdData(N=1000, pA=0.1, pB=0.2, pAB=0.02)
ldd <- makeLdData(N=1000, pA=0.1, pB=0.2, D=0)
d4c <- cbind(
  matrix(c(0,0,1,0,1,1),ncol=2)[ldd[,1]+1,],
      matrix(c(0,0,1,0,1,1),ncol=2)[ldd[,2]+1,]
      )
pegas::LD2(alleles2loci(d4c))


gt1 <- as.genotype.allele.count(ldd[,1])
gt2 <- as.genotype.allele.count(ldd[,2])
genetics::LD(gt1, gt2)

# LD2 from pegas is based on
# Zaykin, D. V., Pudovkin, A. and Weir, B. S. (2008)
# Correlation-based inference for linkage disequilibrium with multiple alleles.
# Genetics, 180, 533–545.
