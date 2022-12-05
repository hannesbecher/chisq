
setwd("~/git_repos/chisq/")

# The chi-squared distribution arises from summing squared normals

# four sample from the std normal
a01 <- rnorm(100)
a02 <- rnorm(100)
a03 <- rnorm(100)
a04 <- rnorm(100)


# sums of their squares
s01 <- a01^2                         # chisq, 1 df
s02 <- a01^2 + a02^2                 # chisq, 2 df
s03 <- a01^2 + a02^2 + a03^2         # chisq, 3 df
s04 <- a01^2 + a02^2 + a03^2 + a04^2 # chisq, 4 df

# plot histograms, show relative freqs
png("plot.png", width=7, height=5, units = "in", res=150)
hist(s01, freq = F, col="#00000000", border = 1, xlab="histograms - samples, curves - theoretical density",
     main="Some samples following the chi-squared")
hist(s02, add=T, freq = F, col="#00000000", border = 2)
hist(s03, add=T, freq = F, col="#00000000", border = 3)
hist(s04, add=T, freq = F, col="#00000000", border = 4)
legend("topright",
       col=1:4,
       lty=1,
       legend=c("1", "2", "3", "4"),
       title = "df")
# add theoretical density curves
curve(dchisq(x,1),0,10, col=1, add=T, lty=2)
curve(dchisq(x,2),0,10, add=T, col=2, lty=2)
curve(dchisq(x,3),0,10, add=T, col=3, lty=2)
curve(dchisq(x,4),0,10, add=T, col=4, lty=2)
dev.off()
