library(lme4)
library(Matrix)
library(MASS)
library(ggplot2)


coffebra = read.table("simBRAx1.txt", header = T)
coffebra.cont <- coffebra
# new: nn & grpunt
nn=dim(coffebra.cont)[1]
grpunt=12*5


righe.cont <- 1:2
f1 <- 2.5
f2 <- 2

coffebra.cont$logV[righe.cont]=f1*coffebra$logV[righe.cont]
coffebra.cont$logQ[righe.cont]=f2*coffebra$logQ[righe.cont]

plot(coffebra.cont$logQ, coffebra.cont$logV, cex = 1.5, pch = 20, xlab = "W", ylab = "V")
abline(lm(logV ~ logQ, data=coffebra.cont))

fit.coffebra.cont = lmer(logV ~ logQ + t +m1 + m2+m3+m4+m5+ m6+m7+m8+m9+m10+m11+(1 | country), data = coffebra.cont, REML=F)


y = coffebra.cont$logV
ngroups = max(fit.coffebra.cont@Gp)
X = fit.coffebra.cont@pp$X

Z = t(fit.coffebra.cont@pp$Zt)
beta = fit.coffebra.cont@beta
ttest = beta/diag(summary(fit.coffebra.cont)$vcov)^0.5
sigma2U = summary(fit.coffebra.cont)$varcor[[1]][1]
sigma2eps = summary(fit.coffebra.cont)$sigma^2
loglikhat = summary(fit.coffebra.cont)$logLik[1]
Vhat = sigma2U * Z %*% t(Z) + sigma2eps * diag(length(y))
Ghat = sigma2U * diag(max(fit.coffebra.cont@Gp))
### as.matrix(Vhat) to het the matrix

i1 = t(X)%*% ginv(as.matrix(Vhat)) %*% X
i1inv = ginv(i1)
i2 = t(X)%*% ginv(as.matrix(Vhat)) %*% y
### The product i1inv %*% i2 is identical to beta

# ranef : stima dei random effects
### identical to ranef(ranef(fit.coffebra.cont)$pigname[[1]])
uhat = Ghat%*%t(Z)%*%ginv(as.matrix(Vhat)) %*% (y - X%*%beta)
yhat = X%*%beta + Z%*%uhat
res = y - yhat

### pointer = matrix(rep(1:7, ngroups), nrow = 7, ncol = ngroups)
### Generate a matrix of starting time-index
positions = matrix(1:length(y), nrow = max(coffebra.cont$t), ncol = ngroups)

#temp = matrix(NA, nrow = 10, ncol = length(y)/(ncol(X)+1))
temp = NULL
### for(l in 1:nrow(temp)) {
for(l in 1:50) {
    subset=matrix(NA, nrow=ncol(X), ncol = ngroups)
    set.seed(l)
    sss = sample(1:(max(coffebra.cont$t) - ncol(X)+1), ngroups, replace = T)
    for (ii in 1:ncol(positions)) {
        subset[,ii] = positions[sss[ii]:(sss[ii]+ncol(X)-1), ii]
    }
    temp = rbind(as.vector(subset), temp)
}

# temp = matrix(NA, nrow = length(y)/3, ncol = 7)
# subset = 1:3
# for(l in 1:7){
	# if(l == 1) subset = subset
	# else subset = subset+1
# for(i in 1:ngroups) {
	# if (i == 1) out = subset
	# else out = c(out, subset+(i-1)*9)
# }
# temp[, l] = out
# }


resorder =matrix(NA, nrow = length(y), ncol = nrow(temp))
for(j in 1:nrow(temp)) {
    cat(j, "\n")
    mysub = temp[j, ]
    fit.coffebra.cont.sub = lmer(logV ~ logQ + t +m1 + m2+m3+m4+m5+ m6+m7+m8+m9+m10+m11+(1 | country), data = coffebra.cont, REML=F, subset = mysub)
    ngroups = max(fit.coffebra.cont.sub@Gp)
    y.sub = fit.coffebra.cont.sub@resp$y
    X.sub = fit.coffebra.cont.sub@pp$X
    Z.sub = t(fit.coffebra.cont.sub@pp$Zt)
    beta.sub = fit.coffebra.cont.sub@beta
    sigma2U.sub = summary(fit.coffebra.cont.sub)$varcor[[1]][1]
    sigma2eps.sub = summary(fit.coffebra.cont.sub)$sigma^2
    loglikhat.sub = summary(fit.coffebra.cont.sub)$logLik[1]
    Ghat.sub = sigma2U.sub * diag(ngroups)

### This is to check if beta is obtained
    Vhat.sub = sigma2U.sub * Z.sub %*% t(Z.sub) + sigma2eps.sub * diag(length(y.sub))
### as.matrix(Vhat) to het the matrix
    i1 = t(X.sub)%*% ginv(as.matrix(Vhat.sub)) %*% X.sub
    i1inv = ginv(i1)
    i2 = t(X.sub)%*% ginv(as.matrix(Vhat.sub)) %*% y.sub
### The product i1inv %*% i2 is identical to beta

### Ste equal to zero rows unused
    not.mysub = setdiff(1:length(y), mysub)
    Z.reduced = Z
    Z.reduced[not.mysub, ] = 0

    Vhat.sub1 = sigma2U.sub * Z.reduced %*% t(Z.reduced) + sigma2eps.sub * diag(length(y))
    Ghat.sub1 = sigma2U.sub * diag(ngroups)


### identical to ranef(ranef(fit.coffebra.cont.sub)$pigname[[1]])
    uhat.sub = Ghat.sub1%*%t(Z.reduced)%*%ginv(as.matrix(Vhat.sub1)) %*% (y - X%*%beta.sub)
###    yhat.sub = X%*%beta.sub + Z.reduced%*%uhat.sub
    yhat.sub = X%*%beta.sub + Z%*%uhat.sub
    res.sub = y - yhat.sub

    resorder[,j] = sort(res.sub^2)
}

bsb = temp[which.min(apply(resorder, 2, median)), ]

### Storage of quantities
res.fwd=matrix(NA, nrow = length(y), ncol = length(y)-length(bsb) + 1)
beta.fwd = matrix(NA, nrow = ncol(X), ncol = length(y)-length(bsb) + 1)
ttest.fwd = beta.fwd
sigma2U.fwd = matrix(NA, nrow = 1, ncol = length(y)-length(bsb)+ 1)
sigma2eps.fwd = matrix(NA, nrow = 1, ncol = length(y)-length(bsb)+ 1)
row.inside.fwd = matrix(NA, nrow = length(y),ncol = length(y)-length(bsb)+ 1)
uhat.fwd = matrix(NA, nrow=ngroups, ncol = length(y)-length(bsb)+ 1)


############################  START OF THE FS

# new: seq code and for code
fscycle=seq(70,nn,5)



for(m in fscycle) {
    cat("Subset size m of the FWD = ", m, "\n")
    if (m == length(bsb)) { ### first step of the forward search
	step = 1
	ttt = bsb
        init.fit = lmer(logV ~ logQ + t +m1 + m2+m3+m4+m5+ m6+m7+m8+m9+m10+m11+(1 | country), data = coffebra.cont, REML=F, subset = ttt)
###        init.fit = lmer(w ~ t +(1 | pigname), data = pig, REML=F, subset = ttt)
        ngroups = max(init.fit@Gp)
        y.sub = init.fit@resp$y
        X.sub = init.fit@pp$X
        Z.sub = t(init.fit@pp$Zt)
        beta.sub = init.fit@beta
        ttest.sub = beta.sub/diag(summary(init.fit)$vcov)^0.5
        sigma2U.sub = summary(init.fit)$varcor[[1]][1]
        sigma2eps.sub = summary(init.fit)$sigma^2
        loglikhat.sub = summary(init.fit)$logLik[1]
        Ghat.sub = sigma2U.sub * diag(ngroups)

### Ste equal to zero rows unused
	not.mysub = setdiff(1:length(y), ttt)
        Z.reduced = Z
        Z.reduced[not.mysub, ] = 0

        Vhat.sub1 = sigma2U.sub * Z.reduced %*% t(Z.reduced) + sigma2eps.sub * diag(length(y))
	Ghat.sub1 = sigma2U.sub * diag(ngroups)


### identical to ranef(ranef(fit.coffebra.cont.sub)$pigname[[1]])
        uhat.sub = Ghat.sub1%*%t(Z.reduced)%*%ginv(as.matrix(Vhat.sub1)) %*% (y - X%*%beta.sub)
###        yhat.sub = X%*%beta.sub + Z.reduced%*%uhat.sub
        yhat.sub = X%*%beta.sub + Z%*%uhat.sub
        res.sub = y - yhat.sub
        posm = which(order(res.sub^2)<= m)
        res.fwd[, step] = as.vector(res.sub)
        beta.fwd[, step] = as.vector(beta.sub)
        ttest.fwd[, step] = as.vector(ttest.sub)
        sigma2U.fwd[1, step] = as.vector(sigma2U.sub)
        sigma2eps.fwd[1, step] = as.vector(sigma2eps.sub)
        row.inside.fwd[1:length(bsb), step] = bsb
        uhat.fwd[, step] = as.vector(uhat.sub)
    }
    else{
        step = step+1
        
        # new: ordering a data matrix with countries
        tmpres=as.data.frame(cbind(1:nn,as.vector(coffebra.cont$country),as.vector(res.sub^2)))
        names(tmpres)=c("unit","country","res2")
        
        # new: ordering two levels
        otmpres=tmpres[order(tmpres$country, tmpres$res2),]
        # we take the first m/5 units for each country
        ordbsb=as.vector(otmpres$unit)
        ordind=c((1:(m/5))+0*grpunt, (1:(m/5))+1*grpunt, (1:(m/5))+2*grpunt, (1:(m/5))+3*grpunt, (1:(m/5))+4*grpunt)
        
        nbsb=ordbsb[ordind]
        
        # original
        # nbsb= which(order(res.sub^2)<= m+1)
        
        
        ttt = nbsb
        init.fit = lmer(logV ~ logQ + t +m1 + m2+m3+m4+m5+ m6+m7+m8+m9+m10+m11+(1 | country), data = coffebra.cont, REML=F, subset = ttt)
###        init.fit = lmer(w ~ t +(1 | pigname), data = pig, REML=F, subset = ttt)
        check.ngroups = max(init.fit@Gp)
        if(check.ngroups < ngroups) {
        	#m = m+1
        	}
        else {
        cat("\n", "Progression = ", step, "\n")

        y.sub = init.fit@resp$y
        X.sub = init.fit@pp$X
        Z.sub = t(init.fit@pp$Zt)
        beta.sub = init.fit@beta
        ttest.sub = beta.sub/diag(summary(init.fit)$vcov)^0.5
        sigma2U.sub = summary(init.fit)$varcor[[1]][1]
        sigma2eps.sub = summary(init.fit)$sigma^2
        loglikhat.sub = summary(init.fit)$logLik[1]
        Ghat.sub = sigma2U.sub * diag(ngroups)

### Ste equal to zero rows unused

	    not.mysub = setdiff(1:length(y), ttt)
	    Z.reduced = Z

        if(length(not.mysub)>0)
            Z.reduced[not.mysub, ] = 0

        Vhat.sub1 = sigma2U.sub * Z.reduced %*% t(Z.reduced) + sigma2eps.sub * diag(length(y))
	    Ghat.sub1 = sigma2U.sub * diag(ngroups)


### identical to ranef(ranef(fit.coffebra.cont.sub)$pigname[[1]])
        uhat.sub = Ghat.sub1%*%t(Z.reduced)%*%ginv(as.matrix(Vhat.sub1)) %*% (y - X%*%beta.sub)
###        yhat.sub = X%*%beta.sub + Z.reduced%*%uhat.sub
        yhat.sub = X%*%beta.sub + Z%*%uhat.sub
        res.sub = y - yhat.sub
        res.fwd[, step] = as.vector(res.sub)
        beta.fwd[, step] = as.vector(beta.sub)
        ttest.fwd[, step] = as.vector(ttest.sub)
        sigma2U.fwd[1, step] = as.vector(sigma2U.sub)
        sigma2eps.fwd[1, step] = as.vector(sigma2eps.sub)
        row.inside.fwd[1:length(nbsb), step] = nbsb
        uhat.fwd[, step] = as.vector(uhat.sub)

        }
    }
}


steps=1:step

plot(length(bsb):length(y), (res.fwd[1,])^2, type = "n", ylim = range(res.fwd^2, na.rm = T), xlab = "Subset size m", ylab = "Squared Residuals")
for (i in 1:nrow(res.fwd)) 
  lines(fscycle, (res.fwd[i, steps])^2, lty = i, col = i)

dev.print(device = postscript, file = "C:/ASIM/multilevelFS/figures/squared.resid.sim.eps", height = 18, width = 30, onefile = T)


plot(length(bsb):length(y), beta.fwd[1,], type = "n", ylim = range(beta.fwd, na.rm = T), xlab = "Subset size m", ylab = "Estimates of fixed effects")
for (i in 1:nrow(beta.fwd)) 
  lines(fscycle, beta.fwd[i, steps], lty = i, col = i)

dev.print(device = postscript, file = "C:/ASIM/multilevelFS/figures/bhat.sim.eps", height = 18, width = 30, onefile = T)



plot(length(bsb):length(y), ttest.fwd[1,], type = "n", ylim = range(ttest.fwd, na.rm = T), xlab = "Subset size m", ylab = "T Test")
for (i in 1:nrow(beta.fwd)) 
  lines(fscycle, ttest.fwd[i, steps], lty = i, col = i)
abline(h = c(-2,2), lwd = 3, lty = 2, col = "blue")

dev.print(device = postscript, file = "C:/ASIM/multilevelFS/figures/ttest.bhat.sim.eps", height = 18, width = 30, onefile = T)


plot(length(bsb):length(y), uhat.fwd[1,], type = "n", ylim = range(uhat.fwd, na.rm = T), xlab = "Subset size m", ylab = "Predicted random effects")
for (i in 1:nrow(uhat.fwd)) 
  lines(fscycle, uhat.fwd[i, steps], lty = i, col = i)
abline(h = 0, lwd = 5, lty = 2, col = "blue")

dev.print(device = postscript, file = "C:/ASIM/multilevelFS/figures/predicted.reff.sim.eps", height = 18, width = 30, onefile = T)


########################## 

plot(length(bsb):length(y), sigma2eps.fwd[1,], type = "n", ylim = range(sigma2eps.fwd, na.rm = T), xlab = "Subset size m", ylab = "Estimated variances")
for (i in 1:nrow(sigma2eps.fwd)) lines(length(bsb):length(y), sigma2eps.fwd[i, ], lty = i, col = i)

plot(length(bsb):length(y), sigma2U.fwd[1,], type = "n", ylim = range(sigma2U.fwd, na.rm = T), xlab = "Subset size m", ylab = "Estimated variances")
for (i in 1:nrow(sigma2U.fwd)) lines(length(bsb):length(y), sigma2U.fwd[i, ], lty = i, col = i)

