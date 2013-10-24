load.data <- function() {
  
  d <- read.table("data2.txt", header=T, sep="\t", stringsAsFactors=F)
  d$powt <- factor(d$powt)
  d$population <- factor(d$population)
  d$genotyp <- factor(d$genotyp)
  d$rwness <- factor(d$rwness)
  d$shape <- factor(d$shape)
  d <<- cbind(d, label=paste(d[,3],d[,18],d[,19], sep=":"))
}

generate.name.table <- function() {

  u <- unique(d[,c(3,18,19)])
  u <- u[order(u[,1], u[,2], u[,3]),]
  u <- cbind(u, label=paste(u[,1], u[,2], u[,3], sep=":"))
  u <- cbind(u, labels=paste(u[,1], u[,2], u[,3], sep=":"))
  levels(u$labels)[1] <- "MARESI"
  levels(u$labels)[10] <- "POMO"
  u <<- u
}


# Visualisation of data points
visualise.data <- function() {

  old.par <- par(mfrow=c(4, 4))
  for (i in 5:17) {
    plot(y=d[,i], x=d$label, main=names(d)[i], las=2, at=c(1, 3,4,5,6, 8,9,10,11, 13), col=c("indianred","olivedrab")[u$shape], names=u$labels)
    #mtext("Population", side = 1, line = 5)
  }
  par(old.par)
  
}

load.data()
generate.name.table()
visualise.data()

stderr <- function(x) sqrt(var(x)/length(x))

means.lines <- aggregate(d[5:17], list(Line=d$genotyp), mean)
ses.lines <- aggregate(d[5:17], list(Line=d$genotyp), stderr)

ms.lines <- data.frame(line=means.lines[,1])
for (i in 2:14) {
  ms.lines <- cbind(ms.lines, means.lines[,i], ses.lines[,i])
  colnames(ms.lines)[2*i-2] <- names(means.lines)[i]
  colnames(ms.lines)[2*i-1] <- paste(names(means.lines)[i], "(se)", sep="")
}
write.table(format(ms.lines, digits=3, nsmall=2), file="out/means.lines.txt", sep="\t", quote=F)

means <- aggregate(d[5:17], list(Shape=d$shape, Rwness=d$rwness, Pop=d$population), mean)
ses <- aggregate(d[5:17], list(Shape=d$shape, Rwness=d$rwness, Pop=d$population), stderr)

ms <- means[,1:3]
for (i in 4:16) {
  ms <- cbind(ms, means[,i], ses[,i])
  colnames(ms)[2*i-4] <- names(means)[i]
  colnames(ms)[2*i-3] <- paste(names(means)[i], "(se)", sep="")
}
write.table(format(ms, digits=3, nsmall=2), file="out/means.factors.txt", sep="\t", quote=F)



means.p <- means[c(1,10),]
means.p$Pop <- factor(means.p$Pop)
means.c <- means[-c(1,10),]
means.c$Pop <- factor(means.c$Pop)

ses.c <- ses[-c(1,10),]
ses.c$Pop <- factor(ses.c$Pop)
ses.p <- ses[c(1,10),]
ses.p$Pop <- factor(ses.p$Pop)

error.bar <- function(x, y, upper, lower=upper, length=0.1,...){
  if(length(x) != length(y) | length(y) !=length(lower) | length(lower) != length(upper))
    stop("vectors must be same length")
  arrows(x,y+upper, x, y-lower, angle=90, code=3, length=length, ...)
}
library(gplots)


# Print means for parents:
t <- means.p
te <- ses.p
old.par <- par(mfrow=c(4, 4))
for (i in 4:16) {
  max <- max(t[,i]) + 2*max(te[,i])
  b <- barplot2(t[,i], main=names(t)[i], xlab="Population", names=c("MARESI", "POMO"), ylim=c(0, max), 
                col=c("indianred","olivedrab"), plot.grid=T)
  error.bar(b, t[,i], 2*te[,i])
}
par(old.par)


# Print means for children:
t <- means.c
te <- ses.c
old.par <- par(mfrow=c(4, 4))
for (i in 4:16) {
  max <- max(t[,i]) + 2*max(te[,i])
  fxt <- ftable(xtabs(t[,i]~t$Shape+t$Rwness+t$Pop))
  fxte <- ftable(xtabs(te[,i]~t$Shape+t$Rwness+t$Pop))
  lo <- fxt-fxte
  up <- fxt+fxte
  barplot2(fxt, beside=T, space=c(0.2, 0.8), names=paste(t$Pop, t$Rwness, t$Shape, sep=":"), las=2, font.lab=2,
           plot.grid=T, ylim=c(0,max(up)), col=c("indianred","olivedrab"), border="black", main=names(t)[i],
           plot.ci=TRUE, ci.l=lo, ci.u=up)
  #mtext("Population", side = 1, line = 5)
}
par(old.par)


library(lme4)

for (j in c("MP", "MPS")) {
  
  dd <- d[d$population==j,]
  
  loopname <- paste("\n\n ----- ", j, " ----- \n", sep="")
  write(loopname, file="out/means.txt", append=T)
  write(loopname, file="out/means_random.txt", append=T)
  
  for (i in 5:17) {
    
    trait <- names(dd)[i]
    loopname <- paste("\n\n", j, " : ", trait, "\n", sep="")
    fname <- paste("out/", j, "_", trait, ".txt", sep="")
    print(loopname)
    write(loopname, file="out/means.txt", append=T)
    write(loopname, file="out/means_random.txt", append=T)
    
    formula <- paste(trait, "~rwness*shape+(1|powt)", sep="")
    m <- lmer(formula, dd)
    c <- coef(summary(m))
    a <- anova(m)
    a$Fpr <- 1 - pf(a[,4], 1, dim(m@frame)[1])
    
    
    errorVar <- attr(VarCorr(m), "sc")^2
    randomVar <- VarCorr(m)[["powt"]][1]
    
    dat <- with(dd, expand.grid(shape=unique(shape), rwness=unique(rwness)))
    dat[trait] <- predict(m, dat, REform=NA)
    
    ef <- model.matrix(~shape, dat)
    ef[,1] <- ef[,1] - ef[,2]
    means <- solve(t(ef) %*% ef) %*% t(ef) %*% dat[,trait]
    mean.shape <- cbind(shape=levels(dat$shape), means)
    colnames(mean.shape)[2] <- trait
  
    ef <- model.matrix(~rwness, dat)
    ef[,1] <- ef[,1] - ef[,2]
    means <- solve(t(ef) %*% ef) %*% t(ef) %*% dat[,trait]
    mean.rwness <- cbind(rwness=levels(dat$rwness), means)
    colnames(mean.rwness)[2] <- trait
    
    datR <- with(dd, expand.grid(shape=unique(shape), rwness=unique(rwness), powt=unique(powt)))
    datR[trait] <- predict(m, datR, REform=NA) #same values for diff powt
    datR$plusRandom <- predict(m, datR, REform=~(1|powt)) #diff values for diff powt
    #datR$preNULL <- predict(m, datR, REform=NULL) #diff values for diff powt
    #print(datR)
    
    write.table( format(dat, digits=4, nsmall=2), file="out/means.txt", append=T, quote=F, row.names=F, sep="\t")
    write.table(format(datR, digits=4, nsmall=2), file="out/means_random.txt", append=T, quote=F, row.names=F, sep="\t")
    
    write.table( round(a, digits=3), file=fname, quote=F, sep="\t")
    write("\n", file=fname, append=T)
    write.table( round(c, digits=3), file=fname, append=T, quote=F, sep="\t")
    write("\n", file=fname, append=T)
    write(paste("Error variance: ", errorVar, "\n"), file=fname, append=T)
    write(paste("Random effect (powt) variance: ", randomVar, "\n"), file=fname, append=T)
    write.table( format(mean.shape, digits=4, nsmall=2), file=fname, append=T, quote=F, row.names=F, sep="\t")
    write("\n", file=fname, append=T)
    write.table( format(mean.rwness, digits=4, nsmall=2), file=fname, append=T, quote=F, row.names=F, sep="\t")
    write("\n", file=fname, append=T)
    write.table( format(dat, digits=4, nsmall=2), file=fname, append=T, quote=F, row.names=F, sep="\t")
    write("\n", file=fname, append=T)
    write.table(format(datR, digits=4, nsmall=2), file=fname, append=T, quote=F, row.names=F, sep="\t")    
    write("\n", file=fname, append=T)
  }
  
}

analysis.for.populations <- function() {
  
  dd <- d[grep(d$population, pattern="MP"),]
  for (i in 5:17) {
    
    trait <- names(dd)[i]
    formula <- paste(trait, "~population*rwness*shape+(1|powt)", sep="")
    m <- lmer(formula, dd)
    c <- coef(summary(m))
    a <- anova(m)
    a$Fpr <- 1 - pf(a[,4], 1, dim(m@frame)[1])
   
    errorVar <- attr(VarCorr(m), "sc")^2
    randomVar <- VarCorr(m)[["powt"]][1]
    
    fname <- paste("out/", "all_", trait, ".txt", sep="")
    write.table( round(a, digits=3), file=fname, quote=F, sep="\t")
    write("\n", file=fname, append=T)
    write.table( round(c, digits=3), file=fname, append=T, quote=F, sep="\t")
    write("\n", file=fname, append=T)
    write(paste("Error variance: ", errorVar, "\n"), file=fname, append=T)
    write(paste("Random effect (powt) variance: ", randomVar, "\n"), file=fname, append=T)
    
    write(trait, file="out/all.txt", append=T)
    write.table( format(a, digits=3, nsmall=2), file="out/all.txt", append=T, quote=F, sep="\t")
    write("\n\n", file="out/all.txt", append=T)
    
    write(trait, file="out/all_c.txt", append=T)
    write.table( format(c, digits=3, nsmall=2), file="out/all_c.txt", append=T, quote=F, sep="\t")
    write("\n", file="out/all_c.txt", append=T)
    
    
    
    dat <- with(dd, expand.grid(population=unique(population), shape=unique(shape), rwness=unique(rwness)), powt=unique(powt))
    dat[trait] <- predict(m, dat, REform=NA)
    
    ef <- model.matrix(~shape, dat)
    ef[,1] <- ef[,1] - ef[,2]
    means <- solve(t(ef) %*% ef) %*% t(ef) %*% dat[,trait]
    mean.shape <- cbind(shape=levels(dat$shape), means)
    colnames(mean.shape)[2] <- trait
    
    ef <- model.matrix(~rwness, dat)
    ef[,1] <- ef[,1] - ef[,2]
    means <- solve(t(ef) %*% ef) %*% t(ef) %*% dat[,trait]
    mean.rwness <- cbind(rwness=levels(dat$rwness), means)
    colnames(mean.rwness)[2] <- trait
    
    datR <- with(dd, expand.grid(shape=unique(shape), rwness=unique(rwness), powt=unique(powt)))
    datR[trait] <- predict(m, datR, REform=NA) #same values for diff powt
    datR$plusRandom <- predict(m, datR, REform=~(1|powt)) #diff values for diff powt
    #datR$preNULL <- predict(m, datR, REform=NULL) #diff values for diff powt
    #print(datR)
    
    write.table( format(dat, digits=4, nsmall=2), file="out/means.txt", append=T, quote=F, row.names=F, sep="\t")
    write.table(format(datR, digits=4, nsmall=2), file="out/means_random.txt", append=T, quote=F, row.names=F, sep="\t")
    
  }
  
}

  

calc.means.manually <- function() {
  
  nd <- with(dd, expand.grid(shape=unique(shape), rwness=unique(rwness), l.ziarn=0))
  mm <- model.matrix(terms(m),nd)
  nd$l.ziarn <- mm %*% fixef(m)
  
  pvar1 <- diag(mm %*% tcrossprod(vcov(m),mm))
  tvar1 <- pvar1+VarCorr(m)$powt[1]  ## must be adapted for more complex models
  nd <- data.frame(nd
    , plo = nd$l.ziarn-2*sqrt(pvar1), phi = nd$l.ziarn+2*sqrt(pvar1)
    , tlo = nd$l.ziarn-2*sqrt(tvar1), thi = nd$l.ziarn+2*sqrt(tvar1)
  )
  
  library(ggplot2) # Plotting
  #plot confidence
  g0 <- ggplot(nd, aes(x=shape, y=l.ziarn, colour=rwness))+geom_point()
  g0 + geom_errorbar(aes(ymin = plo, ymax = phi))+
    labs(title="CI based on fixed-effects uncertainty ONLY")
  #plot prediction
  g0 + geom_errorbar(aes(ymin = tlo, ymax = thi))+
    ggtitle("CI based on FE uncertainty + RE variance")
}

nd <- with(dd, expand.grid(shape=unique(shape), rwness=unique(rwness), l.ziarn=0))
mm <- model.matrix(terms(m),nd)
nd$l.ziarn <- mm %*% fixef(m)




data("Orthodont",package="MEMSS")
fm1 <- lmer(
  formula = distance ~ age*Sex + (age|Subject)
  , data = Orthodont
)
newdat <- expand.grid(
  age=c(8,10,12,14)
  , Sex=c("Male","Female")
  , distance = 0
)


dd <- data.frame(a = gl(3,4), b = gl(4,1,12)) # balanced 2-way
options("contrasts")
mm1 <- model.matrix(~ a + b, dd)
mm2 <- model.matrix(~ a + b, dd, contrasts = list(a = "contr.sum"))
mm3 <- model.matrix(~ a + b, dd, contrasts = list(a = "contr.sum", b = "contr.poly"))
m.orth <- model.matrix(~a+b, dd, contrasts = list(a = "contr.helmert"))
crossprod(m.orth) # m.orth is  ALMOST  orthogonal

#contrasts test
{
utils::example(factor)
fff <- ff[, drop = TRUE]  # reduce to 5 levels.
contrasts(fff) # treatment contrasts by default
contrasts(C(fff, sum))
contrasts(fff, contrasts = FALSE) # the 5x5 identity matrix

contrasts(fff) <- contr.sum(5); contrasts(fff)  # set sum contrasts
contrasts(fff, 2) <- contr.sum(5); contrasts(fff)  # set 2 contrasts
# supply 2 contrasts, compute 2 more to make full set of 4.
contrasts(fff) <- contr.sum(5)[, 1:2]; contrasts(fff)

## using sparse contrasts: % useful, once model.matrix() works with these :
ffs <- fff
contrasts(ffs) <- contr.sum(5, sparse = TRUE)[, 1:2]; contrasts(ffs)
stopifnot(all.equal(ffs, fff))
contrasts(ffs) <- contr.sum(5, sparse = TRUE); contrasts(ffs)
}

x <- unique(getME(m, "X"))
est <- x %*% fixef(m)
est.cov <- x %*% vcov(m) %*% t(x)

xf <- model.matrix(~shape, dd)
xf[,1] <- xf[,1] - xf[,2]
xf

ms <- solve(t(xf) %*% xf) %*% t(xf)
means <- ms %*% est





library(gplots)   # Load the gplots graphics library  