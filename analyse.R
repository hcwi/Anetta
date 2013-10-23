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

means <- aggregate(d[5:17], list(Shape=d$shape, Rwness=d$rwness, Pop=d$population), mean)
ses <- aggregate(d[5:17], list(Shape=d$shape, Rwness=d$rwness, Pop=d$population), stderr)

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


d.mp <- d[d$population=="MP",]
d.mps <- d[d$population=="MPS",]

for (dd in list(d.mp, d.mps)) {
  i = 6
  library(lme4)
  formula <- paste(names(dd)[i], "~rwness*shape+(1|powt)", sep="")
  m <- lmer(formula, dd)
  
  dat <- with(dd, expand.grid(shape=unique(shape), rwness=unique(rwness)))
  dat$preNA2 <- predict(m, dat, REform=NA)
  print(dat)
  
  datR <- with(dd, expand.grid(shape=unique(shape), rwness=unique(rwness), powt=unique(powt)))
  datR$prePowt <- predict(m, datR, REform=~(1|powt))
  datR$preNA <- predict(m, datR, REform=NA)
  datR$preNULL <- predict(m, datR, REform=NULL)
  print(datR)
  
  dat2 <- with(dd, expand.grid(shape=unique(shape), rwness=unique(powt)))
  dat2 <- with(dd, expand.grid(shape=unique(shape), powt=unique(powt), rwness=unique(rwness)))
  dat2$preNA <- predict(m, dat2, REform=~(1|powt))
  print(dat2)
  
  print("-------")
}



dd <- d.mps
i = 6
library(lme4)
formula <- paste(names(dd)[i], "~rwness*shape+(1|powt)", sep="")
lmer(formula, dd)

coef(summary(m))



x <- unique(getME(model, "X"))
est <- x %*% fixef(model)
est.cov <- x %*% vcov(model) %*% t(x)


data(Orthodont, package="MEMSS")
newdat <- expand.grid(age=c(8,10,12,14), Sex=c("Male","Female"), distance = 0)
fm1 = lmer(formula = distance ~ age*Sex + (age|Subject), data = Orthodont)
mm = model.matrix(terms(fm1),newdat)
newdat$distance = mm %*% fixef(fm1)
pvar1 <- diag(mm %*% tcrossprod(vcov(fm1),mm))
tvar1 <- pvar1+VarCorr(fm1)$Subject[1]
newdat2 <- data.frame(newdat, 
                      plo = newdat$distance-2*sqrt(pvar1), 
                      phi = newdat$distance+2*sqrt(pvar1), 
                      tlo = newdat$distance-2*sqrt(tvar1),
                      thi = newdat$distance+2*sqrt(tvar1))

newdat2$pred <- predict(fm1, newdat2, distance=0)









Age = c(50, 55, 60, 65, 70)                                     # Age groups  
Male = c(15.4, 24.3, 37.0, 54.6, 71.1)                          # Death rates for males  
Female = c(8.4, 13.6, 19.3, 35.1, 50.0)                         # Death rates for females  
Deathrate = matrix(c(Male,Female), nrow=length(Age), ncol=2, dimnames=list(Age, c("Male","Female")))         
Deathrate2 = t(Deathrate)                                       # Transpose the Deathrate matrix  

barplot(Deathrate2,                             # Data (bar heights) to plot  
        beside=TRUE,                            # Plot the bars beside one another; default is to plot stacked bars  
        space=c(0.2,0.8),                       # Amount of space between i) bars within a group, ii) bars between groups  
        names.arg=c("65-69", "60-64", "55-59", "50-54", "70-74"),    #Names for the bars  
        col=c("blue", "red"),                   # Color of the bars  
        border="black",                         # Color of the bar borders  
        main=c("Death rates in Virginia"),      # Main title for the plot  
        xlab="Age group",                       # X-axis label  
        ylab="Death rate",                      # Y-axis label  
        font.lab=2)    

legend("topleft",                               # Add a legend to the plot  
       legend=c("Male", "Female"),             # Text for the legend  
       fill=c("blue", "red"))                  # Fill for boxes of the legend  


library(gplots)   # Load the gplots graphics library  

# Generate (fake) confidence intervals (CI should be derived from the underlying data)  
cil <- Deathrate2 * 0.85  
ciu <- Deathrate2 * 1.15  

barplot2(Deathrate2,                            # Data (bar heights) to plot  
          beside=TRUE,                            # Plot the bars beside one another; default is to plot stacked bars  
          space=c(0.2,0.8),                       # Amount of space between i) bars within a group, ii) bars between groups  
          names.arg=c("65-69", "60-64", "55-59", "50-54", "70-74"),    #Names for the bars  
          col=c("blue", "red"),                   # Color of the bars  
          border="black",                         # Color of the bar borders  
          main=c("Death rates in Virginia"),      # Main title for the plot  
          xlab="Age group",                       # X-axis label  
          ylab="Death rate",                      # Y-axis label  
          font.lab=2,                             # Font to use for the axis labels: 1=plain text, 2=bold, 3=italic, 4=bold italic  
          plot.ci=TRUE,                           # Plot confidence intervals  
          ci.l=cil,                               # Lower values for the confidence interval  
          ci.u=ciu,                               # Upper values for the confidence interval  
          plot.grid=TRUE) 
             