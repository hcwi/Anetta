d <- read.table("data2.txt", header=T, sep="\t", stringsAsFactors=F)
d$powt <- factor(d$powt)
d$population <- factor(d$population)
d$genotyp <- factor(d$genotyp)
d$rwness <- factor(d$rwness)
d$shape <- factor(d$shape)

d <- cbind(d, label=paste(d[,3],d[,18],d[,19], sep=":"))


u <- unique(d[,c(3,18,19)])
u <- u[order(u[,1], u[,2], u[,3]),]
u <- cbind(u, label=paste(u[,1], u[,2], u[,3], sep=":"))
u <- cbind(u, labels=paste(u[,1], u[,2], u[,3], sep=":"))
levels(u$labels)[1] <- "MARESI"
levels(u$labels)[10] <- "POMO"

old.par <- par(mfrow=c(4, 4))
for (i in 5:17) {
  plot(y=d[,i], x=d$label, main=names(d)[i], las=2, at=c(1, 3,4,5,6, 8,9,10,11, 13), col=c("indianred","olivedrab")[u$shape], names=u$labels)
  #mtext("Population", side = 1, line = 5)
}
par(old.par)






means <- aggregate(d[5:17], list(Pop=d$population, Rwness=d$rwness, Shape=d$shape), mean)
sds <- aggregate(d[5:17], list(Pop=d$population, Rwness=d$rwness, Shape=d$shape), sd)
means.p <- means[c(1,10),]
means.p$Pop <- factor(means.p$Pop)
means.c <- means[-c(1,10),]
means.c$Pop <- factor(means.c$Pop)

sds.c <- sds[-c(1,10),]
sds.c$Pop <- factor(sds.c$Pop)
sds.p <- sds[c(1,10),]
sds.p$Pop <- factor(sds.p$Pop)

for (i in 4:16) {
  plot(y=means.c[,i], x=means.c[,1]:means.c[,2]:means.c[,3], ylab=names(means.c)[i], )
}


             