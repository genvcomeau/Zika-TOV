install.packages("ggplot2")
install.packages("psych")
install.packages("lsmeans")
install.packages("prevalence")
install.packages("plyr")
install.packages("aod")
install.packages("binom")
install.packages("wordcloud")
install.packages("stringr")

library(ggplot2)
library(psych)
library(lsmeans)
library(prevalence)
library(plyr)
library(aod)
library(binom)
library("wordcloud")
library(stringr)

tov <- read.csv("C:/Users/Genevieve/Documents/R/TOVR1.csv")
tov.noNA <- read.csv("C:/Users/Genevieve/Documents/R/TOVnoNA.csv")

tov$pos <- ifelse(tov$Ct.Mean < 35.00, 1, 0)
tov$sampfull <- tov$Sample
tov$Sampfull <- paste(tov$sampfull, " ", tov$Stage)
tov.uniq <- tov[!duplicated(tov$Sampfull),]

tov.noNA$sampfull <- tov.noNA$Sample
tov.noNA$Sampfull <- paste(tov.noNA$sampfull, " ", tov.noNA$Stage)
tov.noNA.uniq <- tov.noNA[!duplicated(tov.noNA$Sampfull),]

tov.noNA.uniq.1 <- sapply(strsplit(as.character(tov.noNA.uniq$Sample), ""), tail, 3)
tov.noNA.uniq$SampType <- tov.noNA.uniq.1[1,]

#Sample distribution counts
xtabs(~ Cycle + SampType, tov.noNA.uniq)


mean(tov.uniq$Ct.Mean, na.rm = TRUE)
# 38.00604

mean(tov.uniq$Ct.Dev, na.rm = TRUE)
#0.8540

mean(tov.uniq$R2, na.rm = TRUE)
#0.9211

pos <- subset(tov, tov$Ct.Mean < 35.01)
tov.pos.uniq <- pos[!duplicated(pos$Sample),]

binom.confint()
#G2 im 33 pup probs artifact

#freq barplot of pos
barplot(table(tov.pos.uniq$Temp), main = "Positive Cases by Temperature", xlab = "Temperature")
barplot(table(tov.pos.uniq$Cycle), main = "Positive Cases by Gonotropic Cycle", xlab = "Cycle")
barplot(table(tov.pos.uniq$Day), main = "Positive Cases by Day", xlab = "Day")

#plot of cases by day
pos.day <- matrix(c(0,3,7,14,1,10,2,4), ncol = 2)
colnames(pos.day) <- c("Day", "Positive Cases")
plot(pos.day[,1], pos.day[,2])

#Subsetting, graphing positives
#,text(Ct.Mean, R2,labels = tov.G1.3.uniq$Sample, cex = 0.7, pos = 2))

tov.G1im <- subset(tov, tov$Stage == "G1 im")
tov.G1im.uniq <- tov.G1im[!duplicated(tov.G1im$Sample),]

plot(y = tov.G1im.uniq$R2, x =tov.G1im.uniq$Ct.Mean, pch = 1,cex=0.75, col ="grey", 
     xlim=c(20,45), ylim=c(0.7,1), ylab="R2 value", xlab="Ct Mean value", main="G1 immatures")
with(subset(tov.G1im.uniq, tov.G1im.uniq$Ct.Mean>=35.01), 
     points(y = R2, x = Ct.Mean,
            pch=20, cex = 0.75, col="blue"))
with(subset(tov.G1im.uniq, tov.G1im.uniq$Ct.Mean<35.01), 
     points(y = R2, x = Ct.Mean, 
            pch=20, cex = 0.75, col="red"),
     textxy(Ct.Mean, R2,labels = Sample, cex = 0.7, pos = 2))
abline(v=35.0, col="black", lty=2, lwd=2)


tov.G1.3 <- subset(tov, tov$Stage == "G1 3")
tov.G1.3.uniq <- tov.G1.3[!duplicated(tov.G1.3$Sample),]
tov.G1.3.uniq.pos <- subset(tov.G1.3.uniq, tov.G1.3.uniq$Ct.Mean < 35.01)

with(tov.G1.3.uniq, plot(y = tov.G1.3.uniq$R2, x =tov.G1.3.uniq$Ct.Mean, pch = 1,cex=0.75, col ="grey", 
     xlim=c(20,45), ylim=c(0.7,1), ylab="R2 value", xlab="Ct Mean value", main="G1 3 day"))
text(tov.G1.3.uniq.pos[,c(4,9)], labels = tov.G1.3.uniq.pos[,2], pos = 2)
with(subset(tov.G1.3.uniq, tov.G1.3.uniq$Ct.Mean>=35.01), 
     points(y = R2, x = Ct.Mean,pch=20, cex = 0.75, col="blue"))
with(subset(tov.G1.3.uniq, tov.G1.3.uniq$Ct.Mean<35.01), 
     points(y = R2, x = Ct.Mean,pch=20, cex = 0.75, col="red"))
abline(v=35.0, col="black", lty=2, lwd=2)


textplot(tov.G1.3.uniq$Ct.Mean, tov.G1.3.uniq$R2, tov.G1.3.uniq$Sample)
with(subset(tov.G1.3.uniq, tov.G1.3.uniq$Ct.Mean>=35.01), 
     points(y = R2, x = Ct.Mean,pch=20, cex = 0.75, col="blue"))
with(subset(tov.G1.3.uniq, tov.G1.3.uniq$Ct.Mean<35.01), 
     points(y = R2, x = Ct.Mean,pch=20, cex = 0.75, col="red"))
abline(v=35.0, col="black", lty=2, lwd=2)



tov.G1.7 <- subset(tov, tov$Stage == "G1 7")
tov.G1.7.uniq <- tov.G1.7[!duplicated(tov.G1.7$Sample),]

plot(y = tov.G1.7.uniq$R2, x =tov.G1.7.uniq$Ct.Mean, pch = 1,cex=0.75, col ="grey", 
     xlim=c(20,45), ylim=c(0.7,1), ylab="R2 value", xlab="Ct Mean value", main="G1 7 day")
with(subset(tov.G1.7.uniq, tov.G1.7.uniq$Ct.Mean>=35.01), 
     points(y = R2, x = Ct.Mean,pch=20, cex = 0.75, col="blue"))
with(subset(tov.G1.7.uniq, tov.G1.7.uniq$Ct.Mean<35.01), 
     points(y = R2, x = Ct.Mean,pch=20, cex = 0.75, col="red"))
abline(v=35.0, col="black", lty=2, lwd=2)


tov.G1.14 <- subset(tov, tov$Stage == "G1 14")
tov.G1.14.uniq <- tov.G1.14[!duplicated(tov.G1.14$Sample),]

plot(y = tov.G1.14.uniq$R2, x =tov.G1.14.uniq$Ct.Mean, pch = 1,cex=0.75, col ="grey", 
     xlim=c(20,45), ylim=c(0.7,1), ylab="R2 value", xlab="Ct Mean value", main="G1 14 day")
with(subset(tov.G1.14.uniq, tov.G1.14.uniq$Ct.Mean>=35.01), 
     points(y = R2, x = Ct.Mean,pch=20, cex = 0.75, col="blue"))
with(subset(tov.G1.14.uniq, tov.G1.14.uniq$Ct.Mean<35.01), 
     points(y = R2, x = Ct.Mean,pch=20, cex = 0.75, col="red"))
abline(v=35.0, col="black", lty=2, lwd=2)

tov.G2im <- subset(tov, tov$Stage == "G2 im")
tov.G2im.uniq <- tov.G2im[!duplicated(tov.G2im$Sample),]

plot(y = tov.G2im.uniq$R2, x =tov.G2im.uniq$Ct.Mean, pch = 1,cex=0.75, col ="grey", 
     xlim=c(20,45), ylim=c(0.7,1), ylab="R2 value", xlab="Ct Mean value", main="G2 immatures")
with(subset(tov.G2im.uniq, tov.G2im.uniq$Ct.Mean>=35.01), 
     points(y = R2, x = Ct.Mean,
            pch=20, cex = 0.75, col="blue"))
with(subset(tov.G2im.uniq, tov.G2im.uniq$Ct.Mean<35.01), 
     points(y = R2, x =Ct.Mean,
            pch=20, cex = 0.75, col="red"))
abline(v=35.0, col="black", lty=2, lwd=2)


tov.G2.3 <- subset(tov, tov$Stage == "G2 3")
tov.G2.3.uniq <- tov.G2.3[!duplicated(tov.G2.3$Sample),]

plot(y = tov.G2.3.uniq$R2, x =tov.G2.3.uniq$Ct.Mean, pch = 1,cex=0.75, col ="grey", 
     xlim=c(20,45), ylim=c(0.7,1), ylab="R2 value", xlab="Ct Mean value", main="G2 3 day")
with(subset(tov.G2.3.uniq, tov.G2.3.uniq$Ct.Mean>=35.01), 
     points(y = R2, x = Ct.Mean,pch=20, cex = 0.75, col="blue"))
with(subset(tov.G2.3.uniq, tov.G2.3.uniq$Ct.Mean<35.01), 
     points(y = R2, x = Ct.Mean,pch=20, cex = 0.75, col="red"))
abline(v=35.0, col="black", lty=2, lwd=2)


tov.G2.7 <- subset(tov, tov$Stage == "G2 7")
tov.G2.7.uniq <- tov.G2.7[!duplicated(tov.G2.7$Sample),]

plot(y = tov.G2.7.uniq$R2, x =tov.G2.7.uniq$Ct.Mean, pch = 1,cex=0.75, col ="grey", 
     xlim=c(20,45), ylim=c(0.7,1), ylab="R2 value", xlab="Ct Mean value", main="G2 7 day")
with(subset(tov.G2.7.uniq, tov.G2.7.uniq$Ct.Mean>=35.01), 
     points(y = R2, x = Ct.Mean,pch=20, cex = 0.75, col="blue"))
with(subset(tov.G2.7.uniq, tov.G2.7.uniq$Ct.Mean<35.01), 
     points(y = R2, x = Ct.Mean,pch=20, cex = 0.75, col="red"))
abline(v=35.0, col="black", lty=2, lwd=2)


tov.G2.14 <- subset(tov, tov$Stage == "G2 14")
tov.G2.14.uniq <- tov.G2.14[!duplicated(tov.G2.14$Sample),]

plot(y = tov.G2.14.uniq$R2, x =tov.G2.14.uniq$Ct.Mean, pch = 1,cex=0.75, col ="grey", 
     xlim=c(20,45), ylim=c(0.7,1), ylab="R2 value", xlab="Ct Mean value", main="G2 14 day")
with(subset(tov.G2.14.uniq, tov.G2.14.uniq$Ct.Mean>=35.01), 
     points(y = R2, x = Ct.Mean,pch=20, cex = 0.75, col="blue"))
with(subset(tov.G2.14.uniq, tov.G2.14.uniq$Ct.Mean<35.01), 
     points(y = R2, x = Ct.Mean,pch=20, cex = 0.75, col="red"))
abline(v=35.0, col="black", lty=2, lwd=2)


tov.G3.3 <- subset(tov, tov$Stage == "G3 3")
tov.G3.3.uniq <- tov.G3.3[!duplicated(tov.G3.3$Sample),]

plot(y = tov.G3.3.uniq$R2, x =tov.G3.3.uniq$Ct.Mean, pch = 1,cex=0.75, col ="grey", 
     xlim=c(20,45), ylim=c(0.7,1), ylab="R2 value", xlab="Ct Mean value", main="G3 3 day")
with(subset(tov.G3.3.uniq, tov.G3.3.uniq$Ct.Mean>=35.01), 
     points(y = R2, x = Ct.Mean,pch=20, cex = 0.75, col="blue"))
with(subset(tov.G3.3.uniq, tov.G3.3.uniq$Ct.Mean<35.01), 
     points(y = R2, x = Ct.Mean,pch=20, cex = 0.75, col="red"))
abline(v=35.0, col="black", lty=2, lwd=2)


tov.G3.7 <- subset(tov, tov$Stage == "G3 7")
tov.G3.7.uniq <- tov.G3.7[!duplicated(tov.G1.7$Sample),]

plot(y = tov.G3.7.uniq$R2, x =tov.G3.7.uniq$Ct.Mean, pch = 1,cex=0.75, col ="grey", 
     xlim=c(20,45), ylim=c(0.7,1), ylab="R2 value", xlab="Ct Mean value", main="G3 7 day")
with(subset(tov.G3.7.uniq, tov.G3.7.uniq$Ct.Mean>=35.01), 
     points(y = R2, x = Ct.Mean,pch=20, cex = 0.75, col="blue"))
with(subset(tov.G3.7.uniq, tov.G3.7.uniq$Ct.Mean<35.01), 
     points(y = R2, x = Ct.Mean,pch=20, cex = 0.75, col="red"))
abline(v=35.0, col="black", lty=2, lwd=2)


tov.G3.14 <- subset(tov, tov$Stage == "G3 14")
tov.G3.14.uniq <- tov.G3.14[!duplicated(tov.G3.14$Sample),]

plot(y = tov.G3.14.uniq$R2, x =tov.G3.14.uniq$Ct.Mean, pch = 1,cex=0.75, col ="grey", 
     xlim=c(20,45), ylim=c(0.7,1), ylab="R2 value", xlab="Ct Mean value", main="G3 14 day")
with(subset(tov.G3.14.uniq, tov.G3.14.uniq$Ct.Mean>=35.01), 
     points(y = R2, x = Ct.Mean,pch=20, cex = 0.75, col="blue"))
with(subset(tov.G3.14.uniq, tov.G3.14.uniq$Ct.Mean<35.01), 
     points(y = R2, x = Ct.Mean,pch=20, cex = 0.75, col="red"))
abline(v=35.0, col="black", lty=2, lwd=2)

#dummy var for cycles

tov$G1 <- ifelse(tov$Cycle == "G1", 0,0)
tov$G2 <- ifelse(tov$Cycle == "G2", 1,0)
tov$G3 <- ifelse(tov$Cycle == "G3", 1,0)

tov$tw7 <- ifelse(tov$Temp == 27, 0,0)
tov$th0 <- ifelse(tov$Temp == 30, 1,0)
tov$th3 <- ifelse(tov$Temp == 33, 1,0)

tov$im <- ifelse(tov$Day == 0, 0,0)
tov$three <- ifelse(tov$Day == 3, 1,0)
tov$seven <- ifelse(tov$Day == 7, 1,0)
tov$fourteen <- ifelse(tov$Day == 14, 1,0)

#logit regress + poisson

logit.tov <-
  glm(pos ~  three + seven + fourteen, 
      data = tov.uniq , family = "poisson", na.action = na.exclude)
summary(logit.tov)
