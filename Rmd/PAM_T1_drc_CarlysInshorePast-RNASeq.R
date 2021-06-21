setwd("/Users/danbarshis/dansstuff/Projeks/ODU/Projeks/Field/2021_Keys/StressExperiments/2021-06-08_IPast-RNASeq-Inshore/scripts")
library(Hmisc)
library(lsmeans)
library(emmeans)
library(car)
library(drc)

Allpam<-read.delim("../data/PAM/2021-06-08_CarlysInshore_Past_RNASeq.txt")

Allpam$Geno<-as.factor(Allpam$Geno)
aggregate(T3PAM ~ Temp, data=Allpam, summary)

Allpam$FacTemp<-as.factor(Allpam$Temp)
aggregate(T3PAM ~ Temp, data=Allpam, length)

bartlett.test(T3PAM ~ FacTemp, data=Allpam)
leveneTest(T3PAM ~ FacTemp, data=Allpam)
aggregate(T3PAM ~ FacTemp, data=Allpam, FUN= function(x) shapiro.test(x)$p.value)

n1<-aov(T3PAM ~ FacTemp + Error(Geno/FacTemp), Allpam)
summary(n1)

print(lsmeans(n1, list(pairwise ~ FacTemp)), adjust = c("tukey"))

#### Full curve fit ####
DRCpam = drm(T3PAM ~ Temp, data = Allpam,
                 fct = LL.3(names = c('hill', 'max', 'ed50')))
summary(DRCpam)
plot(DRCpam)

DRCpamgeno = drm(T3PAM ~ Temp, data = Allpam, curveid = Geno,
             fct = LL.3(names = c('hill', 'max', 'ed50')))
summary(DRCpamgeno)
compParm(DRCpamgeno, 'ed50')
compParm(DRCpamgeno, 'ed50', "-")
plot(DRCpamgeno)
ED(DRCpamgeno, c(50))[,1]

#### Plotting ####
temp_x<- seq(30, 40.5, length = 100)

pdf("../plots/2021-06-08_CarlysInshorePastRNASeq_DRC.pdf",10,7)
line_width=2
Colorz="blue"
Syms=16
Denscity <- 45
matplot(temp_x, predict(DRCpamgeno, data.frame(Temp = temp_x), interval="confidence"),
        type="n",col=Colorz,lty=c(1,3,3),lwd=line_width,ylab="Fv/Fm",xlab="Temperature °C", xlim=c(29.5,41),ylim=c(0,0.65), cex.axis=1.5, cex.lab=1.5)
polygon(c(temp_x, rev(temp_x)),c(predict(DRCpamgeno, data.frame(Temp = temp_x), interval="confidence")[,2],rev(predict(DRCpamgeno, data.frame(Temp = temp_x), interval="confidence")[,3])), col=Colorz, density = Denscity)
matpoints(temp_x, predict(DRCpamgeno, data.frame(Temp = temp_x), interval="confidence"),
          type="l",col=Colorz,lty=c(1,3,3),lwd=line_width)
with(Allpam,matpoints(Temp,T3PAM,pch=Syms, col=Colorz, cex=1.5))

title(main="Carly's Inshore Porites astreoides")
abline(v=mean(ED(DRCpamgeno, c(50))[,1]), col=Colorz)
text(mean(ED(DRCpamgeno, c(50))[,1])+0.6, 0.65,labels=as.character(round(mean(ED(DRCpamgeno, c(50))[,1]), digits=2)),col=Colorz, cex=1.5)
dev.off()

pdf("../plots/2021-06-08_CarlysInshorePastRNASeq_DRC_genofits.pdf",10,7)
par(mar=c(5.1,5.1,4.1,2.1))
Line_width=2.5
with(Allpam, plot(Temp,T3PAM,xlim=c(29.5,41),ylim=c(0,0.60), type='n', cex.axis=2.5, cex.lab=2.5, ylab="Fv/Fm", xlab="Temperature °C"))
Colours<-c('#4575b4','#74add1','#fee090','#f46d43','#d73027')
for(i in order(DRCpamgeno$coefficients[11:15])){
  Geno=names(table(Allpam$Geno))[i]
  print(Geno)
  print(summary(mod1<-drm(T3PAM ~ Temp, data=Allpam[Allpam$Geno==Geno,],fct = LL.3(names = c('hill', 'max', 'ed50')))))
  points(Allpam[Allpam$Geno==Geno,"Temp"],Allpam[Allpam$Geno==Geno,"T3PAM"], col=Colours[i], pch=16)
  lines(temp_x, predict(mod1, data.frame(Temp = temp_x)), col=Colours[i], lwd=Line_width)
}
abline(v=DRCpamgeno$coefficients[11:15][order(DRCpamgeno$coefficients[11:15])], col=Colours, lwd=Line_width)
legend("bottomleft", col=Colours, lty=1, cex=1.75, lwd=Line_width, paste0(names(DRCpamgeno$coefficients)[11:15][order(DRCpamgeno$coefficients[11:15])],": ",round(DRCpamgeno$coefficients[11:15][order(DRCpamgeno$coefficients[11:15])], digits=2)))
dev.off()


