############################################################################
############################################################################
#Species co-occurrence analysis in the Grand River, Ontario#################
#Code prepared by Dr. Karl A. Lamothe, Fisheries and Oceans Canada, GLLFAS## 
#karl.lamothe@dfo-mpo.gc.ca#################################################
############################################################################
############################################################################
library(pacman)
p_load(psych)
p_load(vegan)
p_load(RPresence)
p_load(ggplot2)
p_load(RCurl)

#Set working directory
setwd("...")

#Load in data sets
#community presence absence data
Comm.data <- read.csv("Multivariate Analyses/Fish.sites.comm.abundance.csv",header=T) 
head(Comm.data)
colSums(Comm.data) #check to ensure data looks ok

#habitat covariates
site.habitat <- read.csv("Multivariate Analyses/Community.Site.Cov.csv",header=T)
head(site.habitat)

#########################################################################################################
#########################################################################################################
#Multivariate analysis
#Gather habitat variables to work with
habitat<-as.data.frame(cbind(Mean.Substrate = site.habitat$Mean.substrate..mm.,
                                 Sand.FineG= site.habitat$Sand.Finegravel..proportion.,
                                 Turbidity = site.habitat$Turbidity..cm.,
                                 Velocity = site.habitat$Flow.06D..m.s.,
                                 Depth = site.habitat$Mean.depth..cm.))

#Perform correlation analysis
x<-corr.test(habitat, method="pearson")
y<-corr.test(habitat, method="spearman")

#Look at output
round(x$r,2) #corr values
round(x$p,2) #p values
round(y$r,2) #corr values
round(y$p,2) #p values

#plot pairwise relationship between habitat variables
pairs(habitat,pch=20) 

#Scale and center the variables
habitat.scaled.1<-data.frame(scale(habitat[,c(1,3,5)],scale=T,center=T))
head(habitat.scaled.1)
habitat.scaled<-cbind(habitat.scaled.1,habitat$Velocity, habitat$Sand.FineG)
colnames(habitat.scaled)<-c("Mean.Sub","Turbi","Depth","Flow","sfg")
head(habitat.scaled)

#PCA
par(mfrow=c(1,1))
site.sp.trans<-decostand(Comm.data, 'hellinger') #hellinger transformation
PCA.Species<-rda(site.sp.trans, cor=T, scale = T) #performs PCA
biplot(PCA.Species, scaling=1, pch=20, display = "species") #plot
summary(PCA.Species)

#fuction to assess significance of the principal components.
sign.pc<-function(x,R=9999,s=10, cor=T,...){
  pc.out<-princomp(x,cor=cor,...)  # run PCA
  pve=(pc.out$sdev^2/sum(pc.out$sdev^2))[1:s]  # the proportion of variance of each PC
  pve.perm<-matrix(NA,ncol=s,nrow=R)  # a matrix with R rows and s columns
  for(i in 1:R){
    x.perm<-apply(x,2,sample)# permutation each column
    pc.perm.out<-princomp(x.perm,cor=cor,...)# run PCA
    pve.perm[i,]=(pc.perm.out$sdev^2/sum(pc.perm.out$sdev^2))[1:s] # the proportion of variance of each PC.perm
  }
  pval<-apply(t(pve.perm)>pve,1,sum)/R # calcalute the p-values
  return(list(pve=pve,pval=pval))
}

# apply the function
sign.pc(site.sp.trans,cor=F)

#envfit
ef<-envfit(PCA.Species, habitat.scaled, na.rm=T, permutations=9999, choices=c(1:4))
ef

#plot
par(mfrow=c(3,1))
summary(PCA.Species)
species<-scores(PCA.Species, choices=c(1:4))
species.df<-as.data.frame(species$species)

plot(species.df$PC1,species.df$PC2, xlab="Component 1 (10.24%)", ylab="Component 2 (6.01%)",
     pch=20,las=1,asp=1,type="n")
abline(v=0, h=0, lty="dashed")
text(species.df$PC1,species.df$PC2, labels=row.names(species.df))
plot(ef, p.max=0.05, choices=c(1,2), col = "red")

plot(species.df$PC2,species.df$PC3, xlab="Component 2 (6.01%)", ylab="Component 3 (4.97%)",
     pch=20,las=1,asp=1,type="n")
abline(v=0, h=0, lty="dashed")
text(species.df$PC2,species.df$PC3, labels=row.names(species.df))
plot(ef, p.max=0.05, choices=c(2,3), col = "red")

plot(species.df$PC3,species.df$PC4, xlab="Component 3 (4.97%)", ylab="Component 4 (4.46%)",
     pch=20,las=1,asp=1,type="n")
abline(v=0, h=0, lty="dashed")
text(species.df$PC3,species.df$PC4, labels=row.names(species.df))
plot(ef, p.max=0.05, choices=c(3,4), col = "red")

plot(PCA.Species$x[,1], PCA.Species$x[,2], xlab=paste("PCA 1 (", round(s$importance[2]*100, 1), "%)", sep = ""), 
     ylab=paste("PCA 2 (", round(s$importance[5]*100, 1), "%)", sep = ""), 
     pch=pch.group, col="black", bg=col.group, cex=2, las=1, asp=1)

par(mai=c(0.5,0.5,0.1,0.1))
biplot(PCA.Species, scaling=-1, pch=20, display = c("species"), col="black",  las=1)
plot(ef, p.max=0.05, choices=c(1,2), col = "red")
biplot(PCA.Species, scaling=-1, pch=20, display = c("species"), las=1,choices = c(2:3), col="black")
plot(ef, p.max=0.05, choices=c(2,3), col="red")
biplot(PCA.Species, scaling=-1, pch=20, display = c("species"), choices = c(3:4), col="black",  las=1)
plot(ef, p.max=0.05, choices=c(3,4), col = "red")

#########################################################################
#Occupancy Models
rm(list=ls())
{setwd("C:/Users/LAMOTHEK/Desktop/Lamothe/Eastern Sand Darter/Species Cooccurrence/Data/Occupancy Models")
  data<-read.csv("ESD.SS.Grand.csv",header=T) #Species occurrence data
  unit.cov<-read.csv("Covariate.Grand.csv",header=T) #Covariate
  surv.cov<-read.csv("Removal.Grand.csv",header=T)} #Removal
head(data)
str(data)
data<-data[,-5]

#Save non transformed variables for each site for plotting
Turbidity.non<-unit.cov$Turbidity..cm.[1:151]
Depth.non<-unit.cov$Mean.depth..cm.[1:151]
Sfg.non<-unit.cov$Sand.Finegravel..proportion.[1:151]
Subs.non<-unit.cov$Mean.substrate..mm.[1:151]
Velocity.non<-unit.cov$Flow.06D..m.s.[1:151]

#Prepare site specific data for occ models
head(unit.cov)
unit.cov1<-data.frame(scale(unit.cov[,c(2,4,12)],scale=T,center=T))
head(unit.cov1)
unit.cov1<-cbind(unit.cov1,unit.cov$Flow.06D..m.s., unit.cov$Sand.Finegravel..proportion.)
colnames(unit.cov1)<-c("Depth","Mean.Sub","Turbi","Flow","sfg")
head(unit.cov1)

#Prepare survey specific data for occ models
head(surv.cov)
surv.cov<-surv.cov[,-c(1,5)]
surv.cov1<-as.data.frame(t(surv.cov))
surv.cov2<-stack(surv.cov1)
surv.cov2<-data.frame(surv.cov2[,1])
colnames(surv.cov2)<-c("Rem")
head(surv.cov2)

#check to make sure the variable is formatted properly
length(surv.cov2$Rem)==nrow(data[,-1])*ncol(data[,-1])

#Create data set up for occupancy model
WT<-createPao(data=data[,-1], unitcov=unit.cov1, survcov=surv.cov2)

#Detection models
m1<-occMod(model=list(psi~1,p~1),data=WT,type="so.2sp.2") #intercept model
m2<-occMod(model=list(psi~1,p~SP),data=WT,type="so.2sp.2") #SP = species level effect, occupancy varies by species
m3<-occMod(model=list(psi~1,p~Rem),data=WT,type="so.2sp.2")
m4<-occMod(model=list(psi~1,p~Mean.Sub),data=WT,type="so.2sp.2")
m5<-occMod(model=list(psi~1,p~SP+Rem),data=WT,type="so.2sp.2")
m6<-occMod(model=list(psi~1,p~SP+Mean.Sub),data=WT,type="so.2sp.2")
m7<-occMod(model=list(psi~1,p~Rem+Mean.Sub),data=WT,type="so.2sp.2")
m8<-occMod(model=list(psi~1,p~SP*Rem),data=WT,type="so.2sp.2")
m9<-occMod(model=list(psi~1,p~SP*Mean.Sub),data=WT,type="so.2sp.2")

results<-createAicTable(list(m1,m2,m3,m4,m5,m6,m7,m8,m9))
summary(results)
#m9 is the best model for detection which includes mean substrate and species specific detection

summary(m1);summary(m2);summary(m3);summary(m4)
summary(m5);summary(m6);summary(m7);summary(m8)
summary(m9)

mean(modAvg(results, param = "pA")$est)
mean(modAvg(results, param = "pB")$est)
sd(modAvg(results, param = "pA")$est)/sqrt(151)
sd(modAvg(results, param = "pB")$est)/sqrt(151)

#Occupancy
m10<-occMod(model=list(psi~SP,p~SP*Mean.Sub),data=WT,type="so.2sp.2")#Species specific occupancy

m11<-occMod(model=list(psi~SP+Turbi,p~SP*Mean.Sub),data=WT,type="so.2sp.2")#Species specific occupancy + turbidity
m12<-occMod(model=list(psi~SP+Depth,p~SP*Mean.Sub),data=WT,type="so.2sp.2")#Species specific occupancy + depth
m13<-occMod(model=list(psi~SP+sfg,p~SP*Mean.Sub),data=WT,type="so.2sp.2")#Species specific occupancy + sfg
m14<-occMod(model=list(psi~SP+Flow,p~SP*Mean.Sub),data=WT,type="so.2sp.2")

m15<-occMod(model=list(psi~SP*Turbi,p~SP*Mean.Sub),data=WT,type="so.2sp.2")#Species specific occupancy*turbidity
m16<-occMod(model=list(psi~SP*Depth,p~SP*Mean.Sub),data=WT,type="so.2sp.2")#Species specific occupancy*depth
m17<-occMod(model=list(psi~SP*sfg,p~SP*Mean.Sub),data=WT,type="so.2sp.2")#Species specific occupancy*sfg
m18<-occMod(model=list(psi~SP*Flow,p~SP*Mean.Sub),data=WT,type="so.2sp.2")#Species specific occupancy*sfg

m19<-occMod(model=list(psi~SP+INT+Turbi,p~SP*Mean.Sub),data=WT,type="so.2sp.2")
m19$warnings
coef(m19, "psi")
unique(fitted(m19, "psiA"))
unique(fitted(m19, "psiBa"))
unique(fitted(m19, "psiBA"))
unique(fitted(m19, "nu"))#this is why we have a warning - unable to estimate nu, otherwise model is fine
unique(fitted(m19, "pA"))
unique(fitted(m19, "pB"))

m19$dmat

m20<-occMod(model=list(psi~SP+INT+Depth,p~SP*Mean.Sub),data=WT,type="so.2sp.2")
m20$warnings
unique(fitted(m20, "nu"))#this is why we have a warning - unable to estimate nu, otherwise model is fine

#m21<-occMod(model=list(psi~sfg+SP+INT,p~SP*Mean.Sub),data=WT,type="so.2sp.2") 
#m21$warnings#convergence warning > 3, likely not an issue but try random starting values
#set.seed(29534)
#m21<-occMod(model=list(psi~sfg+SP+INT,p~SP*Mean.Sub),data=WT,type="so.2sp.2", randinit = 500)
#m21$warnings
#coef(m21, "psi") #clearly the interaction term is causing the issue
#coef(m21,"p")
#m22<-occMod(model=list(psi~Flow+SP+INT,p~SP*Mean.Sub),data=WT,type="so.2sp.2") 

m23<-occMod(model=list(psi~INT+Turbi*SP,p~SP*Mean.Sub),data=WT,type="so.2sp.2") 
m24<-occMod(model=list(psi~INT+Depth*SP,p~SP*Mean.Sub),data=WT,type="so.2sp.2")
m25<-occMod(model=list(psi~INT+sfg*SP,p~SP*Mean.Sub),data=WT,type="so.2sp.2") 
m25$warnings

#m26<-occMod(model=list(psi~INT+Flow*SP,p~SP*Mean.Sub),data=WT,type="so.2sp.2")
#m26$warnings
#coef(m26, "psi")
#unique(fitted(m26, "psiA"))
#unique(fitted(m26, "psiBa"))
#unique(fitted(m26, "psiBA"))
#unique(fitted(m26, "nu"))#this is why we have a warning - unable to estimate nu, otherwise model is fine
#unique(fitted(m26, "pA"))
#unique(fitted(m26, "pB"))

results<-createAicTable(list(m1,m2,m3,m4,m5,m6,m7,m8,m9,
                             m10,m11,m12,m13,m14,
                             m15,m16,m17,m18,m19,
                             m20,#m21,m22,
                             m23,m24,m25))
summary(results)

summary(m1);summary(m2);summary(m3);summary(m4)
summary(m5);summary(m6);summary(m7);summary(m8)
summary(m9);summary(m10);summary(m11);summary(m12);summary(m13)
summary(m14);summary(m15);summary(m16);summary(m17);summary(m18)
summary(m19);summary(m20);summary(m21);summary(m22);summary(m23);
summary(m24);summary(m25);summary(m26)

coef(m25,"psi")
coef(m25,"p")
coef(m25,"psi.VC")
coef(m25,"p.VC")

unique(fitted(m25, "psiA"))
unique(fitted(m25, "psiBa"))
unique(fitted(m25, "psiBA"))
unique(fitted(m25, "nu"))
unique(fitted(m25, "pA"))
unique(fitted(m25, "pB"))

mean(modAvg(results, param = "psiA")$est)
mean(modAvg(results, param = "psiBa")$est)
mean(modAvg(results, param = "psiBA")$est)
sd(modAvg(results, param = "psiA")$est)/sqrt(151)
sd(modAvg(results, param = "psiBa")$est)/sqrt(151)
sd(modAvg(results, param = "psiBA")$est)/sqrt(151)

mean(modAvg(results, param = "psiA")$est)*mean(modAvg(results, param = "psiBA")$est)+
       (1-mean(modAvg(results, param = "psiA")$est))*mean(modAvg(results, param = "psiBa")$est)
sd(modAvg(results, param = "psiA")$est)*sd(modAvg(results, param = "psiBA")$est)+
       (1-sd(modAvg(results, param = "psiA")$est))*sd(modAvg(results, param = "psiBa")$est)

0.1003252/sqrt(151)

#SIF
(mean(modAvg(results, param = "psiA")$est)*mean(modAvg(results, param = "psiBA")$est))/
  (mean(modAvg(results, param = "psiA")$est)*(mean(modAvg(results, param = "psiA")$est)*
     mean(modAvg(results, param = "psiBA")$est)+(1-mean(modAvg(results, param = "psiA")$est))*
        mean(modAvg(results, param = "psiBa")$est)))
   
(mean(fitted(m25, param = "psiA")$est)*mean(fitted(m25, param = "psiBA")$est))/
  (mean(fitted(m25, param = "psiA")$est)*(mean(fitted(m25, param = "psiA")$est)*
                                                mean(fitted(m25, param = "psiBA")$est)+(1-mean(fitted(m25, param = "psiA")$est))*
                                                mean(fitted(m25, param = "psiBa")$est)))
(sd(fitted(m25, param = "psiA")$est)*sd(fitted(m25, param = "psiBA")$est))/
  (sd(fitted(m25, param = "psiA")$est)*(sd(fitted(m25, param = "psiA")$est)*
                                            sd(fitted(m25, param = "psiBA")$est)+(1-sd(fitted(m25, param = "psiA")$est))*
                                            sd(fitted(m25, param = "psiBa")$est)))
1.273088/sqrt(151)

(mean(fitted(m26, param = "psiA")$est)*mean(fitted(m26, param = "psiBA")$est))/
  (mean(fitted(m26, param = "psiA")$est)*(mean(fitted(m26, param = "psiA")$est)*
                                            mean(fitted(m26, param = "psiBA")$est)+(1-mean(fitted(m26, param = "psiA")$est))*
                                            mean(fitted(m26, param = "psiBa")$est)))
(mean(fitted(m19, param = "psiA")$est)*mean(fitted(m19, param = "psiBA")$est))/
  (mean(fitted(m19, param = "psiA")$est)*(mean(fitted(m19, param = "psiA")$est)*
                                            mean(fitted(m19, param = "psiBA")$est)+(1-mean(fitted(m19, param = "psiA")$est))*
                                            mean(fitted(m19, param = "psiBa")$est)))
(mean(fitted(m23, param = "psiA")$est)*mean(fitted(m23, param = "psiBA")$est))/
  (mean(fitted(m23, param = "psiA")$est)*(mean(fitted(m23, param = "psiA")$est)*
                                            mean(fitted(m23, param = "psiBA")$est)+(1-mean(fitted(m23, param = "psiA")$est))*
                                            mean(fitted(m23, param = "psiBa")$est)))
(mean(fitted(m20, param = "psiA")$est)*mean(fitted(m20, param = "psiBA")$est))/
  (mean(fitted(m20, param = "psiA")$est)*(mean(fitted(m20, param = "psiA")$est)*
                                            mean(fitted(m20, param = "psiBA")$est)+(1-mean(fitted(m20, param = "psiA")$est))*
                                            mean(fitted(m20, param = "psiBa")$est)))
(mean(fitted(m24, param = "psiA")$est)*mean(fitted(m24, param = "psiBA")$est))/
  (mean(fitted(m24, param = "psiA")$est)*(mean(fitted(m24, param = "psiA")$est)*
                                            mean(fitted(m24, param = "psiBA")$est)+(1-mean(fitted(m24, param = "psiA")$est))*
                                            mean(fitted(m24, param = "psiBa")$est)))
SIF<-c(0.104,0,0)
weight<-c(0.894,0.058,0.047)
weighted.mean(SIF,weight)
###############################################################
##################################
#plotting functions
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  require(grid)
  plots <- c(list(...), plotlist)
  numPlots = length(plots)
  if (is.null(layout)) {
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  if (numPlots==1) {
    print(plots[[1]])
  } else {
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    for (i in 1:numPlots) {
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

theme_me <- theme_bw() +
  theme(axis.title = element_text(size = 11, family = "sans", colour = "black"),
        axis.text.x = element_text(size = 10, family = "sans", colour = "black"),
        axis.text.y = element_text(size = 10, family = "sans", colour = "black"),
        legend.title = element_text(size = 10, family = "sans", colour = "black"),
        legend.text = element_text(size = 8, family = "sans", colour = "black"),
        strip.text = element_text(size = 11, family = "sans", colour = "black"))

############################Figure 3#############################################
ma.psiA<-modAvg(results, param = "psiA",index=1)
ma.psiBa<-modAvg(results, param = "psiBa",index=1)
ma.psiBA<-modAvg(results, param = "psiBA",index=1)

ma.pA<-modAvg(results, param = "pA",index=1)
ma.pB<-modAvg(results, param = "pB",index=1)
ma.prA<-modAvg(results, param = "rA",index=1)
ma.prBa<-modAvg(results, param = "rBa",index=1)
ma.prBA<-modAvg(results, param = "rBA",index=1)

mean((ma.psiA$est*ma.psiBA$est)/(ma.psiA$est*((ma.psiA$est*ma.psiBA$est)+(1-ma.psiA$est)*ma.psiBa$est)))
sd((ma.psiA$est*ma.psiBA$est)/(ma.psiA$est*((ma.psiA$est*ma.psiBA$est)+(1-ma.psiA$est)*ma.psiBa$est)))

idx<-order(Sfg.non)
w<-as.data.frame(cbind(Sfg.non[idx],ma.psiA$est[idx]))
x<-as.data.frame(cbind(Sfg.non[idx],ma.psiBa$est[idx]))
y<-as.data.frame(cbind(Sfg.non[idx],ma.psiBA$est[idx]))
z<-as.data.frame(cbind(Sfg.non[idx],ma.psiA$est[idx]*ma.psiBA$est[idx]+(1-ma.psiA$est[idx])*ma.psiBa$est[idx]))

Sp<-c(rep("psiA",151),rep("psiBa", 151),rep("psiBA",151),rep("psiB",151))
colnames(w)<-c("Sfg.non","Occupancy")
colnames(x)<-c("Sfg.non","Occupancy")
colnames(y)<-c("Sfg.non","Occupancy")
colnames(z)<-c("Sfg.non","Occupancy")
wxyz<-rbind(w,x,y,z)
wxyz<-cbind(wxyz,Sp)
head(wxyz)

sfg.gg<-ggplot()+
  geom_line(data=wxyz, aes(x=Sfg.non,y=Occupancy,color=factor(Sp)),lwd=1)+
  ylim(0,1)+
  ylab("Probability of occupancy")+
  xlab("Proportion of sand and fine gravel")+
  scale_color_manual(values=c("#000000", "#E69F00", "#56B4E9", "#009E73"),
                     name="",
                     breaks=c("psiA", "psiBa","psiBA","psiB"),
                     labels=c("ESD","SS | ESD = 0","SS | ESD = 1", "SS"))+
  theme_me+
  theme(legend.position="none")+
  annotate("text",x=0,y=1,label="A",size=4)
sfg.gg

#Water clarity
summary(results)
ma.psiA<-modAvg(results, param = "psiA",index=6)
ma.psiBa<-modAvg(results, param = "psiBa",index=6)
ma.psiBA<-modAvg(results, param = "psiBA",index=6)

idx<-order(Turbidity.non)
w<-as.data.frame(cbind(Turbidity.non[idx],ma.psiA$est[idx]))
x<-as.data.frame(cbind(Turbidity.non[idx],ma.psiBa$est[idx]))
y<-as.data.frame(cbind(Turbidity.non[idx],ma.psiBA$est[idx]))
z<-as.data.frame(cbind(Turbidity.non[idx],ma.psiA$est[idx]*ma.psiBA$est[idx]+(1-ma.psiA$est[idx])*ma.psiBa$est[idx]))

Sp<-c(rep("psiA",151),rep("psiBa", 151),rep("psiBA",151),rep("psiB",151))
colnames(w)<-c("Turbidity.non","Occupancy")
colnames(x)<-c("Turbidity.non","Occupancy")
colnames(y)<-c("Turbidity.non","Occupancy")
colnames(z)<-c("Turbidity.non","Occupancy")
wxyz<-rbind(w,x,y,z)
wxyz<-cbind(wxyz,Sp)
head(wxyz)

turbidgg<-ggplot()+
  geom_line(data=wxyz, aes(x=Turbidity.non,y=Occupancy,color=factor(Sp)),lwd=1)+
  ylim(0,1)+
  ylab("Probability of occupancy")+
  xlab("Water clarity (cm)")+
  scale_color_manual(values=c("#000000", "#E69F00", "#56B4E9", "#009E73"),
                     name="",
                     breaks=c("psiA", "psiBa","psiBA","psiB"),
                     labels=c("ESD","SS | ESD = 0","SS | ESD = 1", "SS"))+
  theme_me+
  theme(legend.position="none")+
  annotate("text",x=22,y=1,label="B",size=4)
turbidgg

#depth
summary(results)
ma.psiA<-modAvg(results, param = "psiA",index=10)
ma.psiBa<-modAvg(results, param = "psiBa",index=10)
ma.psiBA<-modAvg(results, param = "psiBA",index=10)

idx<-order(Depth.non)
w<-as.data.frame(cbind(Depth.non[idx],ma.psiA$est[idx]))
x<-as.data.frame(cbind(Depth.non[idx],ma.psiBa$est[idx]))
y<-as.data.frame(cbind(Depth.non[idx],ma.psiBA$est[idx]))
z<-as.data.frame(cbind(Depth.non[idx],ma.psiA$est[idx]*ma.psiBA$est[idx]+(1-ma.psiA$est[idx])*ma.psiBa$est[idx]))

Sp<-c(rep("psiA",151),rep("psiBa", 151),rep("psiBA",151),rep("psiB",151))
colnames(w)<-c("Depth.non","Occupancy")
colnames(x)<-c("Depth.non","Occupancy")
colnames(y)<-c("Depth.non","Occupancy")
colnames(z)<-c("Depth.non","Occupancy")
wxyz<-rbind(w,x,y,z)
wxyz<-cbind(wxyz,Sp)
head(wxyz)

depthgg<-ggplot()+
  geom_line(data=wxyz, aes(x=Depth.non,y=Occupancy,color=factor(Sp)),lwd=1)+
  ylim(0,1)+
  ylab("Probability of occupancy")+
  xlab("Depth (cm)")+
  scale_color_manual(values=c("#000000", "#E69F00", "#56B4E9", "#009E73"),
                     name="", breaks=c("psiA", "psiBa","psiBA","psiB"),
                     labels=c("ESD","SS | ESD = 0","SS | ESD = 1", "SS"))+
  theme_me+
  theme(legend.position="none")+
  annotate("text",x=20,y=1,label="C",size=4)
depthgg

#velocity
summary(results)
ma.psiA<-modAvg(results, param = "psiA",index=2)
ma.psiBa<-modAvg(results, param = "psiBa",index=2)
ma.psiBA<-modAvg(results, param = "psiBA",index=2)

idx<-order(Velocity.non)
w<-as.data.frame(cbind(Velocity.non[idx],ma.psiA$est[idx]))
x<-as.data.frame(cbind(Velocity.non[idx],ma.psiBa$est[idx]))
y<-as.data.frame(cbind(Velocity.non[idx],ma.psiBA$est[idx]))
z<-as.data.frame(cbind(Velocity.non[idx],ma.psiA$est[idx]*ma.psiBA$est[idx]+(1-ma.psiA$est[idx])*ma.psiBa$est[idx]))

Sp<-c(rep("psiA",151),rep("psiBa", 151),rep("psiBA",151),rep("psiB",151))
colnames(w)<-c("Velocity.non","Occupancy")
colnames(x)<-c("Velocity.non","Occupancy")
colnames(y)<-c("Velocity.non","Occupancy")
colnames(z)<-c("Velocity.non","Occupancy")
wxyz<-rbind(w,x,y,z)
wxyz<-cbind(wxyz,Sp)
head(wxyz)

Velocitygg<-ggplot()+
  geom_line(data=wxyz, aes(x=Velocity.non,y=Occupancy,color=factor(Sp)),lwd=1)+
  ylim(0,1)+
  ylab("Probability of occupancy")+
  xlab("Velocity (m/s)")+
  scale_color_manual(values=c("#000000", "#E69F00", "#56B4E9", "#009E73"),
                     name="", breaks=c("psiA", "psiBa","psiBA","psiB"),
                     labels=c("ESD","SS | ESD = 0","SS | ESD = 1", "SS"))+
  theme_me+
  theme(legend.position="none")+
  annotate("text",x=0.75,y=1,label="D",size=4)
Velocitygg

multiplot(sfg.gg,turbidgg,depthgg,Velocitygg,cols=2)

###################################
#Figure 2
ma.pA<-modAvg(results, param = "pA",index=1)
ma.prBa<-modAvg(results, param = "rBa",index=1)

idx<-order(Subs.non)
w<-as.data.frame(cbind(Subs.non[idx],ma.pA$est[idx]))
z<-as.data.frame(cbind(Subs.non[idx],ma.prBa$est[idx]))
Sp<-c(rep("Eastern sand darter",151),rep("Silver shiner", 151))
colnames(w)<-c("Subs.non","Detection")
colnames(z)<-c("Subs.non","Detection")
wz<-rbind(w,z)
wz<-cbind(wz,Sp)
head(wz)

detection<-ggplot()+
  geom_line(data=wz, aes(x=Subs.non,y=Detection,lty=factor(Sp)),lwd=1)+
  
  ylim(0,1)+
  ylab("Probability of detection")+
  xlab("Mean substrate size (mm)")+
  scale_linetype_discrete(name="",
                          breaks=c("Eastern sand darter", "Silver shiner"),
                          labels=c(expression(italic("Ammocrypta pellucida"),italic("Notropis photogenis"))))+
  theme_me+
  theme(legend.position=c(.6,.8),
        legend.text = element_text(size=12),
        legend.background=element_blank())

detection

