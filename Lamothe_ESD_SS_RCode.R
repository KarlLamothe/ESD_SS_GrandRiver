############################################################################
############################################################################
#Community and species co-occurrence analysis in the Grand River, Ontario###
#Code prepared by Dr. Karl A. Lamothe, Fisheries and Oceans Canada, GLLFAS## 
#karl.lamothe@dfo-mpo.gc.ca#################################################
############################################################################
############################################################################
rm(list=ls())
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
#There are 150 sites, as opposed to 151, because at 1 site there were no 
#fish captured
Comm.data <- read.csv("Fish.sites.comm.PA_Lamothe.csv",header=T) 
head(Comm.data)
colSums(Comm.data[3:ncol(Comm.data)])

#habitat covariates
site.habitat<-read.csv("Covariates_Comm_Analysis_Lamothe.csv",
                         header=T)
head(site.habitat)

############################################################################
############################################################################
#Multivariate analysis
############################################################################
############################################################################
#Gather habitat variables to work with
habitat<-as.data.frame(
  cbind(Mean.Substrate=site.habitat$Mean.substrate..mm.,
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
habitat.scaled<-cbind(habitat.scaled.1,habitat$Velocity, habitat$Sand.FineG)
colnames(habitat.scaled)<-c("Mean.Sub","Turbi","Depth","Flow","sfg")
head(habitat.scaled)

#PCA
#hellinger transformation
site.sp.trans<-decostand(Comm.data[3:ncol(Comm.data)], 'hellinger') 
PCA.Species<-rda(site.sp.trans, scale = F) #performs PCA
summary(PCA.Species)

#fuction to assess significance of the principal components.
sign.pc<-function(x,R=9999,s=10, cor=T,...){
  pc.out<-princomp(x,cor=cor,...)  # run PCA
  # the proportion of variance of each PC
  pve=(pc.out$sdev^2/sum(pc.out$sdev^2))[1:s]  
  # a matrix with R rows and s columns
  pve.perm<-matrix(NA,ncol=s,nrow=R)  
  for(i in 1:R){
    x.perm<-apply(x,2,sample)# permutation each column
    pc.perm.out<-princomp(x.perm,cor=cor,...)# run PCA
    # the proportion of variance of each PC.perm
    pve.perm[i,]=(pc.perm.out$sdev^2/sum(pc.perm.out$sdev^2))[1:s] 
  }
  pval<-apply(t(pve.perm)>pve,1,sum)/R # calcalute the p-values
  return(list(pve=pve,pval=pval))
}

# apply the function
sign.pc(site.sp.trans,cor=F)

#envfit
ef<-envfit(PCA.Species, habitat.scaled, na.rm=T, permutations=9999, 
           choices=c(1:3))
ef

#plot
biplot(PCA.Species, scaling=-1, pch=20, display = c("species"), col="black")
plot(ef, p.max=0.05, choices=c(1,2), col = "red")
biplot(PCA.Species, scaling=-1, pch=20, display = c("species"), 
       choices = c(2:3), col="black",xlim=c(-1,1))
plot(ef, p.max=0.05, choices=c(2,3), col="red", xlim=c(-1,1))

#########################################################################
#########################################################################
#Occupancy Models
#########################################################################
#########################################################################
rm(list=ls())
{#Species occurrence data
  data<-read.csv("ESD_SS_Hauls_Lamothe.csv",header=T) 
  #Covariate
  unit.cov<-read.csv("Covariate_Modelling_Lamothe.csv",header=T)
  #Removal
  surv.cov<-read.csv("Removal_Lamothe.csv",header=T)}
head(data)
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
unit.cov1<-cbind(unit.cov1,unit.cov$Flow.06D..m.s., 
                 unit.cov$Sand.Finegravel..proportion.)
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
length(surv.cov2$Rem)
nrow(data[,-1])*ncol(data[,-1])

#Create data set up for occupancy model
WT<-createPao(data=data[,-1], unitcov=unit.cov1, survcov=surv.cov2)

#Detection models
#intercept model
m1<-occMod(model=list(psi~1,p~1),data=WT,type="so.2sp.2") 
#SP = species level effect, occupancy varies by species
m2<-occMod(model=list(psi~1,p~SP),data=WT,type="so.2sp.2") 
m3<-occMod(model=list(psi~1,p~Rem),data=WT,type="so.2sp.2")
m4<-occMod(model=list(psi~1,p~Mean.Sub),data=WT,type="so.2sp.2")
m5<-occMod(model=list(psi~1,p~SP+Rem),data=WT,type="so.2sp.2")
m6<-occMod(model=list(psi~1,p~SP+Mean.Sub),data=WT,type="so.2sp.2")
m7<-occMod(model=list(psi~1,p~Rem+Mean.Sub),data=WT,type="so.2sp.2")
m8<-occMod(model=list(psi~1,p~SP*Rem),data=WT,type="so.2sp.2")
m9<-occMod(model=list(psi~1,p~SP*Mean.Sub),data=WT,type="so.2sp.2")

results<-createAicTable(list(m1,m2,m3,m4,m5,m6,m7,m8,m9))
summary(results)
#m9 is the best model for detection which includes mean substrate 
#and species specific detection

summary(m1)
summary(m2)
summary(m3)
summary(m4)
summary(m5)
summary(m6)
summary(m7)
summary(m8)
summary(m9)

mean(modAvg(results, param = "pA")$est)
mean(modAvg(results, param = "pB")$est)
sd(modAvg(results, param = "pA")$est)
sd(modAvg(results, param = "pB")$est)

#Occupancy
#Intercept model
m10<-occMod(model=list(psi~1,p~SP*Mean.Sub),data=WT,type="so.2sp.2")
#Species specific occupancy
m11<-occMod(model=list(psi~SP,p~SP*Mean.Sub),data=WT,type="so.2sp.2")
#Species specific occupancy + turbidity
m12<-occMod(model=list(psi~SP+Turbi,p~SP*Mean.Sub),data=WT,type="so.2sp.2")
#Species specific occupancy + depth
m13<-occMod(model=list(psi~SP+Depth,p~SP*Mean.Sub),data=WT,type="so.2sp.2")
#Species specific occupancy + sfg
m14<-occMod(model=list(psi~SP+sfg,p~SP*Mean.Sub),data=WT,type="so.2sp.2")
#Species specific occupancy*turbidity
m15<-occMod(model=list(psi~SP*Turbi,p~SP*Mean.Sub),data=WT,type="so.2sp.2")
#Species specific occupancy*depth
m16<-occMod(model=list(psi~SP*Depth,p~SP*Mean.Sub),data=WT,type="so.2sp.2")
#Species specific occupancy*sfg
m17<-occMod(model=list(psi~SP*sfg,p~SP*Mean.Sub),data=WT,type="so.2sp.2")
#add interaction term
m18<-occMod(model=list(psi~SP+INT,p~SP*Mean.Sub),data=WT,type="so.2sp.2")

m18$warnings #no issues with convergence or variance-covariance matrix
coef(m18, "psi") #psi coefficients look ok
coef(m18,"p") #p coefficients lok ok
exp(2.631133) #suggesting aggregation
unique(fitted(m18, "psiA"))
unique(fitted(m18, "psiBa"))
unique(fitted(m18, "psiBA"))

unique(fitted(m18, "nu")) #this is why we have a warning - high nu, 
#otherwise model is fine
unique(fitted(m18, "pA"))
unique(fitted(m18, "pB"))

m19<-occMod(model=list(psi~SP+INT+Turbi,p~SP*Mean.Sub),data=WT,type="so.2sp.2")
m19$warnings
coef(m19, "psi")
exp(2.724885)
unique(fitted(m19, "psiA"))
unique(fitted(m19, "psiBa"))
unique(fitted(m19, "psiBA"))
unique(fitted(m19, "nu"))
unique(fitted(m19, "pA"))
unique(fitted(m19, "pB"))

m20<-occMod(model=list(psi~INT+Turbi*SP,p~SP*Mean.Sub),data=WT,type="so.2sp.2") 
m20$warnings
m21<-occMod(model=list(psi~SP+INT+Depth,p~SP*Mean.Sub),data=WT,type="so.2sp.2")
m21$warnings
m22<-occMod(model=list(psi~INT+Depth*SP,p~SP*Mean.Sub),data=WT,type="so.2sp.2")
m22$warnings
m23<-occMod(model=list(psi~sfg*SP+INT,p~SP*Mean.Sub),data=WT,type="so.2sp.2") 
m23$warnings
m24<-occMod(model=list(psi~sfg+SP+INT,p~SP*Mean.Sub),data=WT,type="so.2sp.2") 
m24$warnings#convergence warning > 3, likely not an issue but try random starting values
set.seed(29534)
m24<-occMod(model=list(psi~sfg+SP+INT,p~SP*Mean.Sub),data=WT,type="so.2sp.2", randinit = 500)
m24$warnings

coef(m24, "psi") #clearly the interaction term is causing the issue
coef(m24,"p")

unique(fitted(m24, "psiA"))
unique(fitted(m24, "psiBa"))
unique(fitted(m24, "psiBA"))
unique(fitted(m24, "nu"))
unique(fitted(m24, "pA"))
unique(fitted(m24, "pB"))

results<-createAicTable(list(m10,m11,m12,m13,m14,
                             m15,m16,m17,m18,m19,
                             m20,m21,m22,m23,m24))
summary(results)

summary(m10)
summary(m11)
summary(m12)
summary(m13)
summary(m14)
summary(m15)
summary(m16)
summary(m17)
summary(m18)
summary(m19)
summary(m20)
summary(m21)
summary(m22)
summary(m23)
summary(m24)

coef(m23,"psi")
coef(m23,"p")
coef(m23,"psi.VC")
coef(m23,"p.VC")

unique(fitted(m23, "psiA"))
unique(fitted(m23, "psiBa"))
unique(fitted(m23, "psiBA"))
unique(fitted(m23, "nu"))
unique(fitted(m23, "pA"))
unique(fitted(m23, "pB"))

mean(modAvg(results, param = "psiA")$est)
mean(modAvg(results, param = "psiBa")$est)
mean(modAvg(results, param = "psiBA")$est)
sd(modAvg(results, param = "psiA")$est)
sd(modAvg(results, param = "psiBa")$est)
sd(modAvg(results, param = "psiBA")$est)

mean(modAvg(results, param = "psiA")$est)*mean(modAvg(results, param = "psiBA")$est)+
       (1-mean(modAvg(results, param = "psiA")$est))*mean(modAvg(results, param = "psiBa")$est)
sd(modAvg(results, param = "psiA")$est)*sd(modAvg(results, param = "psiBA")$est)+
       (1-sd(modAvg(results, param = "psiA")$est))*sd(modAvg(results, param = "psiBa")$est)

#Species Interaction Factor
t<-0
sD<-0
for(ii in 1:4){
  ma.psiA<-modAvg(results, param = "psiA",index=ii)
  ma.psiBa<-modAvg(results, param = "psiBa",index=ii)
  ma.psiBA<-modAvg(results, param = "psiBA",index=ii)
  t[ii]<-mean(modAvg(results, param = "psiA")$est)*mean(modAvg(results, param = "psiBA")$est)+
    (1-mean(modAvg(results, param = "psiA")$est))*mean(modAvg(results, param = "psiBa")$est)
  sD[ii]<-sd(modAvg(results, param = "psiA")$est)*sd(modAvg(results, param = "psiBA")$est)+
    (1-sd(modAvg(results, param = "psiA")$est))*sd(modAvg(results, param = "psiBa")$est)
}
t
sD

mean((fitted(m23, "psiA")$est))
sd(((fitted(m23,  "psiA")$est)))
mean(ma.psiA$est*ma.psiBA$est+(1-ma.psiA$est)*ma.psiBa$est)
sd(ma.psiA$est*ma.psiBA$est+(1-ma.psiA$est)*ma.psiBa$est)
############################################################################
############################################################################
#plotting functions
############################################################################
############################################################################
#Personal ggplot theme
theme_me <- theme_bw() +
  theme(axis.title = element_text(size = 11, family = "sans", face = "bold"),
        axis.text.x = element_text(size = 10, family = "sans", colour = "black"),
        axis.text.y = element_text(size = 10, family = "sans", hjust = 0.6,
                                   colour = "black"),
        legend.title = element_text(size = 10, family = "sans"),
        legend.text = element_text(size = 8, family = "sans"),
        strip.text = element_text(size = 11, family = "sans", face = "bold"))
############################################################################
############################################################################
ma.psiA<-modAvg(results, param = "psiA",index=1)
ma.psiBa<-modAvg(results, param = "psiBa",index=1)
ma.psiBA<-modAvg(results, param = "psiBA",index=1)

ma.pA<-modAvg(results, param = "pA",index=1)
ma.pB<-modAvg(results, param = "pB",index=1)
ma.prA<-modAvg(results, param = "rA",index=1)
ma.prBa<-modAvg(results, param = "rBa",index=1)
ma.prBA<-modAvg(results, param = "rBA",index=1)

mean((ma.psiA$est*ma.psiBA$est)/(ma.psiA$est*((ma.psiA$est*ma.psiBA$est)+
                                                (1-ma.psiA$est)*ma.psiBa$est)))
sd((ma.psiA$est*ma.psiBA$est)/(ma.psiA$est*((ma.psiA$est*ma.psiBA$est)+
                                              (1-ma.psiA$est)*ma.psiBa$est)))

idx<-order(Sfg.non)
w<-as.data.frame(cbind(Sfg.non[idx],ma.psiA$est[idx]))
x<-as.data.frame(cbind(Sfg.non[idx],ma.psiBa$est[idx]))
y<-as.data.frame(cbind(Sfg.non[idx],ma.psiBA$est[idx]))
z<-as.data.frame(cbind(Sfg.non[idx],ma.psiA$est[idx]*ma.psiBA$est[idx]+
                         (1-ma.psiA$est[idx])*ma.psiBa$est[idx]))

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
ma.psiA<-modAvg(results, param = "psiA",index=5)
ma.psiBa<-modAvg(results, param = "psiBa",index=5)
ma.psiBA<-modAvg(results, param = "psiBA",index=5)

idx<-order(Turbidity.non)
w<-as.data.frame(cbind(Turbidity.non[idx],ma.psiA$est[idx]))
x<-as.data.frame(cbind(Turbidity.non[idx],ma.psiBa$est[idx]))
y<-as.data.frame(cbind(Turbidity.non[idx],ma.psiBA$est[idx]))
z<-as.data.frame(cbind(Turbidity.non[idx],ma.psiA$est[idx]*ma.psiBA$est[idx]+
                         (1-ma.psiA$est[idx])*ma.psiBa$est[idx]))

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
  theme(legend.position=c(.32,.85),
        legend.text = element_text(size=12),
        legend.background=element_blank())+
  annotate("text",x=22,y=1,label="B",size=4)
turbidgg

#depth
summary(results)
ma.psiA<-modAvg(results, param = "psiA",index=9)
ma.psiBa<-modAvg(results, param = "psiBa",index=9)
ma.psiBA<-modAvg(results, param = "psiBA",index=9)

idx<-order(Depth.non)
w<-as.data.frame(cbind(Depth.non[idx],ma.psiA$est[idx]))
x<-as.data.frame(cbind(Depth.non[idx],ma.psiBa$est[idx]))
y<-as.data.frame(cbind(Depth.non[idx],ma.psiBA$est[idx]))
z<-as.data.frame(cbind(Depth.non[idx],ma.psiA$est[idx]*ma.psiBA$est[idx]+
                         (1-ma.psiA$est[idx])*ma.psiBa$est[idx]))

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
############################################################################
############################################################################
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
############################################################################
############################################################################
#Figure 3
multiplot(sfg.gg,turbidgg,depthgg,cols=3)

############################################################################
############################################################################
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
                          labels=c("Eastern sand darter", "Silver shiner"))+
  theme_me+
  theme(legend.position=c(.8,.93),
        legend.text = element_text(size=12),
        legend.background=element_blank())

detection
