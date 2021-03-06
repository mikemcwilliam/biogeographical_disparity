
rm(list=ls())

library("ggplot2")
library("cowplot")

data <- read.csv("data/traitbiogeography.csv", row.names=1)

info <- read.csv("data/info.csv")
pres_abs <- data[,info$name]
colnames(pres_abs) <- info$region
head(pres_abs)

cols<-info$col 
names(cols)<-info$region

################################## ----- PCA

list <- c("cat_growthrate","cat_skeletaldensity", "cat_corallitesize",
          "cat_colonydiameter","cat_colonyheight", "cat_SA_vol","cat_spacesize")
traits <- na.omit(data[,list])

pca <- prcomp(traits, center=TRUE, scale.=TRUE)
biplot(pca)

vars <- round((pca$sdev^2 / sum(pca$sdev^2)), 3)*100
sum(vars[c(1:4)])

axes <- data.frame(PC1 = pca$x[,1], PC2 = -pca$x[,2], 
PC3=pca$x[,3], PC4 =pca$x[,4])

vecs<-data.frame(varnames=rownames(pca$rotation), pca$rotation)
vecs$label<-c("3","1","7","5","6","2","4")

################################## ----- Convex hulls
library("geometry")

hulls <- NULL
points <- NULL

for (r in info$region){
  df <- axes[!is.na(pres_abs[,r]),]
  h <- df[chull(df$PC1, df$PC2),]
  vol <- convhulln(df,"FA")$vol / convhulln(axes,"FA")$vol *100 
  hulls <-rbind(hulls, cbind(h, region=r, n=nrow(df), vol))
  points <- rbind(points, cbind(df, region=r, n=nrow(df)))
  }
 
ggplot()+
  geom_polygon(data=hulls, aes(PC1, PC2), fill=NA, col="black")+
  geom_point(data=points, aes(PC1, PC2), size=.1)+
  facet_wrap(~region)+
  guides(col="none")

################################## ----- figure 1 
source("R/map.R")
map
source("R/figure_1.R")
fig1
#ggsave("figs/fig1.png", fig1, width=17.8, height=10.5, units="cm", dpi = 300)

################################## ----- 4D hypervolumes

vols <- unique(hulls[,c("region","n","vol")])
dat <- merge(info, vols)

null1 <- NULL #volume of random species sample
reps <- 5 

svol <- function(){
	samp <- axes[sample(nrow(axes), size=n, replace=F),]
	convhulln(samp,"FA")$vol / convhulln(axes,"FA")$vol*100
	}

for (n in seq(10,600,5)){
	df <- replicate(reps, svol(), simplify="matrix") 
	df <- quantile(df, probs = c(0.01,0.5,0.99))
	null1 <- rbind(null1, cbind(n,y1=df[1], y2=df[2], y3=df[3]))
	null1 <- data.frame(null1)
	}

# null1 <- read.csv("output/null1.csv") # 100 reps

ggplot()+
  geom_crossbar(data=null1, aes(x=n, ymin=y1, y=y2, ymax=y3), col="grey")+
  geom_text(data=dat, aes(n, vol, label=abv))+
  labs(y="% of global volume", x="sp richness")

################################## ----- neighbour distances
library("FNN")

neighbs <- 5  

nn <- NULL
for (r in info$region){
  df <- data[!is.na(pres_abs[,r]),list]
  dist <- dist(scale(df))
  knn <- get.knn(dist,k=neighbs)$nn.dist
  tnn <- rowSums(knn)
  nn <- rbind(nn, data.frame(d=tnn/max(tnn)*100, region=r, n=nrow(df)))
  }

nn.av <-aggregate(d~region, nn, mean)
dat$nn <- nn.av$d[match(dat$region, nn.av$region)]

null2 <- NULL #distances of random species sample
reps <- 5

snn <- function(){
	samp <- data[sample(nrow(data), size=n, replace=F),list]
	dist <- dist(scale(samp))
	knn <- get.knn(dist,k=neighbs)$nn.dist
	rowSums(knn)
	}

for (n in seq(10,600,5)){
	df <- replicate(reps, snn(), simplify="matrix") 
	df <- colMeans(apply(df, 2, function(x){x/max(x)*100}))
	df <- quantile(df, probs = c(0.01,0.5,0.99))
	null2 <- rbind(null2, cbind(n,y1=df[1], y2=df[2], y3=df[3]))
	null2 <- data.frame(null2)
	}

#null2 <- read.csv("output/null2.csv") # 100 reps

ggplot()+
  geom_crossbar(data=null2, aes(x=n, ymin=y1, y=y2, ymax=y3), col="grey")+
  geom_text(data=dat, aes(n, nn, label=abv))+
  labs(y="neighbour dissimilarity", x="sp richness")

################################## ----- functional entities
bw <- 0.58
bins <- ggplot(axes, aes(PC1, PC2))+stat_bin2d(binwidth=bw)
count <- ggplot_build(bins)$data[[1]]$count
bins + ggtitle(paste("bins =",length(count)))

ents <- NULL
for (r in info$region){
	df <- axes[!is.na(pres_abs[,r]),]
	p <- ggplot(df, aes(PC1, PC2))+stat_bin2d(binwidth=bw)
	f <- ggplot_build(p)$data[[1]]$count
	n.ent <- length(f)/length(count)*100 # number occupied
	ones <- length(f[f==1])/length(f)*100 # singletons
	ents <- rbind(ents, data.frame(n.ent, ones, region=r))
	}
dat$n.ent <- ents$n.ent[match(dat$region, ents$region)]
dat$ones <- ents$ones[match(dat$region, ents$region)]

null3 <- NULL  # species spread evenly among FEs 
reps <- 100

s.fe <- function(){
	samp <- sample(c(1:length(count)), n, replace=T)
	f <- as.vector(table(samp))
	n.ent <- length(f) / length(count)*100
	ones <- length(f[f==1]) / length(f)*100
	c(n.ent, ones)
	}
	
for (n in seq(10,600,5)){
	df <- replicate(reps, s.fe(), simplify="matrix") 
	null3 <- rbind(null3, data.frame(n, n.ent=df[1,], ones=df[2,]))
	}

ggplot()+
  geom_point(data=melt(setDT(null3), "n"), 
             aes(n, value, col=variable), shape=21, alpha=0.1)+
  geom_point(data=melt(setDT(dat), "n", c("ones","n.ent")), 
             aes(n,value, col=variable), size=2)+
  labs(y="% of entites", x="sp richness")

################################## ----- functional groups
# library("mclust") - disables maps? 
# mod<-mclustBIC(axes[,c(1:2)])
# mod   

nc <- 8
k <- kmeans(axes[,c(1:2)], centers=nc, nstart=25, iter.max=10000000)
clust<-cbind(axes, k=k$cluster)

ggplot(clust, aes(PC1, PC2)) +
  geom_point(colour="grey", alpha=0.9, size=0.15)+
  stat_ellipse(aes(group=factor(k)))

gp <- NULL
for (r in info$region){
	#r <- "Caribbean"
	df <- clust[!is.na(pres_abs[,r]),]
	freq <-data.frame(table(factor(df$k, levels=c(1:nc)))) 
	freq$rank[order(-freq$Freq)] <- 1:nrow(freq)
	props <- freq$Freq /nrow(df)
	gp <-rbind(gp, cbind(freq, prop=props, region=r, n=nrow(df)))
}

ggplot(gp, aes(x=n, y=prop, group=Var1))+
  geom_line()+geom_point()

gp$rich<-ifelse(gp$n > 180, "high", "low")
ggplot(gp[gp$Freq>0,], aes(x=rank, y=Freq, group=region))+
  geom_line()+geom_point()+
  facet_wrap(~rich, scales="free_y")

################################## ----- plot figure 2
source("R/figure_2.R")
fig2
#ggsave("figs/fig2.png", fig2, width=17.8, height=13.5, units="cm", dpi = 300)

################################## ----- KERNEL DENSITY ESTIMATION
library("ks")

doms <- c("Caribbean", "Australia")

par(mfrow=c(1,2))
for (d in doms){
	df <- axes[!is.na(pres_abs[,d]),c(1:2)]
	H <- Hpi(x=df) # optimal bandwidth estimation
	est<- kde(x=df, H=H, compute.cont=TRUE)  # kernel density estimation
	cl<-contourLevels(est, prob=c(.5,.25,.05), approx=TRUE) # contour probs
	
plot(est, cont=seq(1,100,by=1), display="filled.contour2", add=FALSE, 
     ylab="", xlab="", cex.axis=0.75, ylim=c(-4,4), xlim=c(-4, 4.3),las=1, main=d) 
plot(est,abs.cont=cl[1], labels=c(0.5),labcex=0.5, add=TRUE, lwd=0.75)
plot(est,abs.cont=cl[2], labels=c(0.75),labcex=0.5, add=TRUE, lwd=0.5)
plot(est,abs.cont=cl[3], labels=c(0.95),labcex=0.5, add=TRUE, lwd=0.5)
points( df[,], pch=16, cex=0.25, col="black") 
}

################################## ----- Importance for functions

# geometric models
source("R/geometry.R")
geo <- read.csv("data/geometry.csv")

axes$region <- ifelse(pres_abs$Caribbean %in% 1, "Caribbean", 
                      ifelse(pres_abs$Australia %in% 1, "Great Barrier Reef",  NA))

# apply geometric models based on morphology
t0 <- NULL
t1 <- NULL
for (i in 1:nrow(data)){
	gf <- data[i,"raw_growth_form"]
	funct <- geo$function.[geo$morphology==gf]
	pars <- geo[geo$morphology==gf, c("br","bh", "bpa", "th", "tb")]
	pars$d <- data[i,"dat_colonydiameter"]
	t0 <- cbind(t0, do.call(Map, c(f= as.name(funct), pars)))
	growth <- data[i,"dat_growth"]
	pars$d <- pars$d + growth/10
	t1 <- cbind(t1, do.call(Map, c(f= as.name(funct), pars)))
}

funs <- cbind(axes, do.call("rbind", t0))
funs$VOL2 <- cbind(do.call("rbind", t1))$VOL
funs$Accretion <- (funs$VOL2 - funs$VOL) * data[,"dat_skeletal"]
funs$Biomass <- funs$SA * funs$TB

plot_grid(
  ggplot(funs)+
    geom_histogram(aes(Accretion), bins=10)+scale_x_log10(),
  ggplot(funs)+
    geom_histogram(aes(Biomass), bins=10)+scale_x_log10()
  )

################################## ----- plot figure 3
source("R/figure_3.R")
fig3
#ggsave("figs/fig3.png", fig3, width=8.7, height=19, units="cm", dpi = 300)


