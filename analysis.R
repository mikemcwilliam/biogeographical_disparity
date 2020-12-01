
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

list <- c("cat_growthrate","cat_skeletaldensity", "cat_corallitesize","cat_colonydiameter","cat_colonyheight", "cat_SA_vol","cat_spacesize")
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
  hulls <-rbind(hulls, cbind(h, region=r, s=nrow(df), vol))
  points <- rbind(points, cbind(df, region=r, s=nrow(df)))
  }
 
ggplot()+
 geom_polygon(data=hulls, aes(PC1, PC2), fill=NA, col="black")+
 geom_point(data=points, aes(PC1, PC2), size=.1)+
 facet_wrap(~region)+
 guides(col="none")

################################## ----- figure 1 

source("R/map.R")
source("R/figure_1.R")

png("figs/fig1.png", width=17.8, height=10.5, units="cm", res = 300)
fig1
dev.off()

################################## ----- 4D hypervolumes

vols <- unique(hulls[,c("region","s","vol")])
dat <- merge(info, vols)

null1 <- NULL #volume of random species sample
reps <- 5 

svol <- function(){
	samp <- axes[sample(nrow(axes), size=i, replace=F),]
	convhulln(samp,"FA")$vol / convhulln(axes,"FA")$vol*100
	}

for (i in seq(10,600,5)){
	df <- replicate(reps, svol(), simplify="matrix") 
	df <- quantile(df, probs = c(0.01,0.5,0.99))
	null1 <- rbind(null1, cbind(n=i,y1=df[1], y2=df[2], y3=df[3]))
	null1 <- data.frame(null1)
	}

# null1 <- read.csv("output/null1.csv") # 100 reps

ggplot()+
geom_crossbar(data=null1, aes(x=n, ymin=y1, y=y2, ymax=y3))+
geom_point(data=dat, aes(s, vol), col="red")+
labs(x="richness", y="% of global volume")

################################## ----- neighbour distances
library("FNN")

n <- 5  #neighbours

nn <- NULL
for (r in info$region){
  df <- data[!is.na(pres_abs[,r]),list]
  dist <- dist(scale(df))
  knn <- get.knn(dist,k=n)$nn.dist
  tnn <- rowSums(knn)
  nn <- rbind(nn, data.frame(d=tnn/max(tnn)*100, region=r, s=nrow(df)))
  }

nn.av <-aggregate(d~region, nn, mean)
dat$nn <- nn.av$d[match(dat$region, nn.av$region)]

null2 <- NULL #distances of random species sample
reps <- 5

snn <- function(){
	samp <- data[sample(nrow(data), size=i, replace=F),list]
	dist <- dist(scale(samp))
	knn <- get.knn(dist,k=n)$nn.dist
	rowSums(knn)
	}

for (i in seq(10,600,5)){
	df <- replicate(reps, snn(), simplify="matrix") 
	df <- colMeans(apply(df, 2, function(x){x/max(x)*100}))
	df <- quantile(df, probs = c(0.01,0.5,0.99))
	null2 <- rbind(null2, cbind(n=i,y1=df[1], y2=df[2], y3=df[3]))
	null2 <- data.frame(null2)
	}

# null2 <- read.csv("output/null2.csv") # 100 reps

ggplot()+
geom_crossbar(data=null2, aes(x=n, ymin=y1, y=y2, ymax=y3))+
geom_point(data=dat, aes(s, nn), col="red")+
labs(x="richness", y="neighbour dissimilarity")

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
	one <- length(f[f==1])/length(f)*100 # singletons
	ents <- rbind(ents, data.frame(n.ent, one, region=r))
	}
dat$n.ent <- ents$n.ent[match(dat$region, ents$region)]
dat$one <- ents$one[match(dat$region, ents$region)]

null3 <- NULL  # species spread evenly among FEs 
reps <- 100

s.fe <- function(){
	samp <- sample(c(1:length(count)), i, replace=T)
	f <- as.vector(table(samp))
	n.ent <- length(f) / length(count)*100
	ones <- length(f[f==1]) / length(f)*100
	c(n.ent, ones)
	}
	
for (i in seq(10,600,5)){
	df <- replicate(reps, s.fe(), simplify="matrix") 
	null3 <- rbind(null3, data.frame(n=i, n.ent=df[1,], ones=df[2,]))
	}

ggplot()+
geom_point(data=null3, aes(n, n.ent), shape=21)+
geom_point(data=null3, aes(n, ones), shape=21, col="grey")+
geom_point(data=dat, aes(s, n.ent), col="red")+
geom_point(data=dat, aes(s, one), col="blue")

################################## ----- functional groups
library("mclust")

mod<-mclustBIC(axes[,c(1:2)])
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
	gp <-rbind(gp, cbind(freq, prop=props, region=r, s=nrow(df)))
}

ggplot(gp, aes(x=s, y=prop, group=Var1))+geom_line()+geom_point()

gp$rich<-ifelse(gp$s > 180, "high", "low")
ggplot(gp[gp$Freq>0,], aes(x=rank, y=Freq, group=region))+geom_line()+geom_point()+facet_wrap(~rich, scales="free_y")

################################## ----- plot figure 2
source("R/figure_2.R")

png("figs/fig2.png", width=17.8, height=13.5, units="cm", res = 300)
fig2
dev.off()


################################## ----- density plots








#################################






# FIGURE S1

hulls34 <- NULL
for (r in provs){
df <- axes[pres_abs[,r]==1,]
hulls34 <-rbind(hulls34, cbind(df[chull(df$PC3, df$PC4),], province=r))
}

hulls34$province <- factor(hulls34$province, levels=provs)

ggplot()+
geom_point(data=axes, aes(PC3, PC4), col="grey", size=0.8)+
geom_polygon(data=axes[chull(axes$PC3, axes$PC4),], aes(PC3, PC4), col="grey", fill=NA)+
geom_polygon(data=hulls34, aes(PC3, PC4, fill=province), col="black", alpha = 0.6)+
geom_point(data=points, aes(PC3, PC4, fill=province), size=0.8, shape=21, stroke=0.25)+
facet_wrap(~province)+
scale_fill_manual(values=info$col)+guides(fill="none")+
theme_classic()+theme(strip.text=element_text(size=7), strip.background=element_blank(), axis.text=element_text(size=7))
