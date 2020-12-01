
# figure 2 
library("png")

null3.av <- aggregate(.~n, null3, mean) 
 
theme2 <- theme_classic()+
theme(legend.position=c(0.75,0.6), legend.title=element_blank(), legend.text=element_text(size=6),legend.background=element_blank(), 
axis.text=element_text(size=9), axis.title=element_text(size=9))

convex<-readPNG("data/img/Slide3.png")
convex<-rasterGrob(convex, interpolate=TRUE)

A <- ggplot()+
geom_crossbar(data=null1, aes(x=n, ymin=y1, y=y2, ymax=y3), colour=NA, fill="grey", alpha=0.6)+
geom_line(data=null3.av,aes(n, n.ent), col="darkred", linetype="dotted")+
geom_line(data=dat,aes(x=s, y=n.ent), col="darkred")+
geom_point(data=dat,aes(x=s, y=n.ent, shape=domain),col="darkred",fill="white",size=1)+
geom_point(data=dat, aes(s, vol, fill=region, shape=domain), size=2.5)+
labs(x="Species richness", y="Functional diversity (%)")+
scale_fill_manual(values=cols)+
scale_shape_manual(values=c(24,21))+
annotation_custom(convex,xmin=375, xmax=600, ymin=0, ymax=50)+
guides(fill="none", shape=guide_legend(override.aes = list(size=c(2))))+
theme2

neighb<-readPNG("data/img/Slide5.png")
neighb<-rasterGrob(neighb, interpolate=TRUE)

B <- ggplot()+
geom_crossbar(data=null2, aes(x=n, ymin=y1, y=y2, ymax=y3), colour=NA, fill="grey", alpha=0.6)+
geom_line(data=null3.av,aes(n, ones), col="darkred", linetype="dotted")+
geom_line(data=dat,aes(x=s, y=one), col="darkred")+
geom_point(data=dat,aes(x=s, y=one, shape=domain),col="darkred",fill="white",size=1)+
geom_point(data=dat, aes(s, nn, fill=region, shape=domain), size=2.5)+
labs(x="Species richness", y="Neighbour disimilarity (%)")+
scale_fill_manual(values=cols)+
scale_shape_manual(values=c(24,21))+
annotation_custom(neighb,xmin=375, xmax=600, ymin=35, ymax=85)+
guides(fill="none", shape="none")+
theme2

gp$domain <- info$domain[match(gp$region, info$region)]

clustplot<- ggplot(clust, aes(PC1, PC2)) +
geom_point(colour="grey", alpha=0.9, size=0.15)+
stat_ellipse(aes(group=factor(k)), geom="polygon", col="black", alpha=0.1, size=0.25)+
theme_void()+theme(panel.border=element_rect(colour="grey", size=1), plot.margin=unit(c(1,1,1,1), "mm"))

C <- ggplot(gp, aes(x=s, y=prop, group=Var1))+
geom_line(colour="grey")+
geom_point(aes(fill=region, shape=domain),stroke=0.25)+
scale_fill_manual(values=cols)+
scale_shape_manual(values=c(24,21))+ylim(c(0,0.4))+
guides(fill="none", shape="none")+
annotation_custom(ggplotGrob(clustplot),xmin=385, xmax=590, ymin=.25, ymax=.42)+
labs(x="Species richness", y="Proportion of species")+
theme2
  
 gp$rich<-ifelse(gp$s > 180, "richness > 200", "richness < 200")
D <-ggplot(gp[gp$Freq>0,], aes(x=rank, y=Freq, group=region))+
geom_line(aes(col=region))+
geom_point(aes(fill=region, shape=domain),stroke=0.25)+
facet_wrap(~rich, scales="free_y")+
scale_fill_manual(values=cols)+
scale_colour_manual(values=cols)+
scale_shape_manual(values=c(24,21))+
labs(x="Cluster", y="Number of species")+
guides(fill="none", shape="none", col="none")+
theme2+theme(strip.background=element_blank(), strip.text=element_text(size=8, face="bold"))

  
  
# plot

fig2 <- plot_grid(A,B, C, D, labels=c("A","B", "C","D"))


