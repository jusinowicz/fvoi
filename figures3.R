#=============================================================================
#Combined plots for FVOI paper
#=============================================================================
library(gridExtra)
require(grid)
library(tidyverse)
library(optbin)
library(reshape2)
source("../info_theory_functions/info_theory_functions.R")
#=============================================================================
# Load this data here. Other files are loaded within figure blocks due to 
# repeated name usage. 
#=============================================================================
load("./data/fvoi_plot1.var") #Figure 1,5
# fvoi_plot3.var #Figure 3
#"ni_simple.var" #Figure 2
#"env_fit2.var"# #Figure 4
# "dm_simp.var" #Figure 2 


#=============================================================================
# Figure 1 -- This is to generate some plots to be useful in the conceptual
# figure. However, most of the grunt work for this was done in Inkscape. 
#=============================================================================
# Work in progress
ngens = length(env_fit$env)

###Environment, species fitness, germination response
ngens = dim(env_fit$fr)[1]
elt = data.frame( 1:ngens, env_fit$env, env_fit$fr, env_fit$gr,runif(ngens),runif(ngens)) 
names(elt) = c("Time", "env", "fr1","fr2","gr1","gr2","rgr1","rgr2")
scale_gr1 = 0.4
m1=max(elt$fr1)
elt$gr1= elt$gr1*m1*scale_gr1
elt$fr1= elt$fr1+m1*runif(ngens,0,0.5)

#el_long = elt %>% gather(fr, repro, fr1:fr2) %>% gather(gr, germ, gr1:gr2)
#el_long = elt %>% gather(fr, repro, env:gr2) 
# el_long = elt %>% gather(fr, repro, fr1:rgr2) 
# el_long$fr[el_long$fr =="fr1" | el_long$fr =="gr1"|el_long$fr =="rgr1"] = "sp1"
# el_long$fr[el_long$fr =="fr2" | el_long$fr =="gr2"|el_long$fr =="rgr2"] = "sp2"
el_long = elt %>% gather(fr, repro, fr1) %>% gather(gr, germ, gr1)%>% 
gather(rgr, rgerm, rgr1)

el_long$fr[el_long$fr =="fr1" ] = "sp1"
el_long$gr[el_long$gr =="gr1"] = "sp2"
el_long$rgr[el_long$rgr =="rgr1"] = "sp2"

el_long$fr[el_long$fr =="fr2" ] = "sp2"
el_long$gr[el_long$gr =="gr2"] = "sp1"
el_long$rgr[el_long$rgr =="rgr2"] = "sp2"

el2 = subset(el_long, Time<21)
c_use1 = c("#35B779FF","#440154FF","#35B779FF","#440154FF" )

####The mutual information between the cue and the environment: 
# This uses hist() to bin and create breaks first, then calculates 
# frequencies: 
brks = 10
e1 = cut(elt$fr1, hist(elt$fr1,breaks=brks)$breaks )
env_freq = prop.table(table(e1)) #Environment frequency
sE = shannon_D(env_freq) #Shannon entropy

#For species 1: 
c1 = cut(elt$gr1,hist(elt$gr1,breaks=brks)$breaks )
c_and_e = prop.table(table( data.frame ( e =e1, c = c1) ))  #Joint prob between env and cue
sCgivenE = shannon_CE (c_and_e) #Conditional entropy H(C|E)

#Mutual information: 
mI = sE - sCgivenE 

#For text plotting
xpos1 = c(matrix(18,2,1))
ypos1 = c(el2$repro[el2$Time == 20],el2$germ[el2$Time == 20])
ypos1 = ypos1 + (c(1.5,0.1))
suse1 = c("Total rain(E) ","Rain in January(C)")
labels1 = data.frame( 
  label=suse1,
  x = xpos1, y =ypos1)

p1 = ggplot() + geom_line(data=el2, aes(x=Time, y=repro,color =fr )) +
geom_line(data=el2,aes(x=Time, y=germ,color =gr )) +
geom_text( data=labels1, aes(x=x, y=y, label = label,color=label) ) +
scale_color_manual(values=c_use1) +
 #scale_colour_viridis_d()+ 
 	ylab("Rain (cm)")+ 
	theme_bw() + theme(
	text = element_text(size=14),
	panel.border = element_blank(), panel.grid.major = element_blank(),
	panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
	legend.position = "none"
	)
p1



#Table of marginal environmental probabilities
me1 = data.frame(e = melt(env_freq) )
#d1 = d1[d1$value!=0,]
p2 = ggplot( me1, aes(x = 1, y = e.e1)) + 
  geom_raster(aes(fill=e.value)) + 
  scale_fill_gradient(low="grey90", high="red") +
   ylab("Total rain(E)") +
   xlab("")+
   ggtitle("P(E)")+
  theme_bw() + theme(	text = element_text(size=14),
  									axis.text.x = element_blank(), axis.ticks.x = element_blank(),
											 legend.position = "none")
p2

#Table of conditional probabilities
cp1 = melt(c_and_e)
#d1 = d1[d1$value!=0,]
p3 = ggplot( cp1, aes(x = c, y = e)) + 
  geom_raster(aes(fill=value)) + 
  scale_fill_gradient(low="grey90", high="red") +
  ggtitle("P(E|C)")+
  theme_bw() + theme(	text = element_text(size=14), axis.text.y = element_blank(),
  										axis.title.x = element_blank(),axis.title.y = element_blank(),
  										plot.title = element_text(hjust = 0.5),
  										legend.position = "none")
p3



c_use2 = c(color="#440154FF","#35B779FF","#440154FF","#35B779FF"  )


p4 = ggplot() + geom_line(data=el2, aes(x=Time, y=3*repro,color =fr ))+
geom_line(data=el2,aes(x=Time, y=3*germ,color =gr ))+
geom_line(data=el2, aes(x=Time, y=3*m1*scale_gr1*rgerm,color =rgr ),linetype = "dashed")+
scale_color_manual(values=c_use2)+
scale_y_continuous(
	sec.axis = sec_axis (~.*1/(3*m1*scale_gr1), name = "Germination rate" )
	)+
 #scale_colour_viridis_d()+ 
 	ylab("Reproduction")+ 
	theme_bw() + theme(
	text = element_text(size=14),
	panel.border = element_blank(), panel.grid.major = element_blank(),
	panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
	legend.position = "none"
	)
p4

###Lottery model
ll_sub = subset(lott_long, time <= 21)
ll_sub = subset(ll_sub, species==1 | species ==3)

ll_sub$N[ll_sub$species=="1"] = ll_sub$N[ll_sub$species=="1"]*3
ll_sub$species[ll_sub$species=="1"] = "species 1, cue"
ll_sub$species[ll_sub$species=="3"] = "species 1, no cue"

#For text plotting
xpos2 = c(matrix(10,2,1))
ypos2 = c(ll_sub$N[ll_sub$time == 1])
ypos2 = ypos2 + (c(0.2, 0))
suse2 = unique(ll_sub$species)

c_use3 = c("#35B779FF","#35B779FF"  )

p5 = ggplot() + geom_line(data=ll_sub, aes(x=time, y=N,color =species )) +
geom_smooth(data=ll_sub, method="lm", aes(x=time, y=N,color =species), se=FALSE, linetype = 1) +
geom_text( aes(x = xpos2, y = ypos2, label = suse2, color = suse2) ) + +
scale_color_manual(values=c_use3)+
 #scale_colour_viridis_d()+ 
 	ylab("Population")+ xlab("Time")+ 
 	scale_y_log10()+
	theme_bw() + theme(
	text = element_text(size=14),
	panel.border = element_blank(), panel.grid.major = element_blank(),
	panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
	legend.position = "none"
	)
p5

g=grid.arrange( arrangeGrob(p1, ncol=1, nrow=1 ),	
								arrangeGrob(p2, p3, ncol=2, nrow=1,bottom = textGrob("Rain in January(C)",gp = gpar(fontsize = 14)),
							widths=c( unit(0.15, "npc"), unit(0.5, "npc") ) ), 			
				widths=c(unit(0.5, "npc") ), 
				heights=c(  unit(0.25, "npc"),unit(0.5, "npc") )
				)

g=grid.arrange( arrangeGrob(p1, ncol=1, nrow=1 ),	
								arrangeGrob(p4, ncol=1, nrow=1 ) ,	
								arrangeGrob(p2, p3, ncol=2, nrow=1,bottom = textGrob("Rain in January(C)",gp = gpar(fontsize = 14)),
							widths=c( unit(0.15, "npc"), unit(0.5, "npc") ) ),
								arrangeGrob(p5, ncol=1, nrow=1 ),	
				widths=c(unit(0.25, "npc"),unit(0.25, "npc") ), 
				heights=c(  unit(0.25, "npc"),unit(0.5, "npc") )
				)

ggsave(file="figure1.pdf",g)



#########################################################################################################
#Basic histograms of environment


s <- subplot(
  plot_ly(x = elt$fr1, type = "histogram", showlegend=FALSE),
  plotly_empty(),
  plot_ly(x = elt$fr1, y =elt$gr1, type = "histogram2dcontour", showlegend=FALSE),
  plot_ly(y = elt$gr1, type = "histogram", showlegend=FALSE),
  nrows = 2, heights = c(0.2, 0.8), widths = c(0.8, 0.2),
  shareX = TRUE, shareY = TRUE, titleX = FALSE, titleY = FALSE
)

layout(s)

er1 = data.frame(fr1 = rnorm(1000, 20, 5), fr2 = rnorm(1000,25,1) )

xpos2 = c(8)
ypos2 = c(0.3)

hh1 = hist(er1$fr1,breaks=30)
h1 = shannon_D(hh1$counts/sum(hh1$counts))
he1 = paste( "H(E)=", round(h1,2) )

hh2 = hist(er1$fr2,breaks = hh1$breaks)
h2 = shannon_D(hh2$counts/sum(hh2$counts))
he2 = paste( "H(E)=", round(h2,2))


p0=ggplot()+geom_histogram(data=er1,aes(x=fr1,y=(..count..)/sum(..count..)),color="black",fill=c_use[1],alpha=0.4)+
ylab("Frequency")+ xlab("")+scale_x_continuous(limits = c(0, 35)) + 
scale_y_continuous(limits = c(0, 0.5)) +
theme_bw() + geom_text( aes(x = xpos2, y = ypos2, label = he1,size=20) ) + theme(
	text = element_text(size=16),
	panel.border = element_blank(), panel.grid.major = element_blank(),
	panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
	legend.position = "none"
	)
p0


p1=ggplot()+geom_histogram(data=er1, aes(x=fr2,y=(..count..)/sum(..count..)),color="black",fill=c_use[1],alpha=0.4)+
ylab("")+ xlab("Temperature")+scale_x_continuous(limits = c(0, 35))+
scale_y_continuous(limits = c(0, 0.5)) +
theme_bw()+ geom_text( aes(x = xpos2, y = ypos2, label = he2,size=20) )  + theme(
	text = element_text(size=16),
	panel.border = element_blank(), panel.grid.major = element_blank(),
	panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
	legend.position = "none"
	)
p1

g1=grid.arrange( arrangeGrob(ncol=1, nrow = 2 ,p0,p1) ) 


#Data for general KL divergence figure 
er2a = data.frame(fr1= er1$fr2, fr2=er1$fr2-0.5 )
er2b = data.frame(fr1 = er1$fr2-5, fr2= er1$fr2 )
er2aL = gather(er2a, fr)
er2bL = gather(er2b, fr)

c_use = c("#440154FF","#35B779FF" )

xpos2 = c(16)
ypos2 = c(0.3)

#Use the breaks from last part to make this meaningful:
# hk1 = hist(er2a$fr1,breaks = hh1$breaks)
# hk2 = hist(er2a$fr2, breaks = hh1$breaks)

hk1a = c(er2b$fr1, er2b$fr2 )  
hk1b = hist(hk1a,breaks=30) 

hk1 = hist(er2a$fr1,breaks = hk1b$breaks)
hk2 = hist(er2a$fr2, breaks = hk1b$breaks)

p1 = hk1$counts/sum(hk1$counts) + 1/30
p1 = p1/sum(p1)
p2 = hk2$counts/sum(hk2$counts) + 1/30
p2 = p2/sum(p2)
kdif1 = p1*log(p1/p2)
kdif1[!(is.finite(kdif1))] = 0 
n1= data.frame(d1 = c(0,kdif1), y = hk1$breaks)


# hk1 = hist(er2b$fr1,breaks = hh1$breaks)
# hk2 = hist(er2b$fr2, breaks = hh1$breaks )

hk1 = hist(er2b$fr1,breaks = hk1b$breaks)
hk2 = hist(er2b$fr2, breaks = hk1b$breaks )

p1 = hk1$counts/sum(hk1$counts) + 1/30
p1 = p1/sum(p1)
p2 = hk2$counts/sum(hk2$counts) + 1/30
p2 = p2/sum(p2)
kdif2 = p1*log(p1/p2)
#kdif2[!(is.finite(kdif2))] = 0 
n2 = data.frame(d1 = c(0,kdif2), y = hk1$breaks)


k1 = sum(kdif1)
kd1 = paste( "KLD =", round(k1,2))

k2 = sum(kdif2)
kd2 = paste( "KLD =", round(k2,2))

#Just make a prettier smoothed version for plot: 
spline_int1 = as.data.frame(spline(n1$y, n1$d1))
spline_int2 = as.data.frame(spline(n2$y, n2$d1))


p3 = ggplot()+
# geom_ribbon(data=n1, aes(x=y, ymin = 0, ymax = d1), stat="identity"  ) +
geom_bar(data = n1, aes(x=y, y =d1),stat="identity",alpha=0.5 )+
#geom_ribbon(data=spline_int1, aes(x=x, ymin=0, ymax=y), stat="identity",alpha=0.4 ) +
geom_histogram(data=er2aL, aes(x=value, y=(..count..)/sum(..count..), fill = fr ),alpha=0.4, color="black") +
scale_fill_manual(values=c_use)+ 
scale_y_continuous(limits = c(-0.1, 0.5)) +
ylab("Frequency, Divergence (nats) ") + xlab("")+scale_x_continuous(limits = c(15, 30)) +
theme_bw()+ geom_text( aes(x = xpos2, y = ypos2, label = kd1,size=20) )  + theme(
	text = element_text(size=16),
	panel.border = element_blank(), panel.grid.major = element_blank(),
	panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
	legend.position = "none"
	)
p3


p4 = ggplot()+
#geom_ribbon(data=spline_int2, aes(x=x, ymin=0, ymax=y), stat="identity",alpha=0.4 ) +
#geom_line(data=spline_int2, aes(x=x, y=y), stat="identity",size=1  ) +
geom_bar(data = n2, aes(x=y, y =d1),stat="identity",alpha=0.5)+
geom_histogram(data=er2bL, aes(x=value, y=(..count..)/sum(..count..), fill = fr ),alpha=0.4,, color="black") +
scale_fill_manual(values=c_use)+ 
scale_y_continuous(limits = c(-0.1, 0.5)) +
ylab("") + xlab("Temperature")+scale_x_continuous(limits = c(15, 30))+
theme_bw()+ geom_text( aes(x = xpos2, y = ypos2, label = kd2,size=20) )  + theme(
	text = element_text(size=16),
	panel.border = element_blank(), panel.grid.major = element_blank(),
	panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
	legend.position = "none"
	)
p4

g2=grid.arrange( arrangeGrob(ncol=1, nrow = 2, p3,p4) ) 



#The mutual information, as a joint pdf
#Generate correlated random variables: 
mi1 = MASS::mvrnorm(1000, mu=c(25,25), matrix( c(5,2,1,2),2,2 ) )
#Uncorrelated: 
mi2 = MASS::mvrnorm(1000, mu=c(25,25), matrix( c(5,0,0,2),2,2 ) )

mid1 = data.frame(mi1= mi1[,1],mi2 =mi1[,2])
mid1L = gather(mid1, mi)
mid2 = data.frame(mi1= mi2[,1],mi2 =mi2[,2])
mid2L = gather(mid2, mi)

xpos2 = c(8)
ypos2 = c(0.3)

library(RColorBrewer)
rf = colorRampPalette(rev(brewer.pal(11,'Spectral')))
r = rf(32)

kh1 = MASS::kde2d(x=mid1[,1],y=mid1[,2]  )
kh2 = MASS::kde2d(x=mid2[,1],y=mid2[,2] )
image(kh2,col=r)

xy1 = expand.grid(x=kh1$x, y=kh1$y)
xy2 = expand.grid(x=kh1$x, y=kh1$y)

kh1u = data.frame(xy1, z = c(kh1$z))
kh2u = data.frame(xy2, z = c(kh2$z))


#These three are all for the first plot: 
p5 = ggplot()+
geom_histogram(data=mid1, aes(x=mi1, y=(..count..)/sum(..count..), fill = c_use[1] ), color="black", alpha=0.4) +
scale_fill_manual(values=c_use[1])+ 
scale_y_continuous(limits = c(-0.1, 0.5)) +
ylab("") + xlab("")+scale_x_continuous(limits = c(15, 35))+
theme_bw()  + theme(
	text = element_text(size=16),
	panel.border = element_blank(), axis.line =element_blank(), axis.text.x = element_blank(),
											axis.ticks = element_blank(), 
											axis.text.y = element_blank(),
	legend.position = "none"
	)
p5

p6 = ggplot()+
geom_histogram(data=mid1, aes(y=mi2, x=(..count..)/sum(..count..), fill = c_use[2] ), color="black",alpha=0.4) +
scale_fill_manual(values=c_use[2])+ 
scale_x_continuous(limits = c(-0.1, 0.5)) +
ylab("") + xlab("")+scale_y_continuous(limits = c(15, 35))+
theme_bw()  + theme(
	text = element_text(size=16),panel.border = element_blank(),
											axis.line =element_blank(), axis.text.x = element_blank(), 
											axis.text.y = element_blank(),
											axis.ticks = element_blank(), 
	legend.position = "none"
	)
p6


p7 = ggplot()+geom_contour_filled(data = kh1u, aes(x=x,y=y,z = z ))+
scale_x_continuous(limits = c(15, 35))+
scale_y_continuous(limits = c(15, 35))+
ylab("Temperature 2") + xlab("Temperature 1")+
theme_bw()  + theme(
	text = element_text(size=16),
	legend.position = "none"
	)
p7


g3a =grid.arrange( arrangeGrob(p5, ncol=2, nrow=1, 
											widths=c( unit(0.5, "npc"), unit(0.15, "npc") ) ),	
									arrangeGrob(p7, p6, ncol=2, nrow=1,
											widths=c( unit(0.5, "npc"), unit(0.15, "npc") )  ), 		
									heights=c(  unit(0.25, "npc"),unit(0.5, "npc") )
				)


#These three are all for the second plot: 
p8 = ggplot()+
geom_histogram(data=mid2, aes(x=mi1, y=(..count..)/sum(..count..), fill = c_use[1] ), color="black") +
scale_fill_manual(values=c_use[1])+ 
scale_y_continuous(limits = c(-0.1, 0.5)) +
ylab("") + xlab("")+scale_x_continuous(limits = c(15, 35))+
theme_bw()  + theme(
	text = element_text(size=16),
	panel.border = element_blank(), axis.line =element_blank(), axis.text.x = element_blank(),
											axis.ticks = element_blank(), 
											axis.text.y = element_blank(),
	legend.position = "none"
	)
p8

p9 = ggplot()+
geom_histogram(data=mid2, aes(y=mi2, x=(..count..)/sum(..count..), fill = c_use[2] ), color="black") +
scale_fill_manual(values=c_use[2])+ 
scale_x_continuous(limits = c(-0.1, 0.5)) +
ylab("") + xlab("")+scale_y_continuous(limits = c(15, 35))+
theme_bw()  + theme(
	text = element_text(size=16),panel.border = element_blank(),
											axis.line =element_blank(), axis.text.x = element_blank(), 
											axis.text.y = element_blank(),
											axis.ticks = element_blank(), 
	legend.position = "none"
	)
p9


p10 = ggplot()+geom_contour_filled(data = kh2u, aes(x=x,y=y,z = z ))+
scale_x_continuous(limits = c(15, 35))+
scale_y_continuous(limits = c(15, 35))+
ylab("Temperature 2") + xlab("Temperature 1")+
theme_bw()  + theme(
	text = element_text(size=16),
	legend.position = "none"
	)
p10

g3b =grid.arrange( arrangeGrob(p8, ncol=2, nrow=1, 
											widths=c( unit(0.5, "npc"), unit(0.15, "npc") ) ),	
									arrangeGrob(p10, p9, ncol=2, nrow=1,
											widths=c( unit(0.5, "npc"), unit(0.15, "npc") )  ), 		
									heights=c(  unit(0.25, "npc"),unit(0.5, "npc") )
				)

g3 = grid.arrange(g3a,g3b, ncol = 2)

ggsave(file="figure1_HE1.pdf",g1)
ggsave(file="figure1_KLD1.pdf",g2)
ggsave(file="figure1_MIa.pdf",g3a)
ggsave(file="figure1_MIb.pdf",g3b)


#=============================================================================
# Figure 2
#=============================================================================
load("./data/ni_simple2.var") #Figure 2

ngens = dim(Ni_noi)[1]
ni = data.frame(1:ngens, Ni_noi[,1], Ni_i[,1])
names(ni) = c("Time","ni1","ni_i1")
nil = ni%>%gather(ni, pop, ni1:ni_i1) #%>% gather(ni_i, pop_i, ni_i1:ni_i2)
c_use = c("#440154FF","#440154FF","#440154FF","#440154FF" )
ni2 = nil[nil$Time < 100,]

#For text plotting
xpos2 = c(matrix(c(55, 60), 2,1))
ypos2 = c(ni2$pop[ni2$Time == 40]*200,ni2$pop[ni2$Time == 50])
ypos2= ypos2[c(1,4)]
suse2 = c("\u03C1(E~U)", "\u03C1(E|C)")
#suse2 = c("A","B")

p1 = ggplot() + geom_line(data=ni2,aes(x=Time, y=pop,color =ni,linetype = ni ))+
	#geom_line(data=ni2,aes(x=Time, y=pop_i,color =ni_i ),linetype = "dotted")+
	scale_color_manual(values=c_use)+
	 #scale_colour_viridis_d()+ 
	ylab("Population")+ xlab("")+
	geom_text( aes(x = xpos2, y = ypos2, label = suse2,color=suse2)) +
 	scale_y_log10()+
	theme_bw() + theme(
	text = element_text(size=14),
	panel.border = element_blank(), panel.grid.major = element_blank(),
	panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
	legend.position = "none"
	)
p1


mlogr = mean(log(rho_noi))
mlogr_i = mean(log(rhoi_i))
mI_sim =  mlogr_i -mlogr
mI

# infos = data.frame( rho = mlogr, rho_i = mlogr_i, MI_sim = mI_sim, MI =mI, DKL = mI_sim - mI )
# infos = infos %>% gather( type, infos, rho:DKL  )
# infos$type = factor(infos$type, levels = infos$type)

# p1a = ggplot() + geom_bar(data=infos, aes(x = type, y = infos), stat="identity"  ) +
# 	ylab("Fitness value/information")+ xlab("")+ scale_y_continuous(limits = c(0, 2.4))+
# 	theme_bw() + scale_x_discrete(breaks=infos$type,
#                   labels=c("\u03C1(E~U)","\u03C1(E|C)",
#                   			"\u0394 \u03C1", "I(E;C)" ) )+
# 		theme(
# 		text = element_text(size=14),
# 		panel.border = element_blank(), panel.grid.major = element_blank(),
# 		panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
# 		legend.position = "none"
# 		)
# p1a

##A stacked version, where the MI, D, and delta rho are stacked
#infos = data.frame( rho = mlogr, rho_i = mlogr_i, MI_sim = mI_sim, DKL = 0 )
infos = data.frame( rho = mlogr, rho_i = mlogr_i, MI_sim = mI_sim, DKL = 0 )
infos = infos %>% gather( name, infos, rho:DKL )
infos$type1 = c("r1","r1","it1","it2")
infos$type2 = c("rho1","rho2","infos","infos")
infos$tname = factor(infos$name, levels = infos$name)
infos$type1 = factor(infos$type1, levels = unique(infos$type1) )
infos$type2 = factor(infos$type2, levels = unique(infos$type2) )


p1a = ggplot()  + geom_bar(data=infos, aes(x = type2, y = infos,fill = type1),position="stack", stat="identity"  ) +
	ylab("")+ xlab("")+scale_y_continuous(limits = c(0, 2.4))+
	theme_bw() + scale_x_discrete(breaks=unique(infos$type2),
                     labels=c("\u03C1(E~U)","\u03C1(E|C)",
                  			"I(E;C)") )+ #, "\u0394 \u03C1", , expression(D[KL]) ) )+
		theme(
		text = element_text(size=14),
		panel.border = element_blank(), panel.grid.major = element_blank(),
		panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
		legend.position = "none"
		) + scale_fill_grey(start = 0, end = .9)
p1a

load("./data/dm_simp2.var")
ngens = dim(Ni)[1]
ni = data.frame(1:ngens, Ni[,1], No[,1])
names(ni) = c("Time","ni1","no1")
nil = ni%>%gather(ni, pop, ni1:no1) #%>% gather(ni_i, pop_i, ni_i1:ni_i2)
c_use = c("#440154FF","#440154FF","#440154FF","#440154FF" )
ni2 = nil[nil$Time < 100,]

#For text plotting
xpos3 = c(matrix(c(55, 60), 2,1))
ypos3 = c(ni2$pop[ni2$Time == 40],ni2$pop[ni2$Time == 25])
ypos3= ypos3[c(1,4)]
suse3 = c( "\u03C1(E|C)","\u03C1(E~U)")
#suse2 = c("A","B")

p2 = ggplot() + geom_line(data=ni2,aes(x=Time, y=pop,color =ni,linetype = ni ))+
	#geom_line(data=ni2,aes(x=Time, y=pop_i,color =ni_i ),linetype = "dotted")+
	scale_color_manual(values=c_use)+
	 #scale_colour_viridis_d()+ 
	ylab("")+ xlab("")+
	geom_text( aes(x = xpos3, y = ypos3, label = suse3,color=suse3)) +
 	scale_y_log10()+
	theme_bw() + theme(
	text = element_text(size=14),
	panel.border = element_blank(), panel.grid.major = element_blank(),
	panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
	legend.position = "none"
	)
p2


mlogr2 = mean(log(rho_o))
mlogr_i2 = mean(log(rho_i))
mI_sim2 =  mlogr_i2 -mlogr2
mI2 = mI

# infos2 = data.frame( rho = mlogr2, rho_i = mlogr_i2, MI_sim = mI_sim2, MI =mI2 )
# infos2 = infos2 %>% gather( type, infos, rho:MI  )
# infos2$type = factor(infos2$type, levels = infos2$type)

# p2a = ggplot() + geom_bar(data=infos2, aes(x = type, y = infos), stat="identity"  ) +
# 	ylab("")+ xlab("")+scale_y_continuous(limits = c(0, 2.4))+
# 	theme_bw() + scale_x_discrete(breaks=infos2$type,
#                   labels=c("\u03C1(E~U)","\u03C1(E|C)",
#                   			"\u0394 \u03C1", "I(E;C)" ) )+
# 		theme(
# 		text = element_text(size=14),
# 		panel.border = element_blank(), panel.grid.major = element_blank(),
# 		panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
# 		legend.position = "none"
# 		)
# p2a

##A stacked version, where the MI, D, and delta rho are stacked
infos2 = data.frame( rho = mlogr2, rho_i = mlogr_i2, MI_sim = mI_sim2, DKL = abs(mI2-mI_sim2)  )
infos2 = infos2 %>% gather( name, infos, rho:DKL )
infos2$type1 = c("r1","r1","it1","it2")
infos2$type2 = c("rho1","rho2","infos","infos")
infos2$tname = factor(infos2$name, levels = infos2$name)
infos2$type1 = factor(infos2$type1, levels = unique(infos2$type1) )
infos2$type2 = factor(infos2$type2, levels = unique(infos2$type2) )

p2a = ggplot() + geom_bar(data=infos2, aes(x = type2, y = infos,fill = type1),position="stack", stat="identity"  ) +
	ylab("")+ xlab("")+scale_y_continuous(limits = c(0, 2.4))+
	theme_bw() + scale_x_discrete(breaks=unique(infos2$type2),
                     labels=c("\u03C1(E~U)","\u03C1(E|C)",
                  			 "I(E;C)") )+ #,"\u0394 \u03C1", expression(D[KL]) ) )+
		theme(
		text = element_text(size=14),
		panel.border = element_blank(), panel.grid.major = element_blank(),
		panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
		legend.position = "none"
		)+ scale_fill_grey(start = 0, end = .9)
p2a

g=grid.arrange(arrangeGrob(p1,p2, ncol=2, nrow=1, bottom = textGrob("Time",gp = gpar(fontsize = 14)) ),
				arrangeGrob(p1a, p2a, ncol=2, nrow=1),
				widths=c(unit(0.5, "npc") ), 
				heights=c( unit(0.5, "npc"),unit(0.25, "npc") )
				)

ggsave(file="figure2.pdf",g)
#cairo_pdf(file="fvoi_box3.pdf")

#=============================================================================
# Figure 3
#=============================================================================
load("./data/fvoi_plot3.var") #Figure 3
###Draw species intrinsic fitness curve and plot them overtop the realized 
###env distribution: 
ngens = length(env_fit$env)
r1 = seq(0,1,length=ngens)
sp1 = exp(-0.5* ( (r1-env_fit$opt[1])/(env_fit$var[1]) )^2 )
sp2 = exp(-0.5* ( (r1-env_fit$opt[2])/(env_fit$var[1]) )^2 )
elt = data.frame( env_fit$env, r1, sp1, sp2) 
names(elt) = c("env", "r1","sp1","sp2")

p0=ggplot(data=elt)+geom_histogram(aes(x=env,y=..ncount..),color="black",fill="white")+
geom_line(aes(x=r1, y=sp1),color="#440154FF")+
geom_line(aes(x=r1, y=sp2),color="#35B779FF")+
ylab("Fitness value")+ xlab("Environment")+
theme_bw() + theme(
	text = element_text(size=10),
	panel.border = element_blank(), panel.grid.major = element_blank(),
	panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
	legend.position = "none"
	)
p0

###Environment, species fitness, germination response
ngens = dim(env_fit$fr)[1]
elt = data.frame( 1:ngens, env_fit$env, env_fit$fr, env_fit$gr,runif(ngens),runif(ngens)) 
names(elt) = c("Time", "env", "fr1","fr2","gr1","gr2","rgr1","rgr2")
#el_long = elt %>% gather(fr, repro, fr1:fr2) %>% gather(gr, germ, gr1:gr2)
#el_long = elt %>% gather(fr, repro, env:gr2) 
# el_long = elt %>% gather(fr, repro, fr1:rgr2) 
# el_long$fr[el_long$fr =="fr1" | el_long$fr =="gr1"|el_long$fr =="rgr1"] = "sp1"
# el_long$fr[el_long$fr =="fr2" | el_long$fr =="gr2"|el_long$fr =="rgr2"] = "sp2"
el_long = elt %>% gather(fr, repro, fr1:fr2) %>% gather(gr, germ, gr1:gr2)%>% 
gather(rgr, rgerm, rgr1:rgr2)
el_long$fr[el_long$fr =="fr1" ] = "sp1"
el_long$gr[el_long$gr =="gr1"] = "sp1"
el_long$rgr[el_long$rgr =="rgr1"] = "sp1"

el_long$fr[el_long$fr =="fr2" ] = "sp2"
el_long$gr[el_long$gr =="gr2"] = "sp2"
el_long$rgr[el_long$rgr =="rgr2"] = "sp2"

el2 = subset(el_long, Time<21)
c_use = c("#440154FF","#35B779FF" )

p1 = ggplot(data=el2) + geom_line(aes(x=Time, y=repro+2.5,color =fr ))+
geom_line(aes(x=Time, y=3*germ,color =gr ))+
geom_line(aes(x=Time, y=3*rgerm,color =rgr ),linetype = "dotted")+
scale_color_manual(values=c_use)+
scale_y_continuous(
	sec.axis = sec_axis (~.*1/(3), name = "Germination rate" )
	)+
 #scale_colour_viridis_d()+ 
 	ylab("Reproduction")+ 
	theme_bw() + theme(
	text = element_text(size=10),
	panel.border = element_blank(), panel.grid.major = element_blank(),
	panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
	legend.position = "none"
	)
p1

###Lottery model
ll_sub = subset(lott_long, time <= 21)
ll_sub = subset(ll_sub, species == "1" | species =="3")
ll_sub$species[ll_sub$species=="1"] = "species 1, cue"
#ll_sub$species[ll_sub$species=="2"] = "species 2, cue"
ll_sub$species[ll_sub$species=="3"] = "species 1, no cue"
#ll_sub$species[ll_sub$species=="4"] = "species 2, no cue"


#For text plotting
###Lottery model
ll_sub = subset(lott_long, time <= 50)
ll_sub = subset(ll_sub, species==1 | species ==3)

ll_sub$N[ll_sub$species=="1"] = ll_sub$N[ll_sub$species=="1"]*3
ll_sub$species[ll_sub$species=="1"] = "species 1, cue"
ll_sub$species[ll_sub$species=="3"] = "species 1, no cue"

#For text plotting
xpos2 = c(matrix(10,2,1))
ypos2 = c(ll_sub$N[ll_sub$time == 1])
ypos2 = ypos2 + (c(0.2, 0))
suse2 = unique(ll_sub$species)

c_use3 = c("#35B779FF","#35B779FF"  )

p2 = ggplot() + geom_line(data=ll_sub, aes(x=time, y=N,color =species )) +
geom_smooth(data=ll_sub, method="lm", aes(x=time, y=N,color =species), se=FALSE, linetype = 1) +
geom_text( aes(x = xpos2, y = ypos2, label = suse2, color = suse2) ) +
scale_color_manual(values=c_use3)+
 #scale_colour_viridis_d()+ 
 	ylab("Population")+ xlab("Time")+ 
 	scale_y_log10()+
	theme_bw() + theme(
	text = element_text(size=14),
	panel.border = element_blank(), panel.grid.major = element_blank(),
	panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
	legend.position = "none"
	)
p2


# g=grid.arrange(p0, p1, p2, widths=c(unit(0.5, "npc"), unit(0.5, "npc") ),
# 					 heights=unit(0.5, "npc"), ncol = 2,
# 					 bottom = textGrob("Time",gp = gpar(fontsize = 14) ) )

# g= grid.arrange(arrangeGrob(p0,p1, ncol=1, nrow=2),
#          arrangeGrob(p2, ncol=1, nrow=1), heights=c(4,1), widths=c(1,2))

g= grid.arrange(arrangeGrob(p0,p1, ncol=1, nrow=2),
         arrangeGrob(p2, ncol=1, nrow=1), widths=c(unit(0.5, "npc"), 
         	unit(0.75, "npc") ), heights=unit(0.5, "npc"))


ggsave(file="figure3.pdf", g)


#=============================================================================
# Figure 4
#=============================================================================
#load("./data/env_fit2.var") 
load("./data/env_fit3v4.var") 


ngens = dim(env_fit$mc2_all)[1]
nspp =  dim(env_fit$mc2_all)[3]
rhos = NULL

for (s in 1:nspp){ 

	#Convert these to data frames
	mc21 = as.data.frame(env_fit$mc2_all[,,s])
	mr21 = as.data.frame(env_fit$mr2_all[,,s])
	mc31 = as.data.frame(env_fit$mc3_all[,,s])
	mr31 = as.data.frame(env_fit$mr3_all[,,s])

	#Gather all of the runs into long format
	mc21_long = mc21%>%gather(run1, val1,V1:V50) #Competition with info
	mr21_long = mr21%>%gather(run2, val2,V1:V50) #Competition without info
	mc31_long = mc31%>%gather(run3, val3,V1:V50) #No competition, info
	mr31_long = mr31%>%gather(run4, val4,V1:V50) #No competition, no info

	#Add a column to denote species
	gr_use = paste("s",s,sep="")
	mc21_long=cbind(mc21_long, gr = gr_use)
	# mr21_long=cbind(mr21_long, gr = gr_use)
	# mc31_long=cbind(mc31_long, gr = gr_use)
	# mr31_long=cbind(mr31_long, gr = gr_use)

	#Do the subtraction for the fvoi
	rhos_tmp = NULL
	rhos_tmp = cbind(mc21_long,mr21_long,mc31_long,mr31_long )
	rhos_tmp$fvoi_noi = rhos_tmp$val3 - rhos_tmp$val4
	rhos_tmp$fvoi = rhos_tmp$val1 - rhos_tmp$val2

	#Calculate sensitivities, first step in niche and fitness differences
	rhos_tmp$s_noi = (rhos_tmp$val4 - rhos_tmp$val2)/rhos_tmp$val4 
	rhos_tmp$s_i = (rhos_tmp$val3 - rhos_tmp$val1)/rhos_tmp$val3

	niches = seq(1, 0, length = ngens) - 0.29
	rhos_tmp = cbind(Competition = niches, rhos_tmp)

	rhos = rbind(rhos,rhos_tmp)

}


r1 = rhos[rhos$Competition>niches[25],]
dr1= dim(r1)[1]
c_use = c("#440154FF","#35B779FF" )

#r1 = rhos


#For text plotting
xpos2 = c(matrix(c(0.25, 0.55), 2,1))
ypos2 = c( min(r1$val3[r1$Competition == niches[20]]),max(r1$val4[r1$Competition == niches[10]]) )
ypos2 = ypos2 - c(0.1, -0.05)
suse2 = c("Competition", "No competition")
#suse2 = c("A","B")

p0 = ggplot() +
	geom_point(data=r1, aes(x=Competition, y=fvoi, group = gr,color=gr ))+
	geom_smooth(data=r1, method="lm" , formula = y ~ poly(x, 2), aes(x=Competition, y=fvoi, group = gr,color=gr ) )+
	#geom_point(data=r1, aes(x=Competition, y=val4, group = gr,color=gr ),shape=5 )+
	#geom_smooth(data=r1, method="lm" , aes(x=Competition, y=val4, group = gr,color=gr ) )+
	scale_color_manual(values=c_use)+
	ylab("Fitness value of information")+ xlab("")+
	geom_hline(yintercept=0 ,linetype = "dashed")+
	#geom_text( aes(x = xpos2, y = ypos2, label = suse2)) +
	theme_bw() + theme(
	text = element_text(size=14),
	panel.border = element_blank(), panel.grid.major = element_blank(),
	panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
	legend.position = "none"
	)
p0


#For text plotting
xpos3 = c(matrix(c(0.25, 0.35), 2,1))
ypos3 = c( max(r1$val2[r1$Competition == niches[23]]) ,min(r1$val1[r1$Competition == niches[18]]))
ypos3 = ypos3 - c(0.05, -0.4)
suse3 = c( "No information","Information")
#suse2 = c("A","B")

#Niche differences:
y_i = 1-sqrt( subset(r1, r1$gr =="s1")$s_i*subset(r1, r1$gr =="s2")$s_i)
y_noi = 1-sqrt( subset(r1, r1$gr =="s1")$s_noi*subset(r1, r1$gr =="s2")$s_noi)

#Fitness differences (these should be the same): 
f_i = sqrt( subset(r1, r1$gr =="s1")$s_i/subset(r1, r1$gr =="s2")$s_i)
f_noi = sqrt( subset(r1, r1$gr =="s1")$s_noi/subset(r1, r1$gr =="s2")$s_noi)


rn1 = data.frame(Competition = r1$Competition, gr = "s1", y_i=y_i,y_noi=y_noi,f_i = f_i, f_noi=f_noi)

p1 = ggplot() +
	geom_point(data=rn1, aes(x=Competition, y=y_i, group = gr,color=gr ) )+
	geom_smooth(data=rn1, method="lm" , formula = y ~ poly(x, 2), aes(x=Competition, y=y_i, group = gr,color=gr ) )+
	geom_point(data=rn1, aes(x=Competition, y=y_noi, group = gr,color=gr ),shape=5  )+
	geom_smooth(data=rn1, method="lm" , formula = y ~ poly(x, 2),aes(x=Competition, y=y_noi, group = gr,color=gr ) )+
	scale_color_manual(values=c_use) +
	ylab("Niche difference")+ xlab("")+
	geom_hline(yintercept=0 ,linetype = "dashed")+
	geom_text( aes(x = xpos3, y = ypos3, label = suse3)) +
	theme_bw() + theme(
	text = element_text(size=14),
	panel.border = element_blank(), panel.grid.major = element_blank(),
	panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
	legend.position = "none"
	)
p1


# #For text plotting
# xpos3 = c(matrix(c(0.25, 0.35), 2,1))
# ypos3 = c( max(r1$val2[r1$Competition == niches[23]]) ,min(r1$val1[r1$Competition == niches[18]]))
# ypos3 = ypos3 - c(0.25, -0.15)
# suse3 = c( "No information","Information")
# #suse2 = c("A","B")

# p1 = ggplot() +
# 	geom_point(data=r1, aes(x=Competition, y=val1, group = gr,color=gr ) )+
# 	geom_smooth(data=r1, method="lm" , formula = y ~ poly(x, 3), aes(x=Competition, y=val1, group = gr,color=gr ) )+
# 	geom_point(data=r1, aes(x=Competition, y=val2, group = gr,color=gr ),shape=5  )+
# 	geom_smooth(data=r1, method="lm" , formula = y ~ poly(x, 3),aes(x=Competition, y=val2, group = gr,color=gr ) )+
# 	scale_color_manual(values=c_use)+
# 	ylab("Fitness")+ xlab("")+
# 	geom_hline(yintercept=0 ,linetype = "dashed")+
# 	geom_text( aes(x = xpos3, y = ypos3, label = suse3)) +
# 	theme_bw() + theme(
# 	text = element_text(size=14),
# 	panel.border = element_blank(), panel.grid.major = element_blank(),
# 	panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
# 	legend.position = "none"
# 	)
# p1

g=grid.arrange(p0, p1, widths=c(unit(0.5, "npc"), unit(0.5, "npc") ),
					 heights=unit(0.5, "npc"), ncol = 2,
					 bottom = textGrob("Resource overlap",gp = gpar(fontsize = 14) ) )

ggsave(file="fig4.pdf", g)


#=============================================================================
# Figure 5
#=============================================================================
###Social info model
both_long_use = subset(both_long, time <= 60 )
both_long_use$species[both_long_use$species=="1"] = "species 1, no info"
both_long_use$species[both_long_use$species=="2"] = "species 2, no info"
both_long_use$species[both_long_use$species=="3"] = "species 1, social info"
both_long_use$species[both_long_use$species=="4"] = "species 2, social info"

#For text plotting
xpos = c(matrix(40,4,1))
ypos = c(both_long_use$N[both_long_use$time == 60])
ypos = ypos + (c(-1,1,-1,1))
suse = unique(both_long_use$species)

p0 =ggplot()+ geom_line( data = both_long_use, aes ( x = time, y = N, color = species)  )+ 
	geom_text( aes(x = xpos, y = ypos, label = suse, color = suse) ) +
	ylab("Population")+ xlab("Time")+   scale_colour_viridis_d()+
	theme_bw() + theme(
	text = element_text(size=14),
	panel.border = element_blank(), panel.grid.major = element_blank(),
	panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
	legend.position = "none"
	)

p0
ggsave(file="fvoi_box1.pdf", p0)




#=============================================================================
#Misc old plots: 
#=============================================================================
load("env_fit1.var")
rhos = data.frame(1:ngens,  env_fit$mc2_all, env_fit$mr2_all, env_fit$mc2_all-env_fit$mr2_all,env_fit$mc3_all-env_fit$mr3_all )
names(rhos) = c("Competition", "cc1","cc2", "cr1","cr2","crho1","crho2","rho1","rho2")
rhos_long = rhos %>% gather(cc, r_cc, cc1:cc2 )%>%
					gather(cr, r_cr, cr1:cr2 )%>%
					 gather(crho, r, crho1:crho2 )%>% 
					 gather(rho, rr, rho1:rho2 )

r1 = rhos_long[rhos_long$Competition<25,]
p0 = ggplot() +
geom_line(data=r1, aes(x=Competition, y=r, group = crho ) )+
geom_line(data=r1, aes(x=Competition, y=rr, group = rho ) )
p0


p1 = ggplot() + geom_line(data=r1, aes(x=Competition, y=r_cc, group = cc  ) )+
geom_line(data=r1, aes(x=Competition, y=r_cr, group = cr  ) )
p1

#=============================================================================
#Base R plots: 
#=============================================================================

n_plot = 6000
col_use = c("black", "#440154FF")
col_use2 = c("#20A387FF","#FDE725FF")
plot(1:n_plot, cc_noi_out[1:n_plot,2], col = col_use[1], t="l",axes=F,cex.main=1.3,cex.lab=1.3, xlab = "Time", ylab="Population",ylim = c(0,40), xaxs="i",yaxs="i")
for (p in 2:(nPsp+1)){ 
	lines(1:n_plot,cc_noi_out[1:n_plot,p],col= col_use[p-1] )
	lines(1:n_plot,cc_i_out[1:n_plot,p],col= col_use2[p-1] )
}

axis(2, at=seq(0,40,10),cex.axis=1)
axis(1, at=seq(0,n_plot,2000),cex.axis=1)


#Information gain as a function of model
plot(env_fit$mc2_all[,1]-env_fit$mr2_all[,1],col="red")                                                                  
points(env_fit$mc3_all[,1]-env_fit$mr3_all[,1])   

plot(env_fit$mc2_all[,2]-env_fit$mr2_all[,2],col="red")                                                                  
points(env_fit$mc3_all[,2]-env_fit$mr3_all[,2])   

#Gain in rho across competition.
plot(env_fit$mc2_all[,1],col="red", ylim=c(0,2) )                                                                  
points(env_fit$mr2_all[,1])       

plot(env_fit$mc2_all[,2],col="red", ylim=c(0,2) )                                                                  
points(env_fit$mr2_all[,2])       


#=============================================================================
# Figure X in ms? 
#=============================================================================
load("env_fit1.var")
rho_data = cbind(env_fit$mc2_all,env_fit$mr2_all,env_fit$mc3_all,env_fit$mr3_all)
rho_data = as.data.frame(rho_data)
colnames(rho_data) = c("species1_comp_info","species2_comp_info",
						"species1_comp_noinfo","species2_comp_noinfo",
						"species1_info","species2_info",
						"species1_noinfo","species2_noinfo")
#Convert to long: 
rho_long = rho_data %>% gather( species, rho )
p1 =ggplot()+ geom_line( data = rho_long, aes ( x = rho, y = N, color = species)  )+ 
	geom_text( aes(x = xpos2, y = ypos2, label = suse2, color = suse2) ) +
	ylab("")+ xlab("")+   scale_colour_viridis_d()+ 
	theme_bw() + theme(
	text = element_text(size=14),
	panel.border = element_blank(), panel.grid.major = element_blank(),
	panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
	legend.position = "none"
	)
