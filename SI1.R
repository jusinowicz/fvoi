P = c(0.3,0.7) #Environment
X = c(0.5,0.5) #Uniform
X2 = P #Proportionate
R = c(30, 15) #Reproduction

#Uniform 
e1 = P*log(R)       
ce = P*log(X)  
sum(e1+ce)

#Proportionate 
ce2 = P*log(X2) 
sum(e1+ce2) 
Dkl=sum(P*log(P/X) ) 
Dkl2=sum(P*log(P/X2) ) 

#Information section:
pec = matrix( c(0.98,0.02,0.03,0.97),2,2) #Conditional probs of E                                                                                          
xec=pec     #Conditional probs of pop
#Make a joint distribution table so that both sets of marginals = P
jec = matrix(c(0.294, 0.006, 0.021, 0.679),2,2) 
rec = matrix(c(30,15), 2,2) #Reproduction as a table
cc=matrix(colSums(jec),2,2, byrow=T)  #The cue as a matrix

#Calculate this two different ways. From full equation: 
s1 = (pec*log(xec*rec))  
rho1 = sum(colSums(jec)*rowSums(s1) ) 

#From fitness and info components. Note Dkl = 0 because xec = pec
f1 = sum(rowSums(jec)*log(rec[,1]) )
HEC = -sum(jec*(log(jec/cc))) #Conditional entropy
HCC = -sum(cc*pec*(log(xec))) #Conditional cross entropy
rho2 = f1-HEC

#Written out!
cc[1,1]*pec[1,1]*log(xec[1,1]*30) + cc[2,1]*pec[2,1]*log(xec[2,1]*15)+
cc[1,2]*pec[1,2]*log(xec[1,2]*30) + cc[2,2]*pec[2,2]*log(xec[2,2]*15)

#################################
#####Subfair and bet hedging
source("./env_functions.R")

P = c(0.3,0.7) #Environment
X = c(0.5,0.5) #Uniform
X2 = P #Proportionate
R = c(30,1) #Reproduction

gi = seq(0.01,1,0.01)
si = 1

Xbig = matrix(X, length(gi),2)
X2big = matrix(X2, length(gi),2,byrow=T)
Rbig = matrix(R, length(gi),2,byrow=T)
Pbig =matrix(P,length(gi), length(P),byrow=T)

#This will show the optimal germination value: 
pg1=Pbig*log((1-gi)*si + gi*Xbig*Rbig)
pg2=Pbig*log((1-gi)*si + gi*X2big*Rbig)

pg1 = rowSums(pg1) 
pg2 = rowSums(pg2) 

plot(gi,pg2 )
points(gi,pg1,col="red")

gi_opt1 = gi[which(pg1 == max(pg1)) ] 
gi_opt2 = gi[which(pg2 == max(pg2)) ] 

gst = (get_single_opt_KKT( fr=as.matrix(R), ep=as.matrix(P), nspp=1, sr = 1)) #Optimal 
#The percentages stored in gst$bi are out of the total, and must be converted to be 
#out of the amount germinating to extend to the model: 
X2 = c(gst$bi/sum(gst$bi)) 

#Uniform 
# ((1-gi_opt1)+gi_opt1*X*R)^P
# r1 = P*log((1-gi_opt1)+gi_opt1*X*R)
# esub1 = P*log((1-gi_opt1)/X+gi_opt1*R)       

a1 = (1-gi_opt1)+R*gi_opt1*X
r1 = sum(P*log((1-gi_opt1)+R*gi_opt1*X) ) 
esub1 = P*log((1-gi_opt1)/X+R*gi_opt1)       

cesub = P*log(X)  
sum(esub1+cesub)
Dkl=sum(P*log(P/X) )

#Optimal
#X2= c(0.3,0.7) 
r2 = sum(P*log((1-gi_opt1)+R*gi_opt1*X2) ) 
#This matches most closely to the non-bet-hedging example,
#but may not always work
# esub2 = P*log((1-gi_opt2)/X2+R*gi_opt2) 
# cesub2 = P*log(X2) 
# esub2[!is.finite(esub2)] = 0 
# cesub2[!is.finite(cesub2)] = 0 
# sum(esub2+cesub2) 
# Dkl2=sum(P[X2>0]*log(P[X2>0]/X2[X2>0]) )  

#The effective phenotypes convert the optimal X2 to 
#bet-hedging proportions.
#These are the "effective phenotypes"
bh1 = (1-gi_opt1)+R*gi_opt1*X2 
bh2 = bh1/P
esub2 = P*log(bh2)
cesub2 = P*log(P)
sum(esub2+cesub2)
Dkl2=sum(P*log(P/P) ) 


#Laplace or additive smooth?
# al1=0.1
# (X2+al1)/(2+al1*(1:2))
 
#Information section:
pec = matrix( c(0.98,0.02,0.03,0.97),2,2) #Conditional probs of E                                                                                          
xec=pec     #Conditional probs of pop
#Make a joint distribution table so that both sets of marginals = P
jec = matrix(c(0.294, 0.006, 0.021, 0.679),2,2) 
rec = matrix(c(30,0.1), 2,2) #Reproduction as a table
cc=matrix(colSums(jec),2,2, byrow=T)  #The cue as a matrix

gst_ec = matrix(0,2,2)
#gst_ec = list()

for (i in 1:2){
	gst_ec[,i] = (get_single_opt_KKT( fr=as.matrix(rec[,i]), ep=as.matrix(xec[,i]), 
				nspp=1, sr = 1))$bi  #Optimal 
}


# gst_ec = (get_single_opt_KKT( fr=as.matrix(c(rec)), ep=as.matrix(c(jec)), 
# 				nspp=1, sr = 1))$bi  #Optimal 

#Standardize columns
xec = gst_ec/matrix(colSums(gst_ec),2,2,byrow=T)
xec[!is.finite(xec)] = 0
xec[1,2] = 1

#Calculate the growth rate 
s1 = (pec*log((1-gi_opt1)*si + gi_opt1*xec*rec))  
r3 = sum(colSums(jec)*rowSums(s1) ) 

cc * rowSums(  (pec*log( (1-gi_opt1)*si + gi_opt1*xec*rec)  )    )
pec2 = t(pec)
pce = jec/(matrix( rowSums(jec),2,2 )) 
bhec1 = (1-gi_opt1)*si + gi_opt1*xec*rec
# x= -(d*f-b*g)/(b*c-a*d)
# y= (c*f-a*g)/(b*c-a*d)
x= -(pce[2,2]*bhec1[1,1]-pce[1,2]*bhec1[2,1])/(pce[1,2]*pce[2,1]-pce[1,1]*pce[2,2])
y= (pce[2,1]*bhec1[1,1]-pce[1,1]*bhec1[2,1])/(pce[1,2]*pce[2,1]-pce[1,1]*pce[2,2])

x= -(pec2[2,2]*bhec1[1,1]-pec2[1,2]*bhec1[2,1])/(pec2[1,2]*pec2[2,1]-pec2[1,1]*pec2[2,2])
y= (pec2[2,1]*bhec1[1,1]-pec2[1,1]*bhec1[2,1])/(pec2[1,2]*pec2[2,1]-pec2[1,1]*pec2[2,2])

xy = matrix(c(x,y),2,2,byrow=T)
bhec2 = t(pce*xy)

f3 = sum(rowSums(jec)*log(bhec1[1,]) )

bhec2 = bhec1*cc/jec 
f3 = sum(rowSums(jec)*log(rowSums(bhec2)) )
HEC3 = -sum(jec*(log(jec/cc))) #Conditional entropy
rho3 = f3-HEC3

esub3 = pec*log(bhec2)
cesub3 = pec*log(pec)
sum(esub3+cesub3)

f1 = sum(rowSums(jec)*log(rec[,1]) )
HEC = -sum(jec*(log(jec/cc))) #Conditional entropy
HCC = -sum(cc*pec[xec>0]*(log(xec[xec>0]))) #Conditional cross entropy
rho2 = f1-HEC

esub3 = pec*log((1-gi_opt2)/xec+rec) 
cesub3 = pec*log(xec)  
esub3[!is.finite(esub3)] = 0 
cesub3[!is.finite(cesub3)] = 0 
sum(esub3+cesub3) 
Dkl3=sum(pec[xec>0]*log(pec[xec>0]/xec[xec>0]) ) 

pec*( log(xec*( (1-gi_opt2)/xec+rec))) 
pec*( log(( (1-gi_opt2)+rec*xec))) 

#From fitness and info components. Note Dkl = 0 because xec = pec
f1 = sum(rowSums(jec)*log(rec[,1]) )
HEC = -sum(jec*(log(jec/cc))) #Conditional entropy
HCC = -sum(cc*pec[xec>0]*(log(xec[xec>0]))) #Conditional cross entropy
rho2 = f1-HEC

#Written out!
cc[1,1]*pec[1,1]*log((1-gi_opt2)*si + gi_opt2*xec[1,1]*rec[1,1]) + cc[2,1]*pec[2,1]*log((1-gi_opt2)*si + gi_opt2*xec[2,1]*rec[2,1])+
cc[1,2]*pec[1,2]*log((1-gi_opt2)*si + gi_opt2*xec[1,2]*rec[1,2]) + cc[2,2]*pec[2,2]*log((1-gi_opt2)*si + gi_opt2*xec[2,2]*rec[2,2])



#####Subfair and bet hedging with arbitrary number of 
#####environmental states.
num_states = 10 #Environmental bins or states
ngens=1000
nspp = 1
env_states = rbinom(ngens,num_states, 0.4)
env_states = hist(env_states,0:(num_states))$counts + 1e-4
env_states = env_states/sum(env_states)

env = sample(x=(0:(num_states-1)), size=ngens, prob =env_states, replace=T)
#env = make_simple_env(env_states,ngens)
env_prob = prop.table(table(factor(env, levels = 0:(num_states-1))))

P = c(env_prob) #Environment
X = c( matrix( (1/num_states),num_states,1)) #Uniform
X2 = P #Proportionate

#Reproduction
fs_cor = 0.999 
#With method = "constant", this is just a constant value per state. 
#With method = "variable," matches based on fs_cor
#fm_method = "constant"
fm_method = "variable"
fm = matrix(num_states*8,nspp,1) # When this is a constant = num_states, fair odds
#fm = matrix(10,nspp,1) # When this is a constant = num_states, fair odds
fs = matrix(0,num_states,nspp)
for (s in 1:nspp) { fs[,s] = get_species_fit(probs=env_prob, fcor = fs_cor, fm=fm[s], 
	method=fm_method )}
hist(fs) #Show fitness distribution 

mstates=floor(num_states/2)
R = c(get_species_fit_pois(mstates, num_states, nspp,fm ))

gi = seq(0.01,1,0.01)
#gi = 0.1
si = 1
Pbig =matrix(P,length(gi), length(P),byrow=T)

pg1=Pbig*log((1-gi)*si + kronecker(gi, t(X*R)))
pg2=Pbig*log((1-gi)*si + kronecker(gi, t(c(X2)*R)))
pg1 = rowSums(pg1) 
pg2 = rowSums(pg2) 

plot(gi,pg2 )
points(gi,pg1,col="red")

gi_opt1 = gi[which(pg1 == max(pg1,na.rm=T)) ] 
gi_opt2 = gi[which(pg2 == max(pg2,na.rm=T)) ] 

gst = (get_single_opt_KKT( fr=as.matrix(R), ep=as.matrix(P), nspp=1, sr = 1)) #Optimal 

