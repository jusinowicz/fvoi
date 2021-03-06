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
xec=pec     #Conditional probs of phenotype: optimal
xec = matrix(c(0.9,0.1,0.1,0.9),2,2) #To test sub-optimal strategy

#Make a joint distribution table so that both sets of marginals = P
jec = matrix(c(0.294, 0.006, 0.021, 0.679),2,2) 
cc=matrix(colSums(jec),2,2, byrow=T)  #The cue as a matrix

#Calculate this two different ways. From full equation: 
rec = matrix(c(30,0,0,15), 2,2) #Reproduction as a table
s1 = 0 

for (c in 1:2){  
	e1=0
	for (e in 1:2) { 
		p1 = 0 
		for (p in 1:2) { 
			p1 = p1 + xec[p,c]*rec[p,e]
		}
		e1 = e1+ pec[e,c]*log(p1)
	}

	s1 = s1+ colSums(jec)[c]*e1
}
rho1 = s1 

#From fitness and info components. Note Dkl = 0 because xec = pec
rec = matrix(c(30,15), 2,2) #Note the change in the definition of this term. 
f1 = sum(rowSums(jec)*log(rec[,1]) )
HCC = -sum(cc*pec*(log(xec))) #Conditional cross entropy. Breaks into: 
HEC = -sum(jec*(log(jec/cc))) #Conditional entropy, jec/cc
Dkl3=sum(jec*log(pec/xec) ) #Conditional divergence

#Growth rate. Check that these are equal, and equal to rho1
rho2a = f1-HEC-Dkl3
rho2b = f1-HCC

#Written out!
cc[1,1]*pec[1,1]*log(xec[1,1]*30) + cc[2,1]*pec[2,1]*log(xec[2,1]*15)+
cc[1,2]*pec[1,2]*log(xec[1,2]*30) + cc[2,2]*pec[2,2]*log(xec[2,2]*15)

#################################
#####Subfair and bet hedging
source("./env_functions.R")

P = c(0.3,0.7) #Environment
X = c(0.5,0.5) #Uniform
X2 = P #Proportionate
R = c(30,0.1) #Reproduction

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

#This is a second way to find the optimal germination value (in the variable 
# 1-gst$b0)using KKT. This also returns the optimal proportions of phenotypes. 
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
rh1 = a1/X
# esub1 = P*log((1-gi_opt1)/X+R*gi_opt1)       
# cesub = P*log(X)  
esub1 = P*log(rh1)       
cesub = P*log(X)  
rho1= sum(esub1+cesub)
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
rho2 = sum(esub2+cesub2)
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
#Calculate this two different ways. From full equation: 
rec = matrix(c(30,0,0,0.1), 2,2) #Reproduction as a table
ne =2 #Number of environments
s1 = 0 
for (c in 1:2){  
	e1=0
	for (e in 1:2) { 
		p1 = 0 
		for (p in 1:2) { 
			#p1 = p1 +(  (1-gi_opt1)*si + gi_opt1* xec[p,c]*rec[p,e]    )
			p1 = p1+(1-gi_opt1)*si/ne +gi_opt1*xec[p,c]*rec[p,e]
		}

		#e1 = e1+ pec[e,c]*log( (1-gi_opt1)*si + p1)
		
		#Quick note: move the additive constant into the sum via
		#(1-gi_opt1)*si / e, then this maybe becomes d_r and the 
		#y(e) are the bet hedging proportions ala pec[]?
		e1 = e1+ pec[e,c]*log(p1)

	}
	s1 = s1+ colSums(jec)[c]*e1
}

rho3 = s1 

#From fitness and info components. Note Dkl = 0 because xec = pec
rec2 = matrix(c(30,0.1), 2,2) #Note the change in the definition of this term. 
bhec1 = (1-gi_opt1)*si + gi_opt1*xec*rec2
bhec2  = (bhec1/pec)

f4 = sum(jec*log(bhec2) )
HCC4 = -sum(cc*pec*(log(pec))) #Conditional cross entropy. Breaks into: 
HEC4 = -sum(jec*(log(jec/cc))) #Conditional entropy, jec/cc
Dkl4=sum(jec*log(pec/pec) ) #Conditional divergence

#Growth rate. Check that these are equal, and equal to rho1
rho4a = f4 - HCC4
rho4b = f4 - HEC4-Dkl4


#####Showing some work: 
#No storage
sum((jec)*(log(rec)+log(xec) ) )
sum(jec*log(rec))+sum(jec*log(xec))
sum(jec*log(rec))+sum(jec*log(xec)) + sum(jec*log(pec))-sum(jec*log(pec))
sum(jec*log(rec))+sum(jec*log(pec)) - sum(jec*log(pec/xec))

#With storage
sum((jec)*(log(bhec2 )+log(pec ) ) )
sum(jec*log(bhec2))+sum(jec*log(pec))
sum(jec*log(bhec2))+sum(jec*log(pec))- sum(jec*log(pec/pec))


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

