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

#####Subfair and bet hedging
P = c(0.7,0.3) #Environment
X = c(0.5,0.5) #Uniform
X2 = P #Proportionate
R = c(30,0.1) #Reproduction

gi = seq(0.01,1,0.01)
si = 1

Xbig = matrix(X, length(gi),2)
X2big = matrix(X2, length(gi),2,byrow=T)
Rbig = matrix(R, length(gi),2,byrow=T)
Pbig =matrix(P,length(gi), length(P),byrow=T)

pg1=Pbig*log((1-gi)*si + gi*Xbig*Rbig)
pg2=Pbig*log((1-gi)*si + gi*X2big*Rbig)

#This will show the optimal germination value: 
pg1 = rowSums(pg1) 
pg2 = rowSums(pg2) 

plot(gi,pg2 )
points(gi,pg1,col="red")

gi_opt1 = gi[which(pg1 == max(pg1)) ] 
gi_opt2 = gi[which(pg2 == max(pg2)) ] 

gst = (get_single_opt_KKT( fr=as.matrix(R), ep=as.matrix(P), nspp=1, sr = 1)) #Optimal 
X2 = c(gst$bi)

#Uniform 
esub1 = P*log((1-gi_opt1)/X+R)       
cesub = P*log(X)  
sum(esub1+cesub)
Dkl=sum(P*log(P/X) )

#Optimal 
esub2 = P*log((1-gi_opt2)/X2+R) 
cesub2 = P*log(X2)  
esub2[!is.finite(esub2)] = 0 
cesub2[!is.finite(cesub2)] = 0 
sum(esub2+cesub2) 

Dkl2=sum(P*log(P/X2) ) 
 



#####Subfair and bet hedging with arbitrary number of 
#####environmental states.
num_states = 10 #Environmental bins or states
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

gi_opt1 = gi[which(pg1 == max(pg1)) ] 
gi_opt2 = gi[which(pg2 == max(pg2)) ] 