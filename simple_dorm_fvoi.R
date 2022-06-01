#=============================================================================
# R code to measure the fitness value of information for a simple model of 
# germination and dormancy adapted from e.g. Cohen 1967, Ellner 1997. This builds
# from the basic example of Kelly betting on horses (where horses = environments
# in the ecological example)
#
# Each time-step, some proportion of each species' seeds may germinate, grow
# and reproduce. The ultimate fitness payoff of this decision is based on 
# how well the germinating seeds have predicted their environment. The optimal
# strategy is to germinate a proportion of seeds that matches on to the 
# probability of sucessful reproduction. Unlike in pure Kelly betting, the 
# population retains a proportion of its "wealth" (population) as a seed bank. 
#
# The population dynamics are: 
#
# 	Ni[t+1] = Gi[t]*Ni[t]
#
#     Gi[t] = (1-g_i)*s_i + g_i *sum(q_e* f_e)
#	
# Where g_i is the germination rate (this is b_0 in Cover and Thomas),
# f_e is the fitness of a type ni the environment (this is o_i in Cover and Thomas), and
# q_e is the proportion of germinating seeds that developed each type (b_i in C and T) 
#
# If organisms do not use information then we calculate G(E), signifying that the
# population can only respond directly to environmental states (i.e. e in E). 
# In this case there are still optimal g_i and q_e based on
# the distribution of environmental states and the payouts (f_e). 

# If organisms use information in the form of a cue C then we calculate G(E|C). In this 
# case, f_e(E|C) is a conditional probability with a distribution that essentially 
# indicates how reliable a cue C is.  
#
# The equation for FVOI is the difference between population growth with and 
# without information: 
# 
# deltaG(E;C) = G(E|C) - G(E)
# 
#=============================================================================
#=============================================================================
#Load libraries
#=============================================================================
library(tidyverse)
library(RandomFields)
library(vegan)
source("./env_functions.R")
source("./info_theory_functions.R")

#=============================================================================
#Declare variables 
#=============================================================================
ngens = 1000 #Time steps
num_states = 10 #Environmental bins or states
nspp = 2

#Survival rates: 
sr = c(1,1)

#Constant germination rate:
gsu = 0.71 #Constant germination rate for model 
gs_min = 0.01#For control over the germination rate, make this very small. 


#=============================================================================
#Stored values, i.e. population dynamics, information metrics
#=============================================================================
Ni1 = matrix(1, ngens+1,nspp) #Population w/ information
Ni2 = matrix(1, ngens+1,nspp) #Population w/ information
Ni3 = matrix(1, ngens+1,nspp) #Population w/ information
Ni4 = matrix(1, ngens+1,nspp) #Population w/ information
N_noi = matrix(1, ngens+1,nspp) #Population no information
N_noi2 = matrix(1, ngens+1,nspp) #Population no information
No = matrix(1, ngens+1,nspp) #Population optimal germination

rho_i1 = matrix(1, ngens+1,nspp) #Growth rate w/ information
rho_i2 = matrix(1, ngens+1,nspp) #Growth rate w/ information
rho_i3 = matrix(1, ngens+1,nspp) #Growth rate w/ information
rho_i4 = matrix(1, ngens+1,nspp) #Growth rate w/ information

rho_noi = matrix(1, ngens+1,nspp) #Growth rate, no information
rho_noi2 = matrix(1, ngens+1,nspp) #Growth rate, no information,from actual gr
rho_o = matrix(1, ngens+1,nspp) #Growth rate, optimal germination

env_act = matrix(1, ngens+1,1) #Realized environments from sim
env_sensed = matrix(1, ngens+1,nspp) #Realized environments from sim

sp_fit_i = matrix(1, ngens+1,nspp)
sp_fit_o = matrix(1, ngens+1,nspp)
gt_e2 = matrix(1,ngens+1,nspp)

sp_act = matrix(1, ngens+1,nspp) #Realized environments from sim


#array( matrix(1, ngens+1,nspp), dim = c(ngens+1,num_states,nspp) ) #Realized fitness from sim
qi_fit = matrix(1, ngens+1,nspp) #Realized fitness from sim
qo_fit = matrix(1, ngens+1,nspp) #Realized fitness from sim
qnoi_fit = matrix(1, ngens+1,nspp) #Realized fitness from sim
qnoi_fit2 = matrix(1, ngens+1,nspp) #Realized fitness from sim


#=============================================================================
#Make environment
# 	The environment consists of discrete bins or states, each with some 
#	probability p of occuring. These are calculated in 2 steps: 
#	1. Randomly generate a series of ps for the given number of states.
#	2. Generate sequences and determine the winner at each time step based 
#		on the largest value generated. 
#		This creates a distribution that is highly correlated to the one 
#		generated in step 1, but differs slightly. 
#	3. Count the probability of seeing a state from the simulated sequence
#=============================================================================
#1. 
#env_states = make_env_states(num_states)

#1b.
#env_states[env_states<0.02] = 0 ; env_states = env_states/sum(env_states)

#2. Binomial
env_states = rbinom(ngens,num_states, 0.4)
env_states = hist(env_states,0:(num_states))$counts

#3.Poisson
# env_states = rpois(ngens,num_states)
# num_states = max(env_states)
# env_states = hist(env_states,0:(num_states))$counts

#4.Uniform
# env_states = runif(num_states)

env_states = env_states/sum(env_states)

env = sample(x=(0:(num_states-1)), size=ngens, prob = env_states, replace=T)
#env = make_simple_env(env_states,ngens)
env_prob = prop.table(table(factor(env, levels = 0:(num_states-1))))

#=============================================================================
#Make species' responses
# 	Species respond by deciding how much to bet on each state/what proportion
#	to germinate. 
#   Each state is also assigned a value as the payoff/per-capita growth rate
#	(fitness). 
# 	We explore 3 scenarios here:
#	1. No information -- betting proportions drawn from uniform distribution.
#	2. Optimal germination -- solve for the optimal singl-species proportion
#	3. With information -- species respond to a cue that helps them predict
#	   the environmentc (g_i is conditional, g_i(E|C)
#		
#=============================================================================
##################################
###There are two ways to run this, one of which matches the betting example and 
#the other matches the dormancy model example.
#1. With gf_method = variable and fm_method = either variable or  constant: this 
#	matches betting example
#2. With gf_method = constant and fm_method = variable: this matches the dormancy
#	model.

####Fitness
#With fs_cor = 1, the fitness matches the optimal germination scheme. 
#Rare events could be weighted to be more valuable instaed by making
#this negative.
fs_cor = 0.999 
fm_method = "variable"
#fm_method = "constant"


#The conditions for fair/subfair odds are different with this model. Ellner and 
#others have shown that the optimal germination fraction is only <1 when 
#      sr*colMeans(1/fs) > 1
fm = matrix(3,nspp,1) # When this is a constant = num_states, fair odds

fs = matrix(0,num_states,nspp)
# for (s in 1:nspp) { fs[,s] = get_species_fit(probs=env_prob, fcor = fs_cor, fm=fm[s], 
# 	method=fm_method )}

mstates=floor(num_states/2)
fs = get_species_fit_pois(mstates, num_states, nspp,fm )*10

####Germination fraction


####Phenotype fraction

#With gs_cor = 1, the germination fraction is optimal, i.e. fractions 
#match the probability of a state occurring. Decreasing this 
#simulates increasingly poor or mismatched information. cor = 0 
#is no information. 
qs_cor = 0.9999
qs_cor_no = 0
qc = matrix(0.5,nspp,1) #For a constant germination fraction -- matches dormancy model
qs_noi = matrix(0,num_states,nspp) #No info
qs_o = matrix(0,num_states,nspp) #Optimal betting (i.e. proportionate)
qs_method = "variable"

for (s in 1:nspp) { 
	qs_noi[,s] = get_species_fraction(probs = env_prob, gcor = qs_cor, gc = qc[s], 
	method="constant" )
}


#####For the optimum single-species constant rate: 
tsize = 1e4
fr_opt = matrix(1, tsize,nspp)
#fr_opt = array( matrix(1, tsize,nspp), dim = c(tsize,num_states,nspp) ) 
env_opt =  matrix(1, tsize,1)
for (t in 1:tsize){
	fit_tmp = get_fit_one(env_states, fs)
	env_opt[t] = fit_tmp$env_act
	fr_opt[t,] = apply(fit_tmp$sp_fit,2,max)
	# fr_opt[t,,] = fit_tmp$sp_fit
}

gst = get_single_opt_KKT( fr=fs, ep=env_prob, nspp=nspp, sr = sr ) #Optimal 
gst$b0[gst$b0<=0] = gsu
gs_o =  matrix( c( matrix( gst$b0,1,2) ),ngens,nspp,byrow=T)
#qs_o =  matrix( c(get_single_opt( fr=fr_opt, nspp=nspp, sr = sr )$opts),num_states,nspp,byrow=T) #Optimal 
#qs_io = matrix(c(get_multi_opt(fr=fr_opt, gs=gs nspp=nspp, sr = sr ) ),num_states,nspp,byrow=T )

####Conditional germination fraction i.e. germination with information
#This function creates a table of conditional probabilities based on the
#G(E|C). Needs acc, which is a number between 0 and 1 describing how accurate
#the cue is on average (1 is perfect accuracy, 0 is none)
qec = get_cp(env_states, acc=c(0.9,0.9) )

#Make G(C|E) from G(E|C) to simulate population dynamics later:
qj = qec 
qce = qec
for(s in 1:nspp){
	#Joint probability distribution G(E,C) = G(C,E) = G(E|C) * G(C)
	qj[,,s] = qec[,,s]*matrix(env_states, num_states,num_states,byrow=T ) 
	
	####For 0 MI, scramble to Cue: 
	#qj[,,s] = matrix(runif(num_states^2), num_states,num_states,byrow=T ) 
	
	#qi(C|E) = qi(C,E)/G(E)
	qce[,,s] = qj[,,s]/matrix(rowSums(qj[,,s]),num_states,num_states)
	qce[,,s][!is.finite(qce[,,s])] = 0
}

###########
#This is key somehow:
gec = qec
incr=0.01
H1 = seq(0.01,.99,incr)
a1=sr[1]*(1-matrix(H1,length(H1),num_states))+kronecker(H1,t(fs[,1]))
#a1=matrix(sr*(1-H1), length(sr), dim(fs)[1] )+matrix(H1*fs[,1], length(sr),dim(fs)[1],byrow=T)
a2 = array( matrix(1, num_states,num_states), dim = c(num_states,num_states,nspp) )
a2[,,1] = sr[1]*(1-gec[,,1])+  gec[,,1]*matrix(fs[,1], num_states, num_states )
a2[,,2] = sr[2]*(1-gec[,,2])+  gec[,,2]*matrix(fs[,2], num_states, num_states )

g_in_e = H1[apply(a1,2,which.max)]
# g_in_e = matrix(g_in_e,env)
g_in_e2 = a2 
g_in_e2[,,1] =gec[,,1] * matrix( apply(a2[,,1],2,max) >=1, num_states, num_states,byrow=T ) 
g_in_e2[,,2] =gec[,,2] * matrix( apply(a2[,,2],2,max) >=1, num_states, num_states,byrow=T ) 

gst_ec = array( matrix(0,num_states,num_states), dim=c(num_states,num_states,nspp))
xec = gst_ec
for (s in 1:nspp){
	for (i in 1:num_states){
		gst_ec[,i,s] = (get_single_opt_KKT( fr=as.matrix(fs[,i]), ep=as.matrix(qec[,i,s]), 
					nspp=1, sr = 1))$bi  #Optimal 
	}
	#Standardize columns
	xec[,,s] = gst_ec[,,s]/matrix(colSums(gst_ec[,,s]),2,2,byrow=T)
	
}
xec[!is.finite(xec)] = 0


#g_in_e = env_states*g_in_e/(sum(env_states*g_in_e)) 
###########

#Use this for the runif samples: 
qs = cbind(env_states,env_states)	
qs_noi_tmp = seq(0,0.99,1/(length(qs)-1) )

#=============================================================================
#Population dynamics
#=============================================================================	

for (t in 1:ngens){
	fit_tmp = get_fit_one(env_states, fs)
	env_act[t] = fit_tmp$env_act #Store the env state
	sp_fit_o[t,] = fs[(env_act[t]+1), ]

	#No information, v1: totally random
	qs_noi = runif(nspp)
	qs_noi = c(1/num_states, 1/num_states)
	qnoi_fit[t,] = qs_noi 
	rho_noi[t, ] = ( sr*(1-gs_o[t,])   + gs_o[t,]*sp_fit_o[t,]  * qnoi_fit[t,] ) 
	N_noi[t+1,] = N_noi[t, ] * rho_noi[t, ] 
	N_noi[t+1,][N_noi[t+1,]<0] = 0

	#No information, v2: random sampling from "observed" germination rates
	qs_noi2 = apply(qs,2, sample)    
	qnoi_fit2[t,] = qs_noi2[1,] 
	rho_noi2[t, ] = ( sr*(1-gs_o[t,])   + gs_o[t,] * sp_fit_o[t,]  * qnoi_fit2[t,] ) 
	N_noi2[t+1,] = N_noi2[t, ] * rho_noi2[t, ] 
	N_noi2[t+1,][N_noi2[t+1,]<0] = 0
	
	#Optimal proportionate betting
	qs_o = env_prob [(env_act[t]+1)]
	rho_o[t, ] = ( sr*(1- gs_o[t,] )   + gs_o[t,]*sp_fit_o[t,] * qs_o ) 
	No[t+1,] = No[t, ] * rho_o[t, ] 
	No[t+1,][No[t+1,]<0] = 0

	#With information, conditional germination rates
	#The environment that species sensed (i.e. the cue they got), based on G(C|E)
	sp_fit_tmp = matrix((0:(num_states-1)),num_states,nspp)
	for( s in 1:nspp){ 
		env_sensed[t,s] = sample(x=(0:(num_states-1)), size=1, prob =qce[ (env_act[t]+1),,s], replace=T)
		ec = env_sensed[t,s]
		sp_fit_tmp[,s][sp_fit_tmp[,s]!=ec] = -1 #Identify losers
		sp_fit_tmp[,s][sp_fit_tmp[,s]==ec] = fs[,s][sp_fit_tmp[,s]==ec] #Set winning state to its payout
	}
	sp_fit_i[t,] = sp_fit_tmp[sp_fit_tmp>=0]

	#This version is the most similar to the kelly betting example, but does not 
	#make a lot of ecological sense. In particular, it never makes sense to bet
	#everything on a really crappy year under conditions of sub-fair odds. 
	rho_i1[t, ] = ( sr*(1-gs_o[t,] )   + 
				gs_o[t,]*sp_fit_i[t,] * qec[(env_act[t]+1) , , ][ (env_sensed[t,]+1) ]  )

	# #This version uses the germination rate that matches the sensed environment
	rho_i2[t, ] = ( sr*(1-gs_o[t,] )   + 
	  			gs_o[t,]*sp_fit_i[t,] * qs[(env_sensed[t,]+1)] )

	# #This version uses the optimal rates of reproduction from get_single_opt_KKT
	# rho_i3[t, ] = ( sr*(1- gs_o[t,] )   + 
	#   			gs_o[t,]*sp_fit_i[t,] * q_in_e2[ matrix( c( (env_sensed[t,]+1), c(env_act[t,]+1,env_act[t,]+1),1:nspp ), 2, 3 ) ] )
	# gt_e2[t, ] =  g_in_e2[ matrix( c( (env_sensed[t,]+1), c(env_act[t,]+1,env_act[t,]+1),1:nspp ), 2, 3 ) ]

	rho_i3[t, ] = ( sr*(1- gs_o[t,] )   + 
	  			gs_o[t,]*sp_fit_i[t,] * xec[ matrix( c( (env_sensed[t,]+1), c(env_act[t,]+1,env_act[t,]+1),1:nspp ), 2, 3 ) ] )
	gt_e2[t, ] =  xec[ matrix( c( (env_sensed[t,]+1), c(env_act[t,]+1,env_act[t,]+1),1:nspp ), 2, 3 ) ]


	#This version uses the ideal (max) germination rate that matches the sensed environment
	# rho_i4[t, ] = ( sr*(1-gs_o[t,] )   + 
	#   			gs_o[t,]*sp_fit_i[t,] * q_in_e[(env_sensed[t,]+1)] )

	Ni1[t+1,] = Ni1[t, ] * rho_i1[t, ] 
	Ni1[t+1,][Ni1[t+1,]<0] = 0

	Ni2[t+1,] = Ni2[t, ] * rho_i2[t, ] 
	Ni2[t+1,][Ni2[t+1,]<0] = 0

	Ni3[t+1,] = Ni3[t, ] * rho_i3[t, ] 
	Ni3[t+1,][Ni3[t+1,]<0] = 0

	Ni4[t+1,] = Ni4[t, ] * rho_i4[t, ] 
	Ni4[t+1,][Ni4[t+1,]<0] = 0
}


#Plot the population growth
nn=1:ngens
plot(log(No[,1]),t="l", ylim = c(-50,300))
lines(log(N_noi[,1]), col="red")
lines(log(N_noi2[,1]), col="red",lty = 2)
lines(log(Ni1[,1]),col="blue")
lines(log(Ni2[,1]),col="orange")
lines(log(Ni3[,1]),col="green")
lines(log(Ni4[,1]),col="yellow")

#Theoretical prediction based on optimal germination/betting strategy (gs)
lines(log(2^(nn*sum(env_prob*log2(qs[,1]*fs[,1])))),col="red",lty=4)
#Theoretical prediction when optimal germination matches actual probs
Wbp = log2(env_prob*fs[,1])
Wbp[!is.finite(Wbp)] = 0
lines(log(2^(nn*sum(env_prob*Wbp))),col="blue" )

#The theoretical population growth rate:
#log-Rate: 
lGr = colSums(matrix(env_prob,num_states,nspp)*log(qs*fs))
#Rate:
Gr = apply( (qs*fs)^matrix(env_prob,num_states,nspp),2,prod)


####The mutual information between the cue and the environment: 
env_freq = prop.table(table(env_act)) #Environment frequency
sE = shannon_D(env_freq) #Shannon entropy

#For species 1:
#Joint probability p(e,c): 
c_and_e = prop.table(table( data.frame ( e =env_act[1:ngens,1], c = env_sensed[1:ngens,1]) ))  #Joint prob between env and cue
sCgivenE = shannon_CE (c_and_e) #Conditional entropy H(C|E)

#Marginals: 
mar_e = rowSums(c_and_e) 
mar_c = colSums(c_and_e) 

#Mutual information: 
mI = sE - sCgivenE 


#####Divergences: 
#Make some data sets:
breaks = dim(c_and_e)[1]
rho_i1[ngens+1,] = rho_i1[ngens,]
rho_noi[ngens+1,] = rho_noi[ngens,]
rho_i = log(rho_i1)
r_noi = log(rho_noi)

#Make these data sets to mimic actual sampling: 
ds_noi = data.frame(qs = qnoi_fit[,1], envr = sp_fit_o[,1], rho = rho_noi[,1]  )
ds_i = data.frame(qs = qs[ (env_sensed[,1]+1) ], envr = sp_fit_i[,1], rho = rho_i [,1])

#Probability distribution of growth rates with info
b_use_i = seq(min(c(r_noi,rho_i),na.rm=T),max(c(r_noi,rho_i),na.rm=T), length.out=(breaks+1) )
rho_dist_i = hist(rho_i,breaks=b_use_i,plot = FALSE)
rho_i_p = rho_dist_i$counts/sum(rho_dist_i$counts)
rho_i_b = rho_dist_i$mids

#Average log growth rate:
rho_i_d = sum(rho_i_p*rho_i_b )

#Probability distribution of growth rates without info
#b_use_noi = seq(min(rho_noi,na.rm=T),max(rho_noi,na.rm=T), length.out=breaks)
rho_dist_noi = hist(r_noi,breaks=b_use_i,plot = FALSE)
rho_noi_p = rho_dist_noi$counts/sum(rho_dist_noi$counts)
rho_noi_b = rho_dist_noi$mids

#Average log growth rate:
rho_noi_d = sum(rho_noi_p*rho_noi_b )

b_use_i = log(sr[1]*(1-mar_noir)   + fs[1:9,1]  * mar_noir )
names(b_use_i) = as.character(0:8) 
b_use_i = sort(b_use_i) 

a3= cbind(c(0, rho_noi_p), b_use_i)
a3 = a3[sort(rownames(a3)),]  
sum(a3[,1]*a3[,2]+1E-10) 

#####Get a series of KL divergences between distibutions 
#####(KL from philentrop)y:
#First, build probability distributions from data: 
#Conditional of environment given cue p(e|c) 
c_ce = c_and_e/matrix( mar_c, length(mar_e), length(mar_c),byrow=T )

#Conditional of rho in an environment given cue rho(e|c): 
rc_df = data.frame ( r = rho_i[1:ngens,1], 
	c = env_sensed[1:ngens,1])

r_and_c = prop.table(table( rc_df ))  #Joint prob between env and cue 

#We want the rows of r_and_c to match the order of the environments/cues. 
#This is essential for comparing the distributions via KLD. This can be
#achieved by sorting rc_df along environment values so that the ordering
#of rho values now matches the ordering of environment values: 
order_key = unique( rc_df[order(rc_df[,2] ), ] )
#Use match() to achieve the reordering
r_and_c = r_and_c[match(as.character(order_key[,1]), rownames(r_and_c)), ]
#Marginals 
mar_r = rowSums(r_and_c) 
mar_cr = colSums(r_and_c) 
#Conditional
r_ce = r_and_c/matrix(mar_cr, length(mar_r), length(mar_cr),byrow=T)

#Joint distribution of rho_noi and environment
# rnoi_and_e = prop.table(table(  data.frame ( r = r_noi[1:ngens,1], 
# 	e = env_act[1:ngens,1]) ))
# mar_noir = rowSums(rnoi_and_e) 
# mar_noie = colSums(rnoi_and_e) 

# c_rnoi = rnoi_and_e/matrix(mar_noie, length(mar_noir), length(mar_noie),byrow=T)

# #Joint distribution of rho_noi and environment -- use marginals of environment
# bur = seq(0,1,1/(breaks) )
# rni = hist(gnoi_fit[,1],breaks = bur)
# gnuse = gnoi_fit[,1]

# for (b in 1:breaks){ 
# 	gnuse[gnuse >= bur[b] & gnuse <= bur[b+1] ] = bur[b]

# }

# rnoi_and_e = prop.table(table( data.frame (g = gnuse[1:ngens], c = env_sensed[1:ngens,1] ) ) ) 
rnoi_and_e = r_and_c
rnoi_and_e = matrix( mar_cr/num_states, dim(r_ce)[1],dim(r_ce)[2],byrow=T)

mar_noir = rowSums(rnoi_and_e) 
mar_noie = colSums(rnoi_and_e) 

c_rnoi = rnoi_and_e/matrix(mar_noie, length(mar_noir), length(mar_noie),byrow=T)

#KL Divergence of environment, sensed environment
kl_ec = philentropy::KL( rbind(mar_e,mar_c), unit="log" )

#Divergence between environment and rho with information 
kl1 = philentropy::KL( rbind(mar_e, mar_r), unit="log" )

#Divergence between environment and rho without information
klr = philentropy::KL( rbind(mar_e, mar_noir), unit="log" )

#Conditional KL Divergence between environment and rho
kl_egc = 0 
for( s in 1:length(mar_c) ){  
	#Inner sum is across all entries (rows) in a specific column of the conditional 
	#probability table
	kl_egc_tmp = philentropy::KL( rbind(c_ce[,s], r_ce[,s]), unit="log" )
	kl_egc = kl_egc+kl_egc_tmp
}

#This should be the equivalent of rho_noi 
sum(env_prob[1:9]*log( sr[1]*(1-gs_o[1])   + gs_o[1]*fs[1:9,1]  * mar_noir ) )  

f1 = ( sr[1]*(1-gs_o[1])   + gs_o[1]*fs[,1]  * (qs[,1]+1)) /(qs[,1]+1)
f2 = (sr[1]*(1-gs_o[1])   + gs_o[1]*fs[1:9,1]  * mar_noir)/mar_noir   
gp2 = sum(env_prob[1:9]*log(f2)) #This is actually the max achievable gr
gp1 = sum( env_prob*log(f1 ))
gp1 = sum(env_prob*log( sr*(1-gs_o[1])   + gs_o[1]*fs[,1]  *(qs[,1]+1) ) ) #Maximum achievable gr 
#gp1 = (env_prob*log( sr[1]*(1-g_in_e2)   + fs[,1]  * g_in_e2 ) ) #Maximum achievable gr 
rnoi_I = gp1-sE-klr 

#This should be the equivalent of rho_i 
ri_I = gp1 - sCgivenE - kl_egc

####Save stuff for figures
save(file ="dm_simp.var",Ni, No, N_noi, rho_noi, rho_o, rho_i, gs_o, gj, gce, gec, 
		 sE, sCgivenE, mI, mI_sim,env_act,env_sensed)
