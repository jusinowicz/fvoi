#=============================================================================
# This is code that is meant to complement the SI for 
# Usinowicz and O'Connor. It contains functions to run simulations of the 
# simple multiplicative model of populating  growth with and without storage 
# (i.e. seed dormancy), and with and without informaiton (i.e. conditional 
# reproduction strategies). 
#=============================================================================

#=============================================================================
#	Run iterations of a simple multiplicative model of population growth.
#	
#	P			The matrix of environmental probabilities
#	X			The matrix of reproduction probabilities
# 	R 			The reproductive fitness in each environment
#	ngens		Number of generations in pop process
#	nreps		Number of times to rep the sim
#=============================================================================
sim_mult = function(P  =  c(0.3,0.7), X = c(0.5,0.5), R = c(30, 15), 
		ngens= 1000, nreps = 1000) { 

	num_states = length(P)

	#Variables for sim
	Ni = matrix(1, ngens+1,1) #Population
	lami = matrix(1, ngens+1,1) #per-capita growth rate
	rhoi = matrix(1, nreps,1) #each run's log-mean growth rate


	for (r in 1:nreps){
		#Make the environment
		env = sample(x=(0:(num_states-1)), size=ngens, prob = P, replace=T)

		for (t in 1:ngens){

			#Current environment
			ec = env[t]

			#Identify species' reproductive output: 
			sp_fit = matrix((0:(num_states-1)),num_states,1)
			sp_fit[sp_fit!=ec] = -1 #Identify losers
			sp_fit[sp_fit==ec] = R[sp_fit==ec] #Set winning state to its payout
			sp_fit[sp_fit<0] = 0

			lami[t] = colSums(X*sp_fit)
			Ni[t+1] = (colSums(matrix(Ni[t],num_states, 1,byrow=T)*X*sp_fit))
		

		}

		rhoi [r] = mean(log(lami))

	}

	return(rhoi)

}

#=============================================================================
#	Run iterations of a simple multiplicative model of population growth where 
#	populations use information (i.e. access to a predictive cue).
#	
#	P			The matrix of environmental probabilities conditioned on cue
# 	J 			The matrix of joint environmental/cue probabilities
#	X			The matrix of reproduction probabilities conditioned on cue
# 	R 			The reproductive fitness in each environment
#	ngens		Number of generations in pop process
#	nreps		Number of times to rep the sim
#=============================================================================
sim_mult_info = function(P  =  matrix( c(0.98,0.02,0.03,0.97),2,2), 
	X = matrix(c(0.9,0.1,0.1,0.9),2,2), 
	R = matrix(c(30,0,0,15), 2,2), 
	J = matrix(c(0.294, 0.006, 0.021, 0.679),2,2),
		ngens= 1000, nreps = 1000) { 

	num_states = dim(P)[2]

	#Variables for sim
	Ni = matrix(1, ngens+1,1) #Population
	lami = matrix(1, ngens+1,1) #per-capita growth rate
	rhoi = matrix(1, nreps,1) #each run's log-mean growth rate

	# #Get P(E) from joint probability distribution
	PE = rowSums (J)
	# #Get P(C) from joint probability distribution
	PC = colSums (J)

	for (r in 1:nreps){

		#Make the environment
		env_cue = sample(x=(0:(num_states-1)), size=ngens, prob = PC, replace=T)
		ecs = env_cue 
		env = ecs

		for (t in 1:ngens){

			#Make the environment
			env[t] = sample(x=(0:(num_states-1)), size=1, prob =P[, (env_cue[t] +1)], replace=T)

			#Identify species' reproductive output: 
			sp_fit = matrix((0:(num_states-1)),num_states,1)
			sp_fit[sp_fit!=env[t] ] = -1 #Identify losers
			sp_fit[sp_fit==env[t] ] = diag(R)[sp_fit==env[t] ] #Set winning state to its payout
			sp_fit[sp_fit<0] = 0

			lami[t] = colSums(X[,(env_cue[t]+1)]*sp_fit)
			Ni[t+1] = (colSums(matrix(Ni[t,],num_states, 1,byrow=T)*X[,(env_cue[t]+1)]*sp_fit))
		

		}
		
		rhoi [r] = mean(log(lami))

	}

	return(rhoi)

}


#=============================================================================
#	Run iterations of a simple multiplicative model of population growth with
#	bet-hedging (i.e. storage).
#	
#	P			The matrix of environmental probabilities
#	X			The matrix of reproduction probabilities
# 	R 			The reproductive fitness in each environment
#	si 			The survival rate of stored/dormant individuals
#	gi 			The germination rate 
#	ngens		Number of generations in pop process
#	nreps		Number of times to rep the sim
#=============================================================================
sim_stor = function(P  =  c(0.3,0.7), X = c(0.5,0.5), R = c(30, 0.1), si=1, 
		gi = 0.5, ngens= 1000, nreps = 1000) { 

	num_states = length(P)

	#Variables for sim
	Ni = matrix(1, ngens+1,1) #Population
	lami = matrix(1, ngens+1,1) #per-capita growth rate
	rhoi = matrix(1, nreps,1) #each run's log-mean growth rate


	for (r in 1:nreps){
		#Make the environment
		env = sample(x=(0:(num_states-1)), size=ngens, prob = P, replace=T)

		for (t in 1:ngens){

			#Current environment
			ec = env[t]

			#Identify species' reproductive output: 
			sp_fit = matrix((0:(num_states-1)),num_states,1)
			sp_fit[sp_fit!=ec] = -1 #Identify losers
			sp_fit[sp_fit==ec] = R[sp_fit==ec] #Set winning state to its payout
			sp_fit[sp_fit<0] = 0

			lami[t] =  si*(1-gi) + gi*colSums(X*sp_fit)
			Ni[t+1] = lami[t] * Ni[t]
		

		}

		rhoi [r] = mean(log(lami))

	}

	return(rhoi)

}


#=============================================================================
#	Run iterations of a simple multiplicative model of population growth with
#	bet-hedging (i.e. storage) when populations have access to information.
#	
#	P			The matrix of environmental probabilities
#	X			The matrix of reproduction probabilities
# 	R 			The reproductive fitness in each environment
#	si 			The survival rate of stored/dormant individuals
#	gi 			The germination rate 
#	ngens		Number of generations in pop process
#	nreps		Number of times to rep the sim
#=============================================================================
sim_stor_info = function(P  =  matrix( c(0.98,0.02,0.03,0.97),2,2), 
	X = matrix(c(0.9,0.1,0.1,0.9),2,2), 
	R = matrix(c(30,0,0,0.1), 2,2), 
	J = matrix(c(0.294, 0.006, 0.021, 0.679),2,2),
		si=1, gi = 0.5, ngens= 1000, nreps = 1000) { 

	num_states = dim(P)[2]

	#Variables for sim
	Ni = matrix(1, ngens+1,1) #Population
	lami = matrix(1, ngens+1,1) #per-capita growth rate
	rhoi = matrix(1, nreps,1) #each run's log-mean growth rate

	# #Get P(E) from joint probability distribution
	PE = rowSums (J)
	# #Get P(C) from joint probability distribution
	PC = colSums (J)

	for (r in 1:nreps){

		#Make the environment
		env_cue = sample(x=(0:(num_states-1)), size=ngens, prob = PC, replace=T)
		ecs = env_cue 
		env = ecs

		for (t in 1:ngens){

			#Make the environment
			env[t] = sample(x=(0:(num_states-1)), size=1, prob =P[, (env_cue[t] +1)], replace=T)

			#Identify species' reproductive output: 
			sp_fit = matrix((0:(num_states-1)),num_states,1)
			sp_fit[sp_fit!=env[t] ] = -1 #Identify losers
			sp_fit[sp_fit==env[t] ] = diag(R)[sp_fit==env[t] ] #Set winning state to its payout
			sp_fit[sp_fit<0] = 0


			lami[t] =  si*(1-gi) + gi*colSums(X[,(env_cue[t]+1)]*sp_fit)
			Ni[t+1] = lami[t] * Ni[t]

		}
		
		rhoi [r] = mean(log(lami))

	}

	return(rhoi)

}

