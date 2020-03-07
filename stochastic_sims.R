library(deSolve)
library(adaptivetau)

source("/Users/lukestrgar/Documents/lab/params.R")
source("/Users/lukestrgar/Documents/lab/model_funcs_tools.R")
source("/Users/lukestrgar/Documents/lab/stochastic_models.R")
source("/Users/lukestrgar/Documents/lab/deterministic_models.R")

####### MDA STOCHASTIC SIMULATIONS #####################################
################################################################################################

### MDA Intervention
cov = 0.8
eff = 0.99
mda.years = c(2:21)

### Starting State
start = c(S = 5000, # susceptible snails
          E = 2000, # infected snails
          I = 500, # infected snails
          Wt = 72, # worms in treated population
          Wu = 72) # worms in untreated population

### Model Transitions
transitions = list(
  c(S = 1),             #New snail born
  c(S = -1),            #Susceptible snail dies
  c(S = -1, E = 1),     #Susceptible snail becomes exposed
  c(E = -1),            #Exposed snail dies
  c(E = -1, I = 1),     #Exposed snail becomes Infected
  c(I = -1),            #Infected snail dies
  c(Wt = 1, Wu = 1),  #Infected snail emits cercaria that produces an adult worm
  c(Wt = -1),           #Adult worm in the treated population dies
  c(Wu = -1))           #Adult worm in the untreated population dies


#function to simulate transmission over 61 years (1 year transmission spin up, 20 yrs MDA, 40 yrs recovery)
stoch.sim = function(init, cov, sim) {
  
  ## Set coverage as global variable
  assign('cov', cov, envir = .GlobalEnv)
  
  ## Run model to eq. 
  eq_traj = as.data.frame(ode(y=init,times=seq(0,200*365,30),
                              func=schisto_MDA_ODE,parms=params,method="ode23"))
  eq = eq_traj[length(seq(0,200*365,30)), c(2:6)]
  
  ## Save model traj
  w.traj.mda = c(w.traj.mda, (cov*eq$Wt)+(1-cov)*eq$Wu)
  assign('w.traj.mda', w.traj.mda, envir = .GlobalEnv)
  
  init1 = setNames(as.numeric(round(eq)), colnames(eq))

  set.seed(sim)
  
  #simulate 1 year of transmission
  fill[[1]] = ssa.adaptivetau(init1, transitions, 
                              sfx_mda, params, tf=365)
  ## Save model traj
  w.traj.mda = c(w.traj.mda, (cov*as.data.frame(fill[[1]])$Wt) + (1-cov)*as.data.frame(fill[[1]])$Wu)
  assign('w.traj.mda', w.traj.mda, envir = .GlobalEnv)

  for(m in 2:21){    #simulate 20 years of MDA
    init = setNames(as.numeric(fill[[m-1]][dim(fill[[m-1]])[1],c(2:6)]), 
                    colnames(fill[[m-1]])[c(2:6)]) #reset initial states
    
    w.pre = init[4]
    init[4] = round(init[4]* (1-eff))  #apply MDA
    w.post = init[4]
    w.remov[sim] = w.remov[sim] + (w.pre - w.post)
    print(w.remov)
    assign('w.remov', w.remov, envir = .GlobalEnv)
    

    fill[[m]] = ssa.adaptivetau(init, transitions, 
                                sfx_mda, params, tf=365) #stochastic sim for a year
    
    #Save model traj
    w.traj.mda = c(w.traj.mda, (cov*as.data.frame(fill[[m]])$Wt) + (1-cov)*as.data.frame(fill[[m]])$Wu)
    assign('w.traj.mda', w.traj.mda, envir = .GlobalEnv)
    
    fill[[m]][,1] = fill[[m]][,1] + (365*(m-1)+(m-1))    #adjust time
  }
  
  for(f in 22:years){
    init = setNames(as.numeric(fill[[f-1]][dim(fill[[f-1]])[1],c(2:6)]), 
                    colnames(fill[[f-1]])[c(2:6)]) #reset initial states
    
    #init[4] = round(init[4]* (1-eff))  #NO MDA
    
    fill[[f]] = ssa.adaptivetau(init, transitions, 
                                sfx_mda, params, tf=365) #stochastic sim for a year
    
    #Save model traj
    w.traj.mda = c(w.traj.mda, (cov*as.data.frame(fill[[f]])$Wt) + (1-cov)*as.data.frame(fill[[f]])$Wu)
    assign('w.traj.mda', w.traj.mda, envir = .GlobalEnv)
    
    fill[[f]][,1] = fill[[f]][,1] + (365*(f-1)+(f-1))    #adjust time
  }
  
  matfin = do.call(rbind,fill)
  fill.test[c(1:nrow(matfin)), , sim] = cbind(matfin, Wm = cov*matfin[,5] + (1-cov)*matfin[,6])
  assign('fill.test', fill.test, envir = .GlobalEnv)
  assign('w.remov', w.remov, envir = .GlobalEnv)
}


### Number of simulations and parameters to simulate over
cov.range = seq(.5,1,length=par.sims) #coverage values
par.sims = 1
year.days = as.numeric()
for(i in 1:20){
  year.days[i] = 365*i + (i-1)
}
stoch.sims = 5
years = 61
fill = list()
fin.vals = c()

### vector for tracking W traj. NOTE -- this is a global variable and is setup to be written to for a single simulation
  ### (i.e. par.sims = stoch.sims = 1)
w.traj.mda = c()

# vector for tracking avg W removed for each cov. value
w.remov.avg = c()

### Simulation Loop
for(i in 1:par.sims) {

    fill.test = array(data = NA, dim = c(years*365*3, 7, stoch.sims))    #array to fill with all simulations
    w.remov = array(data=0, dim=c(stoch.sims)) #vector to count total number of parasites removed for specific value of cov.
    assign('w.remov', w.remov, envir = .GlobalEnv)
    pe1 = as.numeric()        #Vector for number of chains that go extinct

    #Run sims for parameter set 
    sapply(c(1:stoch.sims), stoch.sim, init = start, cov = cov.range[i], simplify = T)
    
    #Get probability of elimination (P(e)) as number of chains that lead to extinction out of all chains
    for(p in 1:stoch.sims) {
      if(sum(fill.test[max(which(!is.na(fill.test[ , 7, p]))), , p][3:7]) == 0){ 
        pe1[p] = 1   #if no exposed, infected snails and no adult worms, elimination = 1
      } else {
        pe1[p] = 0   #else, elimination = 0
      }
    }
    w.remov.avg[i] = (sum(w.remov) / stoch.sims)
    
    fin.vals[length(fin.vals)+1] = sum(pe1) / stoch.sims
} ## End of Simulation Loop

## Plot model trajectory
mda_times = seq(1, length(w.traj.mda) - 1)
plot(mda_times, w.traj.mda, type='l', col="red")
