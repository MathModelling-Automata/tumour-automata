from typing import no_type_check


import numpy as np
from TumourCA import evolve
from TumourCA import cell_count

def evaluate_single(dictionary):

    print("Evaluating single instance")

    # Variable parameters from dictionary
    mot,mul,mort,agent_ratio =[dictionary[key] for key in ['motility','growth','mortality','agent_ratio']]

    # Fixed parameters
    D=50 # Dimensions of array
    tmax=40 # Timeout
    n_seed = 10 # Number of cancer cells at t=0

    if agent_ratio==0:
        agents=False
    else:
        agents=True
        n_agents=int(agent_ratio*n_seed)
    agent_mot = 1   # Immune cell motility
    agent_mul = 0   # Immune cell growth constant
    agent_mort = 0  # Immune cell death constant
    kill_ratio= 1   # Efficiency of killing

    if agents:
        agent_params=[n_agents,agent_mot,agent_mul,agent_mort,kill_ratio]
    else:
        agent_params=None

    times=np.linspace(0,tmax,tmax)
    trace,kills = evolve(D,n_seed,times,mot,mul,mort,agents=agents,agent_params=agent_params)
    cellcount = cell_count(trace)
    
    return times,trace,kills,cellcount

if __name__=="__main__":

    params_var={
    'motility': 0.001, # Motility
    'growth': 0.023, # Growth constant
    'mortality': 0.036,    # Death constant
    'agent_ratio': 2 }  # Number of immune cells at t=0

    times,trace,kills,cellcount = evaluate_single(params_var)
    print(*[np.shape(x) for x in [times,trace,kills,cellcount]])