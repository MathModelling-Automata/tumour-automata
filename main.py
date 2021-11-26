from TumourCA import run, plot_counts,plate_trace

import multiprocessing as mp
from datetime import datetime
import numpy as np

def parallel_run(function,kwargs,n_iter):
    pool= mp.Pool(mp.cpu_count())
    print("Running %s simulation(s)"%(n_iter))
    jobs=[pool.apply(function,args= kwargs) for i in range(n_iter)]
    pool.close()

    times,traces, kills,cellcounts = [[job[i] for job in jobs] for i in range(0,4)]

    traces=[job[1] for job in jobs]
    kills=np.asarray([job[2] for job in jobs])
    cellcounts=np.asarray([job[3] for job in jobs])
    return times,traces,kills,cellcounts

def get_kwargs(mot,mul,mort,agent_ratio):
    # Bundle arguments for parallel processing    
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

    return [D,n_seed,tmax,mot,mul,mort,agents,agent_params], agents

if __name__=="__main__":

    # Variables
    n_iter = 5

    mot= 0.001 # Motility
    mul= 0.023 # Growth constant
    # mort= 0.036    # Death constant
    mort=0.5    # Death constant
    agent_ratio= 2  # Number of immune cells at t=0

    param_dict= {

        "motility": mot,
        "growth": mul,
        "mortality": mort,
        "agent_ratio": agent_ratio
    }

    # Fixed parameters

    kwargs,agents = get_kwargs(mot,mul,mort,agent_ratio)

    # Run processes across parallel threads
    
    times,traces,kills,cellcounts = parallel_run(run,kwargs,n_iter)
    t = times[0]

    now = datetime.now()
    time = now.strftime("%Y%m%d_%H%M")
    title=" ".join([str(k)+": "+str(v) for (k,v) in param_dict.items()])
    plot_counts(title,time,t,cellcounts,kills,agents,save=True)
    plate_trace(title,time,traces[0],save=True)