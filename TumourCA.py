import numpy as np
import matplotlib.pyplot as plt
import multiprocessing as mp

# Seed 
def seed(D,n_seed,agents=None,agent_params=None):
    L = np.zeros((D,D))
    occ = []
    for n in range(n_seed):
        x,y = np.random.randint(0,D,size=2)
        if (x,y) in occ:
            x,y = np.random.randint(0,D,size=2)
        else:
            occ.append((x,y))
    for coords in occ:
        L[coords[0],coords[1]]=1
    if agents:
        agent_pop=[]
        for a in range(agent_params[0]):
            x,y = np.random.randint(0,D,size=2)
            if (x,y) in agent_pop:
                x,y = np.random.randint(0,D,size=2)
            else:
                agent_pop.append((x,y))
        for agent in agent_pop:
            L[agent[0],agent[1]]=-1

    return L

# Find neighbours for a given cell
def get_nonzero(L):
    pop = np.nonzero(L)
    return [(pop[0][i],pop[1][i]) for i in range(len(pop[0]))]

def get_neighbours(coords,L,D):
    x,y=coords
    
    # Find possible neighbours accounting for edge cases 

    x1 = np.max([0,x-1])
    x2 = np.min([D,x+2])
    y1 = np.max([0,y-1])
    y2= np.min([D,y+2])

    # Get neighbour coordinates
    pos_neighbours = [(a,b) for a in range(x1,x2,1) for b in range(y1,y2,1) if (a,b)!=coords]
    
    # Return neighbours and gaps
    return [n for n in pos_neighbours if L[n]>0],[n for n in pos_neighbours if L[n]==0]

def cell_death(mort,L,cell):
    r = np.random.uniform(0,1)
    if r<mort:
        L[cell]=0
    return L

def decide_fate(mot,mul):
    r2=np.random.uniform(0,1)
    if r2<mot/(mot+mul):
        move=True
        grow=False
    else:
        move=False
        grow=True
    return move, grow

def move_grow(cell,gaps,L,move,grow):
    out=L.copy()
    target=gaps[np.random.randint(0,len(gaps))]
    if move:
        out[target]=out[cell]
        out[cell]=0
        
    elif grow:
        out[target]=out[cell]

    return out

def kill(L,neighbours,kills,kill_ratio):
    out=L.copy()
    if len(neighbours)>0:
        np.random.shuffle(neighbours)
    # target = np.random.randint(0,len(neighbours),len(neighbours))
        for target in neighbours[:kill_ratio*len(neighbours)]:
            out[target]=0
            kills+=1

    return out, kills


def evolve(D,n_seed,times,mot,mul,mort,agents=None, agent_params=None):
    L = seed(D,n_seed,agents,agent_params)
    if agents:
        _,agent_mot,agent_mul,agent_mort,kill_ratio = agent_params

    t=0
    kills=0
    out=[]
    kills=[]
    k=0
    for t in times:
        out.append(L)
        coords = get_nonzero(L)
        for cell in coords:
            neighbours,gaps= get_neighbours(cell,L,D)

            # Find cancer cells
            if L[cell]==1:
                if len(gaps)<=4:
                    L = cell_death(mort,L,cell)
                if len(gaps)>0:
                    move, grow = decide_fate(mot,mul)
                    L=move_grow(cell,gaps,L,move,grow)
                    
            # Find immune cells
            elif L[cell]==-1:
                L,k=kill(L,neighbours,k,kill_ratio)
                # if len(gaps)<=1:
                #     L = cell_death(agent_mort,L,cell)
                if len(gaps)>0:
                    move, grow = decide_fate(agent_mot,agent_mul)
                    L=move_grow(cell,gaps,L,move,grow)
        kills.append(k)

    return out,kills

def cell_count(traces):
    out=[]
    for trace in traces:
        cells = np.count_nonzero(trace>0)
        out.append(cells)
    return out

def run(D,n_seed,tmax,mot,mul,mort,agents,agent_params):
    times=np.linspace(0,tmax,tmax)
    trace,kills = evolve(D,n_seed,times,mot,mul,mort,agents=agents,agent_params=agent_params)
    cellcount = cell_count(trace)
    return times,trace,kills,cellcount

def plot_counts(t,cellcounts,kills,agents=None):

    mean_counts = np.mean(cellcounts,axis=0)

    for c_count in cellcounts:
        plt.plot(t,c_count)
    plt.plot(t,mean_counts,color='black',label="mean_count")
    if agents:
        mean_kills = np.mean(kills,axis=0)
        plt.plot(t,mean_kills,color='black',linestyle='dashed',label="mean_kill_total")
    plt.legend()
    plt.xlabel("Time")
    plt.ylabel("N_cells")
    plt.show()

def plate_trace(trace):
    for image in trace:
        plt.imshow(image)
        plt.show()

if __name__=="__main__":

    # Variables

    mot= 0.001 # Motility
    mul= 0.023 # Growth constant
    # mort= 0.036    # Death constant
    mort=0.5
    agent_ratio= 2  # Number of immune cells at t=0

    # Fixed parameters
    
    D=50 # Dimensions of array
    tmax=5 # Timeout
    n_iter=3    # Repeats
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

    kwargs= [D,n_seed,tmax,mot,mul,mort,agents,agent_params]

    # Run processes across parallel threads

    
    pool= mp.Pool(mp.cpu_count())
    jobs=[pool.apply(run,args= kwargs) for i in range(n_iter)]
    pool.close()

    times,traces, kills,cellcounts = [[job[i] for job in jobs] for i in range(0,4)]
    t = times[0]
    traces=[job[1] for job in jobs]
    kills=np.asarray([job[2] for job in jobs])
    cellcounts=np.asarray([job[3] for job in jobs])

    plot_counts(t,cellcounts,kills,agents)
    plate_trace(traces[0])