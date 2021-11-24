import numpy as np
import matplotlib.pyplot as plt

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
    x2 = np.min([D,x+1])
    y1 = np.max([0,y-1])
    y2= np.min([D,y+1])

    # Get neighbour coordinates
    pos_neighbours = [(a,b) for a in range(x1,x2) for b in range(y1,y2) if (a,b)!=coords]
    
    # Return the neighbour coordinates not occupied by a cell
    return [n for n in pos_neighbours if L[n]==0]

def cell_death(mort,L,cell):
    r = np.random.uniform(0,1)
    if r>mort:
        L[cell]=0
    return L

def move_grow(cell,gaps,L,move,grow):
    out=L.copy()
    target=gaps[np.random.randint(0,len(gaps))]
    if move:
        out[target]=out[cell]
        out[cell]=0
        
    elif grow:
        out[target]=out[cell]

    return out

def evolve(D,n_seed,tmax,mot,mul,mort,agents=None, agent_params=None):

    L = seed(D,n_seed,agents,agent_params)
    if agents:
        _,agent_mot,agent_mul,agent_mort = agent_params

    t=0
    out=[]
    while t<tmax:
        out.append(L)
        coords = get_nonzero(L)
        for cell in coords:
            gaps=  get_neighbours(cell,L,D)

            # Find cancer cells
            if L[cell]==1:
                if len(gaps)<=1:
                    L = cell_death(mort,L,cell)
                if len(gaps)!=0:
                    r2=np.random.uniform(0,1)
                    if r2<mot/(mot+mul):
                        move=True
                        grow=False
                    else:
                        move=False
                        grow=True
                    L=move_grow(cell,gaps,L,move,grow)
                    
            # Find immune cells
            elif L[cell]==-1:
                if len(gaps)<=1:
                    L = cell_death(agent_mort,L,cell)
                if len(gaps)!=0:
                    r2=np.random.uniform(0,1)
                    if r2<agent_mot/(agent_mot+agent_mul):
                        move=True
                        grow=False
                    else:
                        move=False
                        grow=True
                    L=move_grow(cell,gaps,L,move,grow)
        t+=1

    return out

def cell_count(traces):
    out=[]
    for trace in traces:
        cells = np.count_nonzero(trace)
        out.append(cells)
    return np.linspace(0,len(traces),len(traces)),out



D=50
n_seed = 10
mot = 0.001
mul = 0.023
mort = 0.036
n_agents = 5
agent_mot = 1
agent_mul = 0
agent_mort = 0
agents=True

if agents:
    agent_params=[n_agents,agent_mot,agent_mul,agent_mort]
else:
    agent_params=None

tmax=10
n_iter=2


combined=[]
means=[]
for n in range(n_iter):
    print("Runing simulation %s/%s"%(str(n+1),str(n_iter)))
    trace = evolve(D,n_seed,tmax,mot,mul,mort,agents=agents,agent_params=agent_params)
    combined.append(trace)
    t, cellcount = cell_count(trace)
    means.append(cellcount)
    plt.plot(t,cellcount)
plt.plot(t,np.mean(np.asarray(means),axis=0),color='black',label='mean')
plt.xlabel("Time")
plt.ylabel("Cell count")
plt.legend()
plt.show()

for trace in combined:
    plt.imshow(trace[0])
    plt.imshow(trace[-1])
    plt.show()