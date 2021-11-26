import numpy as np
import matplotlib.pyplot as plt

def init_conc(D, C0):
    # Create a nutrient latice
    M = np.zeros((D,D))
    M = M + C0
    return M

def seed(D,n_seed,C0):
    L = np.zeros((D,D))
    occ = []
    for n in range(n_seed):
        x,y = np.random.randint(0,D,size=2)
        if (x,y) in occ:
            x,y = np.random.randint(0,D,size=2)
        else:
            occ.append((x,y))
    for coords in occ:
        L[coords[0],coords[1]]=4
    
    N = init_conc(D,C0)-L

    return N

def diffusion(N, D, diff_rate):
    # Diffusion in the nutrient latice
    N1 = 10*np.ones((D,D))
    for x in range(D):
        for y in range(D):
            x1 = np.max([0,x-1])
            x2 = np.min([D,x+2])
            y1 = np.max([0,y-1])
            y2= np.min([D,y+2])

            # Get neighbour coordinates
            pos_neighbours = [(a,b) for a in range(x1,x2,1) for b in range(y1,y2,1) if (a,b)!=(x,y)]
    
            # Calculate concentration differences
            conc_sum = 0
            for n in pos_neighbours:
                conc_sum += N[n]
    
            conc_diff = conc_sum - len(pos_neighbours)*N[x,y]
            N1[x,y] = N[x,y] + diff_rate*conc_diff

    return N1


def evolve(D, C0, n_seed, diff_rate, tmax):
    times = times=np.linspace(0,tmax,tmax)
    N = seed(D, n_seed, C0)

    for t in times:
        N = diffusion(N, D, diff_rate)
        return N


P = evolve(50, 10, 10, 0.1, 20)
plt.imshow(P)
plt.show()
