import numpy as np
from TumourCA import get_neighbours
D=3
L = np.linspace(1,9,9).reshape(D,D)
print(L)
print("D: ",D)
coords=(2,2)
print("Coords: ",coords)
out = get_neighbours(coords,L,D)
print(out)