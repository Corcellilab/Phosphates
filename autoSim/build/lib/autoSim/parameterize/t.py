import numpy as np

def make_grid(n_mols):
    x = np.linspace(0, 20, n_mols+1, endpoint=False)[1:]

    x = np.array(list(zip(x,x,x)))

    grid = np.zeros((n_mols, 3))
    for i in range(0,3):
        coord = np.random.choice(x[:,i], size=n_mols, replace=False)
        grid[:,i] = coord 

    return grid

######
if __name__ == '__main__':
    make_grid(10)

