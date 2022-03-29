import numpy as np
import matplotlib.pyplot as plt

def potential(x, y):
    D = 0.5
    a = 1
    ATOMS_DIST = 1

    atoms = [[-0.5, 0.5],[0.5, 0.5],[-0.5, -0.5],[0.5,-0.5]]
             
    pot = 0

    for i in range(4):
        atom = atoms[i]
        DIST = np.sqrt((x - atom[0]) ** 2 + (y - atom[1]) ** 2)
        pot += (1 - np.exp(-a * (DIST - ATOMS_DIST))) ** 2

    return D * pot

X = np.linspace(-1, 1,100)
Y = np.linspace(-1, 1,100)

Z = np.array([potential(X, Y[0])])
for x in X[1:]:
    arr = np.array([])
    for y in Y:
        arr = np.append(arr, potential(x, y))
    Z = np.append(Z, [arr], axis=0)


atoms = np.array([
        plt.Circle((-0.5, 0.5),  0.25, fc=(0,0,0,0),ec='r',linewidth=2),
        plt.Circle((0.5, 0.5),   0.25, fc=(0,0,0,0),ec='r',linewidth=2),
        plt.Circle((-0.5, -0.5), 0.25, fc=(0,0,0,0),ec='r',linewidth=2),
        plt.Circle((0.5, -0.5),  0.25, fc=(0,0,0,0),ec='r',linewidth=2),

        plt.Circle((-0.5, 0.5),  0.5, fc=(0,0,0,0),ec=(1,0,0,0.8),linestyle=':'),
        plt.Circle((0.5, 0.5),   0.5, fc=(0,0,0,0),ec=(1,0,0,0.8),linestyle=':'),
        plt.Circle((-0.5, -0.5), 0.5, fc=(0,0,0,0),ec=(1,0,0,0.8),linestyle=':'),
        plt.Circle((0.5, -0.5),  0.5, fc=(0,0,0,0),ec=(1,0,0,0.8),linestyle=':'),
    ])

for atom in atoms:
    plt.gca().add_patch(atom)

plt.imshow(Z, cmap='inferno', extent=(-1,1,-1,1))
plt.colorbar()
plt.show()