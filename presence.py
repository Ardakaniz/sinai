import matplotlib.pyplot as plt
import numpy as np

presence = np.loadtxt("presence.dat")
WIDTH = int(np.sqrt(len(presence)))
presence = presence.reshape((WIDTH,-1))

fig = plt.figure()
ax = fig.add_subplot()

plt.xlim(-0.5,0.5)
plt.ylim(-0.5,0.5)

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

plt.title("Sinai's billard: number of passage through pixels")
plt.imshow(presence, cmap='magma',extent=(-0.5,0.5, 0.5,-0.5))
plt.colorbar()
ax.get_xaxis().set_ticklabels([])
ax.get_yaxis().set_ticklabels([])
ax.set_axisbelow(True)
ax.set_aspect('equal')
plt.show()