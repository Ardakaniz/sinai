import matplotlib.pyplot as plt
import numpy as np

points = np.loadtxt("trajectories.dat", delimiter=';').swapaxes(0,1)

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

moving_atom = None

#for i in range(1, len(points[0])):
#    ps = np.array([ points[0][:i], points[1][:i],  ])
#
#    if moving_atom is not None:
#        moving_atom.remove()
#    moving_atom = plt.Circle((ps[0][-1], ps[1][-1]), 0.25, fc=(0,0,0,0),ec='g')
#    plt.gca().add_patch(moving_atom)
#
#    plt.plot(ps[0], ps[1], c='b')
#    plt.scatter(ps[0][-1], ps[1][-1],c='black')
#    plt.scatter(ps[0][:-1], ps[1][:-1], c='b')
#    ax.set_aspect('equal')
#    plt.savefig(f"figs/{i}.png")

plt.title(f"Sinai's billard : trajectories of {len(points[0]) - 1} bounces")
plt.plot(points[0], points[1], c='b')
if len(points[0]) < 100:
    plt.scatter(points[0], points[1], c='b')
ax.set_aspect('equal')
plt.show()