import matplotlib.pyplot as plt
import numpy as np

points = np.loadtxt("phase.dat", delimiter=';').swapaxes(0,1)

fig = plt.figure()
ax = fig.add_subplot()

xticks = np.linspace(-np.pi, np.pi, 5);
xlabel = [r"$-\pi$", r"$-\frac{\pi}{2}$", r"$0$", r"$+\frac{\pi}{2}$",   r"$+\pi$"]

plt.title(f"Sinai's billard")
plt.plot(points[0], points[1], c='b')
plt.scatter(points[0], points[1], c='b')
plt.xlabel("$\\theta = (\\vec{u}, (0, 1))$")
plt.xticks(xticks)
ax.set_xticklabels(xlabel)
plt.ylabel("$\\vec{u} \\cdot \\vec{t}$")
plt.show()