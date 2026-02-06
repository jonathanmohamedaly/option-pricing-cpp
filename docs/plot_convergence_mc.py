import numpy as np
import matplotlib.pyplot as plt


data = np.loadtxt("convergence_mc.txt", delimiter=";", skiprows=1)

N = data[:, 0]
mc_price = data[:, 1]
error = data[:, 2]


C = error[0] * np.sqrt(N[0])
theoretical = C / np.sqrt(N)

plt.figure()
plt.loglog(N, error, label="Monte Carlo error")
plt.loglog(N, theoretical, "--", label=r"$\sim 1/\sqrt{N}$")

plt.xlabel("Number of simulations")
plt.ylabel("Relative error")
plt.title("Monte Carlo Convergence")
plt.legend()
plt.grid(True, which="both")
plt.show()
