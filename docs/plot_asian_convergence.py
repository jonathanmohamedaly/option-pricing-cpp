import pandas as pd
import matplotlib.pyplot as plt

mc_anti = pd.read_csv("convergence_mc_antithetic_asian.txt", sep=";")
mc_control = pd.read_csv("convergence_mc_control_variate_asian.txt", sep=";")
mc_control_geo = pd.read_csv("convergence_mc_control_variate_geometric.txt", sep=";")

plt.figure(figsize=(10,6))
plt.plot(mc_anti['nSimulations'], mc_anti['Price'], label="Monte Carlo Antithetic", alpha=0.7)
plt.plot(mc_control['nSimulations'], mc_control['Price'], label="Monte Carlo Control Variate", alpha=0.7)
plt.plot(mc_control_geo['nSimulations'], mc_control_geo['Price'], label="Monte Carlo Control Variate Geo", alpha=0.7)

plt.xscale('log')
plt.yscale('log')
plt.xlabel("Number of simulations")
plt.ylabel("Value of the option")
plt.title(" Monte Carlo Antithetic vs Control Variate vs Control Variate Geo")
plt.legend()
plt.grid(True, which="both", ls="--", alpha=0.5)

plt.show()
