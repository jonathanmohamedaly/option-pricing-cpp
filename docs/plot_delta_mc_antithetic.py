import pandas as pd
import matplotlib.pyplot as plt

mc = pd.read_csv("convergence_delta_mc.txt", sep=";")
mc_anti = pd.read_csv("convergence_delta_antithetic.txt", sep=";")

plt.figure(figsize=(10,6))
plt.plot(mc['nSimulations'], mc['Relative_Error'], label="Monte Carlo standard", alpha=0.7)
plt.plot(mc_anti['nSimulations'], mc_anti['Relative_Error'], label="Monte Carlo Antithetic", alpha=0.7)

plt.xscale('log')  
plt.yscale('log')   
plt.xlabel("Number of simulations")
plt.ylabel("Relative error")
plt.title("Delta :  Monte Carlo vs Antithetic vs Control Variate")
plt.legend()
plt.grid(True, which="both", ls="--", alpha=0.5)

plt.show()
