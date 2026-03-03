import pandas as pd
import matplotlib.pyplot as plt

#plot of MC vs Antithetic
mc = pd.read_csv("convergence_delta_asian.txt", sep=";")

plt.figure(figsize=(10,6))
plt.plot(mc['nSim'], mc['RelErr_MC'], label="Monte Carlo standard", alpha=0.7)
plt.plot(mc['nSim'], mc['RelErr_Ant'], label="Monte Carlo Antithetic", alpha=0.7)
plt.plot(mc['nSim'], mc['RelErr_CV'], label="Monte Carlo Control Variate Geo", alpha=0.7)

plt.xscale('log')
plt.yscale('log')
plt.xlabel("Number of simulations")
plt.ylabel("Relative error")
plt.title("Delta :  Monte Carlo vs Antithetic vs Control Variate Geo")
plt.legend()
plt.grid(True, which="both", ls="--", alpha=0.5)

plt.show()

#Histograms
df = pd.read_csv("delta_asian_histogram.txt", sep=";")

plt.figure(figsize=(8, 5))

plt.hist(df["Delta_MC"], bins=30, alpha=0.6, label="Monte Carlo standard")
plt.hist(df["Delta_Ant"], bins=30, alpha=0.6, label="Monte Carlo Antithetic")
plt.hist(df["Delta_CV"], bins=30, alpha=0.6, label="Monte Carlo Control Variate Geo")

plt.xlabel("Value")
plt.ylabel("Frequency")
plt.title("Histogram")
plt.legend()
plt.grid(True)

plt.show()

#delta w.r.t. different S_0
df = pd.read_csv("delta_asian_vs_spot.txt", sep=";")

plt.figure(figsize=(8, 5))

plt.plot(df["S0"], df["Delta_MC"],
         label="Monte Carlo standard",
         linestyle="--",
         marker="o")

plt.plot(df["S0"], df["Delta_Ant"],
         label="Monte Carlo Antithetic",
         linestyle="--",
         marker="s")

plt.plot(df["S0"], df["Delta_CV"],
         label="Monte Carlo Control Variate Geo",
         linestyle="--",
         marker="s")

plt.plot(df["S0"], df["Delta_Geo_Exact"],
         label="Delta Black-Scholes",
         linewidth=2)

plt.xlabel("S0")
plt.ylabel("Delta")
plt.title("Delta : MC vs Antithetic vs Control Variate vs BS")
plt.legend()
plt.grid(True)

plt.show()


