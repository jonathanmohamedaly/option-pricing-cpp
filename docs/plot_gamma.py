import pandas as pd
import matplotlib.pyplot as plt

#plot of MC vs Antithetic
mc = pd.read_csv("convergence_gamma_mc.txt", sep=";")
mc_anti = pd.read_csv("convergence_gamma_antithetic.txt", sep=";")

plt.figure(figsize=(10,6))
plt.plot(mc['nSimulations'], mc['Relative_Error'], label="Monte Carlo standard", alpha=0.7)
plt.plot(mc_anti['nSimulations'], mc_anti['Relative_Error'], label="Monte Carlo Antithetic", alpha=0.7)

plt.xscale('log')
plt.yscale('log')
plt.xlabel("Number of simulations")
plt.ylabel("Relative error")
plt.title("Gamma :  Monte Carlo vs Antithetic")
plt.legend()
plt.grid(True, which="both", ls="--", alpha=0.5)

plt.show()

#Histograms
df = pd.read_csv("gamma_histogram.txt", sep=";")

plt.figure(figsize=(8, 5))

plt.hist(df["Gamma_MC"], bins=30, alpha=0.6, label="Gamma_MC")
plt.hist(df["Gamma_Antithetic"], bins=30, alpha=0.6, label="Gamma_Antithetic")

plt.xlabel("Value")
plt.ylabel("Frequency")
plt.title("Histogram")
plt.legend()
plt.grid(True)

plt.show()

#gamma w.r.t. different S_0
df = pd.read_csv("gamma_vs_spot.txt", sep=";")

plt.figure(figsize=(8, 5))

plt.plot(df["S0"], df["Gamma_MC"],
         label="Gamma MC",
         linestyle="--",
         marker="o")

plt.plot(df["S0"], df["Gamma_Antithetic"],
         label="Gamma Antithetic",
         linestyle="--",
         marker="s")

plt.plot(df["S0"], df["Gamma_BS"],
         label="Gamma Black-Scholes",
         linewidth=2)

plt.xlabel("S0")
plt.ylabel("Gamma")
plt.title("Gamma : MC vs Antithetic vs BS")
plt.legend()
plt.grid(True)

plt.show()

#Boxplots

df = pd.read_csv("gamma_boxplot.txt", sep=";")
plt.figure(figsize=(6, 5))

plt.boxplot(
    [df["Gamma_MC"], df["Gamma_Antithetic"]],
    labels=["MC", "Antithetic"],
    showfliers=False
)

plt.ylabel("Gamma")
plt.title("Boxplot : MC vs Antithetic")
plt.grid(axis="y", linestyle="--", alpha=0.6)

plt.show()
