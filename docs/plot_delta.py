import pandas as pd
import matplotlib.pyplot as plt

#plot of MC vs Antithetic
mc = pd.read_csv("convergence_delta_mc.txt", sep=";")
mc_anti = pd.read_csv("convergence_delta_antithetic.txt", sep=";")

plt.figure(figsize=(10,6))
plt.plot(mc['nSimulations'], mc['Relative_Error'], label="Monte Carlo standard", alpha=0.7)
plt.plot(mc_anti['nSimulations'], mc_anti['Relative_Error'], label="Monte Carlo Antithetic", alpha=0.7)

plt.xscale('log')
plt.yscale('log')
plt.xlabel("Number of simulations")
plt.ylabel("Relative error")
plt.title("Delta :  Monte Carlo vs Antithetic")
plt.legend()
plt.grid(True, which="both", ls="--", alpha=0.5)

plt.show()

#Histograms
df = pd.read_csv("delta_histogram.txt", sep=";")

plt.figure(figsize=(8, 5))

plt.hist(df["Delta_MC"], bins=30, alpha=0.6, label="Delta_MC")
plt.hist(df["Delta_Antithetic"], bins=30, alpha=0.6, label="Delta_Antithetic")

plt.xlabel("Value")
plt.ylabel("Frequency")
plt.title("Histogram")
plt.legend()
plt.grid(True)

plt.show()

#delta w.r.t. different S_0
df = pd.read_csv("delta_vs_spot.txt", sep=";")

plt.figure(figsize=(8, 5))

plt.plot(df["S0"], df["Delta_MC"],
         label="Delta MC",
         linestyle="--",
         marker="o")

plt.plot(df["S0"], df["Delta_Antithetic"],
         label="Delta Antithetic",
         linestyle="--",
         marker="s")

plt.plot(df["S0"], df["Delta_BS"],
         label="Delta Black-Scholes",
         linewidth=2)

plt.xlabel("S0")
plt.ylabel("Delta")
plt.title("Delta : MC vs Antithetic vs BS")
plt.legend()
plt.grid(True)

plt.show()

#Boxplots

df = pd.read_csv("delta_boxplot.txt", sep=";")
plt.figure(figsize=(6, 5))

plt.boxplot(
    [df["Delta_MC"], df["Delta_Antithetic"]],
    labels=["MC", "Antithetic"],
    showfliers=False
)

plt.ylabel("Delta")
plt.title("Boxplot : MC vs Antithetic")
plt.grid(axis="y", linestyle="--", alpha=0.6)

plt.show()
