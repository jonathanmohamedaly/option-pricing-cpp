import pandas as pd
import matplotlib.pyplot as plt

# -------------------------------
# Convergence plot
# -------------------------------
mc = pd.read_csv("convergence_delta_asian_pathwise.txt", sep=";")

plt.figure(figsize=(10,6))
plt.plot(mc['nSim'], mc['RelErr_PW_Ant_CV'], label="Pathwise Antithetic + CV", alpha=0.7)

plt.xscale('log')
plt.yscale('log')
plt.xlabel("Number of simulations")
plt.ylabel("Relative error")
plt.title("Delta Convergence Control Variate")
plt.legend()
plt.grid(True, which="both", ls="--", alpha=0.5)

plt.show()

# -------------------------------
# Histograms
# -------------------------------
df = pd.read_csv("delta_asian_histogram_pathwise.txt", sep=";")

plt.figure(figsize=(8, 5))
plt.hist(df["Delta_PW_Ant_CV"], bins=30, alpha=0.6, label="Pathwise Antithetic + CV")

plt.xlabel("Delta value")
plt.ylabel("Frequency")
plt.title("Histogram of Delta estimates")
plt.legend()
plt.grid(True)
plt.show()

# -------------------------------
# Delta vs Spot
# -------------------------------
df = pd.read_csv("delta_asian_vs_spot_pathwise.txt", sep=";")

plt.figure(figsize=(8, 5))

plt.plot(df["S0"], df["Delta_PW_Ant_CV"], label="Pathwise Antithetic + CV", linestyle="--", marker="^")
plt.plot(df["S0"], df["Delta_Geo_Exact"], label="Delta Black-Scholes Geo", linewidth=2)

plt.xlabel("S0")
plt.ylabel("Delta")
plt.title("Delta vs Spot: Pathwise methods vs Black-Scholes Geometric")
plt.legend()
plt.grid(True)
plt.show()