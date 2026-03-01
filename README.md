# Option Pricing Engine (C++)

## Overview

This project implements a Monte Carlo option pricing engine in C++ with a focus on numerical stability, variance reduction techniques, and sensitivity analysis (Greeks).

The objective is to study convergence properties, improve estimator efficiency, and analyze the robustness of finite-difference Greeks.

---

## Models Implemented

### Black-Scholes (Analytical)
- Closed-form pricing for European call and put options
- Benchmark for Monte Carlo estimators

### Monte Carlo Pricing
- European options
- Asian options (arithmetic and geometric)

The Monte Carlo estimator exhibits the expected convergence rate:

\[
\mathcal{O}(1/\sqrt{N})
\]

---

## Variance Reduction Techniques

To improve estimator efficiency:

- Antithetic Variates
- Control Variates  
- Common Random Numbers (CRN) for Greek estimation

Control variates are shown to be less useful for standard European vanilla options, where the analytical solution is already available.

---

## Greeks Estimation

Finite-difference estimation of:

- Delta
- Gamma
- Vega

Each Greek is analyzed under:
- Standard Monte Carlo
- Antithetic variates
- Antithetic + CRN

The project studies:
- Convergence behavior
- Variance reduction impact
- Numerical stability (especially for Gamma)

Graphical analysis is provided for each sensitivity.

---

## Asian Options

### Arithmetic Asian Options
- Monte Carlo pricing (Call & Put)
- Greek estimation and stability analysis

### Geometric Asian Options
- Analytical benchmark
- Monte Carlo validation
- Greek computation

---

## Key Features

- Full C++ implementation
- Modular structure
- Convergence analysis

---

## Future Improvements

- Pathwise Greek estimators
- Volatility smile integration
- Delta hedging simulation
- Extension to stochastic volatility models (e.g., Heston)

# Moteur de Pricing d’Options (C++)

## Présentation

Ce projet consiste en l’implémentation d’un moteur de pricing Monte Carlo en C++ avec un focus sur la stabilité numérique, la réduction de variance et le calcul des sensibilités (Greeks).

L’objectif est d’analyser la convergence des estimateurs, d’améliorer leur efficacité et d’étudier la robustesse des dérivées numériques.

---

## Modèles Implémentés

### Black-Scholes (Analytique)
- Formule fermée pour options européennes Call et Put
- Benchmark pour les estimateurs Monte Carlo

### Pricing Monte Carlo
- Options européennes
- Options asiatiques (arithmétiques et géométriques)

L’estimateur Monte Carlo présente la convergence théorique attendue :

\[
\mathcal{O}(1/\sqrt{N})
\]

---

## Techniques de Réduction de Variance

Afin d’améliorer l’efficacité :

- Variables antithétiques
- Variables de contrôle
- Common Random Numbers (CRN) pour le calcul des Greeks

Les variables de contrôle sont peu pertinentes pour les options européennes standards disposant d’une solution analytique.

---

## Calcul des Greeks

Estimation par différences finies de :

- Delta
- Gamma
- Vega

Chaque sensibilité est analysée sous :
- Monte Carlo standard
- Variables antithétiques
- Antithétique + CRN

Le projet étudie :
- La convergence
- L’impact de la réduction de variance
- La stabilité numérique (notamment pour Gamma)

Des graphiques illustrent le comportement de chaque Greek.

---

## Options Asiatiques

### Options Asiatiques Arithmétiques
- Pricing Monte Carlo (Call & Put)
- Analyse des Greeks

### Options Asiatiques Géométriques
- Solution analytique de référence
- Validation Monte Carlo
- Calcul des sensibilités

---

## Points Clés

- Implémentation complète en C++
- Structure modulaire
- Analyse de convergence
- Base extensible vers modèles plus complexes

---

## Améliorations Futures

- Estimateurs pathwise des Greeks
- Intégration d’un smile de volatilité
- Simulation de delta-hedging
- Extension vers volatilité stochastique (Heston)