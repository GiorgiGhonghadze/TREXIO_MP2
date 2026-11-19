# Closed-Shell Hartree-Fock & MP2 Energy Calculator

C program for computing the closed-shell restricted Hartree-Fock (RHF) energy and the second-order Møller–Plesset (MP2) correlation energy from the data stored in TREXIO format.

<div align="center">

[Overview](#overview) • [Compilation](INSTALL.md) • [Authors](AUTHORS)

</div>

## Overview

The program reads precomputed molecular orbital (MO) integrals and orbital energies from a TREXIO file (`.h5`) and performs the following steps:

1. **Data Extraction:** Reads the nuclear repulsion energy ($E_{NN}$), number of spin-up electrons ($N_{\text{occ}}$), total number of molecular orbitals ($N_{\text{MO}}$), orbital energies ($\epsilon_p$), one-electron Hamiltonian matrix elements ($\langle i | h | i \rangle$), and sparse two-electron repulsion integrals ($\langle ij | kl \rangle$).

2. **Symmetry Expansion:** Unpacks the sparse two-electron integrals into a flattened array by applying 8-fold permutational symmetry:

$$\langle ij|kl\rangle = \langle il|kj\rangle = \langle kl|ij\rangle = \langle kj|il\rangle = \langle ji|lk\rangle = \langle li|jk\rangle = \langle lk|ji\rangle = \langle jk|li\rangle$$

4. **Computes Hartree-Fock Energy:** 

$$E_{\text{HF}} = E_{NN} + 2\sum_{i=1}^{N_{\text{occ}}} \langle i|h|i\rangle + \sum_{i=1}^{N_{\text{occ}}}\sum_{j=1}^{N_{\text{occ}}} \left[ 2\langle ij|ij\rangle - \langle ij|ji\rangle \right]$$
   
5. **Computes MP2 Correlation Energy:**

$$E_{\text{MP2}} = \sum_{i,j \in \text{occ}} \sum_{a,b \in \text{virt}} \frac{\langle ij|ab\rangle \left(2\langle ij|ab\rangle - \langle ij|ba\rangle\right)}{\epsilon_i + \epsilon_j - \epsilon_a - \epsilon_b}$$

## Repository Structure

```text
.
├── AUTHORS          # List of project contributors
├── INSTALL.md       # Compilation instructions
├── LICENSE          # Project license information
├── README.md        # Project summary and benchmark results
├── data/            # Sample TREXIO molecular input files
├── src/             # C source code
└── tests/           # Validation tests
```

## Computed Benchmark Results

The table below summarizes the Hartree-Fock and MP2 energies obtained across the tested closed-shell molecular systems (energies reported in Hartree atomic units, $\text{a.u.}$):

| Compound | Total MOs | Occupied MOs | Hartree-Fock Energy (a.u.) | MP2 Correlation Energy (a.u.) | Total MP2 Energy (a.u.) |
| :--- | :---: | :---: | :---: | :---: | :---: |
| $\text{H}_2\text{O}$ | 24 | 5 | `-76.0267987082` | `-0.2039599741` | `-76.2307586823` |
| $\text{CH}_4$ | 34 | 5 | `-40.1986733443` | `-0.1640072951` | `-40.3626806394` |
| $\text{HCN}$ | 33 | 7 | `-92.8829467012` | `-0.2892753290` | `-93.1722220301` |
| $\text{C}_2\text{H}_2$ | 38 | 7 | `-76.8257783991` | `-0.2600944278` | `-77.0858728269` |
| $\text{CH}_3\text{OH}$ | 48 | 9 | `-115.0484158187` | `-0.3430913311` | `-115.3915071498` |
| $\text{CO}_2$ | 42 | 11 | `-187.6507108848` | `-0.4886746075` | `-188.1393854923` |



