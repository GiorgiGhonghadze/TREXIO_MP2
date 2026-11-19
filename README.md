Closed-Shell Hartree-Fock & MP2 Energy CalculatorA C program for computing the closed-shell restricted Hartree-Fock (RHF) energy and the second-order Møller–Plesset (MP2) correlation energy from the data stored in TREXIO format.

Overview

The program reads pre-computed molecular orbital (MO) integrals and orbital energies from a TREXIO file (.h5) and performs the following steps:

Data Extraction: Reads the nuclear repulsion energy ($E_{NN}$), number of spin-up electrons ($N_{\text{occ}}$), total number of molecular orbitals ($N_{\text{MO}}$), canonical orbital energies ($\epsilon_p$), one-electron core Hamiltonian matrix elements ($\langle i \vert{} h \vert{} i \rangle$), and sparse two-electron repulsion integrals ($\langle ij \vert{} kl \rangle$).

Symmetry Expansion: Unpacks the sparse two-electron integrals in physicist's notation into a dense 4D tensor by applying 8-fold permutational symmetry:
$$\langle ij\vert{}kl\rangle = \langle il\vert{}kj\rangle = \langle kl\vert{}ij\rangle = \langle kj\vert{}il\rangle = \langle ji\vert{}lk\rangle = \langle li\vert{}jk\rangle = \langle lk\vert{}ji\rangle = \langle jk\vert{}li\rangle$$

Computes Hartree-Fock Energy:

$$E_{\text{HF}} = E_{NN} + 2\sum_{i=1}^{N_{\text{occ}}} \langle i\vert{}h\vert{}i\rangle + \sum_{i=1}^{N_{\text{occ}}}\sum_{j=1}^{N_{\text{occ}}} \left[ 2\langle ij\vert{}ij\rangle - \langle ij\vert{}ji\rangle \right]$$

Computes MP2 Correlation Energy: 

$$E_{\text{MP2}} = \sum_{i,j \in \text{occ}} \sum_{a,b \in \text{virt}} \frac{\langle ij\vert{}ab\rangle \left(2\langle ij\vert{}ab\rangle - \langle ij\vert{}ba\rangle\right)}{\epsilon_i + \epsilon_j - \epsilon_a - \epsilon_b}$$

Repository Structure

├── AUTHORS          # List of project contributors
├── INSTALL.md       # Compilation
├── LICENSE          # Project license information
├── README.md        # Project summary and benchmark results
├── data/            # Sample TREXIO molecular input files
├── src/             # C source code 
└── tests/           # Validation tests

Computed Benchmark Results

The table below summarizes the Hartree-Fock and MP2 energies obtained across the tested systems (energies reported in Hartree atomic units, $\text{a.u.}$):

| System ($N_{\text{occ}} / N_{\text{MO}}$) | Total MOs ($N_{\text{MO}}$) | Occupied MOs ($N_{\text{occ}}$) | Hartree-Fock Energy (a.u.) | MP2 Correlation Energy (a.u.) | Total MP2 Energy (a.u.) |
| :--- | :---: | :---: | :---: | :---: | :---: |
| System 1 ($\text{H}_2\text{O}$) | 24 | 5 | `-76.0267987082` | `-0.2039599741` | `-76.2307586823` |
| System 2 ($\text{CH}_4$) | 34 | 5 | `-40.1986733443` | `-0.1640072951` | `-40.3626806394` |
| System 3 ($\text{HCN}$) | 33 | 7 | `-92.8829467012` | `-0.2892753290` | `-93.1722220301` |
| System 4 ($\text{C}_2\text{H}_2$) | 38 | 7 | `-76.8257783991` | `-0.2600944278` | `-77.0858728269` |
| System 5 ($\text{CH}_3\text{OH}$) | 48 | 9 | `-115.0484158187` | `-0.3430913311` | `-115.3915071498` |
| System 6 ($\text{CO}_2$) | 42 | 11 | `-187.6507108848` | `-0.4886746075` | `-188.1393854923` |



