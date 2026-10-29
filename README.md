# Superconducting Qubit & CR Gate Physics Simulator

A rigorous, research-grade Python simulation suite for superconducting quantum circuits, focusing on Transmon qubits, Cross-Resonance (CR) gate dynamics, and parasitic ZZ-crosstalk analysis using exact Hamiltonian diagonalization and dispersive perturbation theory.

## 🔬 Core Features
- **Transmon Hamiltonian Construction:** Implements full charge/Fock-basis Hamiltonian with self-Kerr non-linearities (charging energies).
- **Exact vs. Approximate ZZ-Crosstalk:** Extracts exact ZZ interaction rates directly from the four lowest eigenenergies ($E_{00}, E_{01}, E_{10}, E_{11}$) and compares them with analytical dispersive perturbation theory.
- **Parameter Sweep Testing:** Automated verification routines across various qubit detunings ($\Delta$) to benchmark physical validity.

## 🛠️ Mathematical Framework
1. **Transmon Hamiltonian:**
   $$H_{\text{transmon}} = \sum_{i} \left( \omega_i a_i^\dagger a_i - \frac{E_{C,i}}{2} n_i(n_i - I) \right)$$
2. **Exact ZZ Interaction:**
   $$\zeta_{ZZ} = (E_{11} + E_{00}) - (E_{10} + E_{01})$$

## 🚀 Quick Start
Run the CR gate and transmon simulation suite locally:
```bash
python cr_gate_simulator.py
