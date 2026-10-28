# transmon-quantum-noise-simulator
Microscopic simulation of transmon qubit T1 relaxation, 1/f flux noise, and MWPM surface code decoding.


A microscopic quantum hardware simulation suite designed to model decoherence mechanisms in superconducting transmon qubits and evaluate topological error correction performance.

## 🔬 Theoretical Framework

### 1. Quasiparticle Poisoning & $T_1$ Relaxation
The relaxation rate $1/T_1$ due to nonequilibrium quasiparticles is modeled using BCS microscopic theory via Fermi's Golden Rule and Type-II coherence factors:

$$\frac{1}{T_1} \propto \int_{\Delta}^{\infty} dE \frac{E(E+\hbar\omega_{01}) + \Delta^2}{\sqrt{E^2-\Delta^2}\sqrt{(E+\hbar\omega_{01})^2-\Delta^2}} f_{qp}(E)$$

### 2. $1/f$ Flux Noise & Dephasing
Incorporates flux dispersion sensitivity $D_\Phi = \frac{\partial \omega_{01}}{\partial \Phi}$ to compute low-frequency flux noise contributions and total effective decoherence rates $\Gamma_2$.

### 3. Surface Code & MWPM Decoding
Simulates physical error probabilities mapped onto a 2D surface code lattice, utilizing Minimum Weight Perfect Matching (MWPM) via NetworkX to evaluate topological threshold resilience.

## 🛠️ Tech Stack
* **Python**
* **NumPy & SciPy** (Numerical integration & physical modeling)
* **NetworkX** (Graph-based syndrome decoding / MWPM)
* **Matplotlib** (Noise dispersion & threshold visualization)

## 📄 License
This project is open-source under the [MIT License](LICENSE).
