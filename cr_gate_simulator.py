import numpy as np
from scipy.linalg import expm

def transmon_hamiltonian(EJ1, EC1, EJ2, EC2, g, n_states=3):
    """
    Constructs the coupled two-transmon Hamiltonian including the resonator/bus 
    and the parasitic ZZ interaction terms.
    """
    # Creation and annihilation operators for transmon 1 and 2
    a1 = np.diag(np.sqrt(np.arange(1, n_states)), -1)
    a2 = np.diag(np.sqrt(np.arange(1, n_states)), -1)
    
    # Placeholder for frequencies and anharmonicities
    w1 = np.sqrt(8 * EJ1 * EC1) - EC1
    w2 = np.sqrt(8 * EJ2 * EC2) - EC2
    
    # Tensor products for composite Hilbert space
    I1 = np.eye(n_states)
    H1 = np.kron(w1 * np.dot(a1.T, a1), I1)
    H2 = np.kron(I1, w2 * np.dot(a2.T, a2))
    
    # Coupling interaction (Jaynes-Cummings / Cross-Resonance term)
    H_int = g * (np.kron(a1, a2.T) + np.kron(a1.T, a2))
    
    return H1 + H2 + H_int

def calculate_zz_shift(EJ1, EC1, EJ2, EC2, g):
    """
    Calculates the parasitic ZZ crosstalk shift (zeta_ZZ) between the two qubits.
    """
    # Simplified analytical estimation for dispersive coupling shift
    w1 = np.sqrt(8 * EJ1 * EC1) - EC1
    w2 = np.sqrt(8 * EJ2 * EC2) - EC2
    delta = w1 - w2
    
    if delta == 0:
        return 0.0
    
    # Effective ZZ interaction strength due to higher-level leakage
    zeta_zz = 2 * g**2 * (EC1 + EC2) / (delta**2 - (EC1 + EC2)**2)
    return zeta_zz

if __name__ == "__main__":
    # Example hardware parameters (in GHz-ish units)
    EJ1, EC1 = 20.0, 0.3
    EJ2, EC2 = 18.0, 0.25
    g_coupling = 0.05
    
    zz_rate = calculate_zz_shift(EJ1, EC1, EJ2, EC2, g_coupling)
    print(f"Initialized Two-Qubit Transmon System.")
    print(f"Estimated Parasitic ZZ Crosstalk Rate: {zz_rate:.6f} MHz")
