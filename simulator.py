import numpy as np
import scipy.linalg as la
import matplotlib.pyplot as plt

def compute_transmon_hamiltonian(Ej=20.0, Ec=0.2, n_states=6, ng=0.0):
    """
    Computes the Superconducting Transmon Qubit Hamiltonian in the charge basis,
    finds energy eigenvalues, and calculates the fundamental transition frequency f01 
    and anharmonicity.
    """
    print(f"Initializing Transmon Hamiltonian with Ej={Ej} GHz, Ec={Ec} GHz, N_states={n_states}")
    
    # Charge operator matrix size
    N = 2 * n_states + 1
    n_vec = np.arange(-n_states, n_states + 1) - ng
    
    # Charging energy matrix (diagonal)
    H_charge = np.diag(4.0 * Ec * (n_vec ** 2))
    
    # Josephson energy matrix (off-diagonal coupling via Cooper pair tunneling)
    H_josephson = -0.5 * Ej * (np.diag(np.ones(N - 1), k=1) + np.diag(np.ones(N - 1), k=-1))
    
    # Total Hamiltonian
    H_transmon = H_charge + H_josephson
    
    # Diagonalize to find energy eigenvalues
    energies, eigenstates = la.eigh(H_transmon)
    
    # Shift energies relative to the ground state (E0 = 0)
    energies = energies - energies[0]
    
    # Calculate transition frequencies (f01, f12) and anharmonicity in GHz
    f01 = energies[1] - energies[0]
    f12 = energies[2] - energies[1]
    anharmonicity = f12 - f01
    
    print(f"Computed energy levels (GHz): {energies[:4]}")
    print(f"Fundamental Transition Frequency (f01): {f01:.4f} GHz")
    print(f"Anharmonicity (f12 - f01): {anharmonicity * 1000:.2f} MHz")
    
    return energies, f01, anharmonicity

if __name__ == "__main__":
    compute_transmon_hamiltonian()
