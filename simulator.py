import numpy as np
import scipy.linalg as la
import matplotlib.pyplot as plt

def compute_transmon_hamiltonian(Ej=20.0, Ec=0.2, n_states=6, ng=0.0):
    """
    Computes the Superconducting Transmon Qubit Hamiltonian in the charge basis
    and returns its energy eigenvalues and transition frequencies.
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
    
    # Shift energies relative to the ground state
    energies = energies - energies[0]
    
    print(f"Computed energy levels (GHz): {energies[:4]}")
    return energies, eigenstates

if __name__ == "__main__":
    compute_transmon_hamiltonian()
