
import numpy as np
import scipy.linalg as la

def construct_full_transmon_hamiltonian(EJ1=20.0, EC1=0.3, EJ2=18.0, EC2=0.25, g_coupling=0.05, n_states=3):
    """
    Constructs the full two-transmon Hamiltonian including individual Kerr non-linearities (charging energies 
    acting as anharmonicities) and capacitive/inductive coupling terms in the tensor product Fock basis.
    """
    print(f"Constructing Full Hamiltonian: Q1(EJ={EJ1}, EC={EC1}), Q2(EJ={EJ2}, EC={EC2}), g={g_coupling}")
    
    # Creation and annihilation operators
    a1 = np.diag(np.sqrt(np.arange(1, n_states)), k=1)
    a2 = np.diag(np.sqrt(np.arange(1, n_states)), k=1)
    
    # Number operators (n = a^\dagger a)
    num1 = np.dot(a1.T, a1)
    num2 = np.dot(a2.T, a2)
    
    # Transmon bare frequencies (omega ~ sqrt(8*EJ*EC) - EC)
    w1 = np.sqrt(8.0 * EJ1 * EC1) - EC1
    w2 = np.sqrt(8.0 * EJ2 * EC2) - EC2
    
    # Identity operators
    I1 = np.eye(n_states, dtype=complex)
    I2 = np.eye(n_states, dtype=complex)
    
    # Base Hamiltonian with Kerr/Anharmonicity terms (- 0.5 * EC * n * (n-1))
    # For simplicity in Fock basis, transmon energy incorporates the non-linear charging term
    H1 = np.kron(w1 * num1 - 0.5 * EC1 * np.dot(num1, num1 - I1), I2)
    H2 = np.kron(I1, w2 * num2 - 0.5 * EC2 * np.dot(num2, num2 - I2))
    
    # Interaction Hamiltonian (Jaynes-Cummings / Cross-Resonance coupling)
    H_int = g_coupling * (np.kron(a1, a2.T) + np.kron(a1.T, a2))
    
    H_total = H1 + H2 + H_int
    return H_total

def calculate_exact_zz_from_energies(H_total):
    """
    Calculates the exact ZZ-crosstalk from the eigenenergies of the four lowest computational states.
    zeta_ZZ = (E_11 - E_10) - (E_01 - E_00) = E_11 + E_00 - E_10 - E_01
    """
    energies, _ = la.eigh(H_total)
    # The lowest 4 energy levels correspond to computational subspace |00>, |01>, |10>, |11>
    E00 = energies[0]
    E01 = energies[1]
    E10 = energies[2]
    E11 = energies[3]
    
    zeta_zz_exact = (E11 + E00) - (E10 + E01)
    return zeta_zz_exact, (E00, E01, E10, E11)

def calculate_approximate_zz(EJ1=20.0, EC1=0.3, EJ2=18.0, EC2=0.25, g_coupling=0.05):
    """
    Calculates the approximate ZZ crosstalk using analytical dispersive perturbation theory.
    """
    w1 = np.sqrt(8.0 * EJ1 * EC1) - EC1
    w2 = np.sqrt(8.0 * EJ2 * EC2) - EC2
    delta = w1 - w2
    
    if delta == 0:
        return 0.0
        
    zeta_zz_approx = 2.0 * (g_coupling ** 2) * (EC1 + EC2) / (delta**2 - (EC1 + EC2)**2)
    return zeta_zz_approx

def run_parameter_sweep_tests():
    """
    Runs unit tests and compares exact vs approximate ZZ crosstalk across different detuning values.
    """
    print("\n--- Running Parameter Sweep & Verification Tests ---")
    detuning_values = [0.5, 1.0, 2.0, 3.0, 5.0]
    
    for d in detuning_values:
        # Varying EJ2 to change detuning
        w1 = np.sqrt(8.0 * 20.0 * 0.3) - 0.3
        target_w2 = w1 - d
        # Approximate reverse engineering of EJ2 for target frequency
        estimated_EJ2 = ((target_w2 + 0.25)**2) / (8.0 * 0.25)
        
        H = construct_full_transmon_hamiltonian(EJ1=20.0, EC1=0.3, EJ2=estimated_EJ2, EC2=0.25, g_coupling=0.05)
        exact_zz, _ = calculate_exact_zz_from_energies(H)
        approx_zz = calculate_approximate_zz(EJ1=20.0, EC1=0.3, EJ2=estimated_EJ2, EC2=0.25, g_coupling=0.05)
        
.        print(f"Detuning Δ = {d:.2f} GHz | Exact ZZ: {exact_zz*1000:.4f} MHz | Approx ZZ: {approx_zz*1000:.4f} MHz")

if __name__ == "__main__":
    H_sys = construct_full_transmon_hamiltonian()
    exact_zz, levels = calculate_exact_zz_from_energies(H_sys)
    approx_zz = calculate_approximate_zz()
    
    print(f"\nExact Lowest 4 Energy Levels (GHz): {levels}")
    print(f"Exact ZZ Interaction (from energy levels): {exact_zz*1000:.4f} MHz")
    print(f"Approximate ZZ Interaction (analytical): {approx_zz*1000:.4f} MHz")
    
    # Run tests across various detunings
    run_parameter_sweep_tests()
