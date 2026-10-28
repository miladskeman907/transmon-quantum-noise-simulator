import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as integrate
import networkx as nx
import random

# ==========================================
# 1. الثوابت الفيزيائية ونموذج نويز 1/f
# ==========================================
hbar = 1.0545718e-34
kB = 1.380649e-23

Delta = 1.76 * kB * 1.2
EJ0 = 20 * 2 * np.pi * 1e9 * hbar
EC = 0.25 * 2 * np.pi * 1e9 * hbar
x_qp = 1e-5
T_env = 0.010
t_gate = 40e-9

# بارامترات نويز التدفق 1/f
A_flux = 3e-6  
ln_factor = 10.0 

def get_qubit_dynamics(flux):
    """حساب الطاقات والحساسية المشتقة للتدفق"""
    d_flux = 1e-5
    
    EJ_f = 2 * EJ0 * np.abs(np.cos(np.pi * flux / 2))
    omega_01 = (np.sqrt(8 * EJ_f * EC) - EC) / hbar
    
    EJ_f_plus = 2 * EJ0 * np.abs(np.cos(np.pi * (flux + d_flux) / 2))
    omega_01_plus = (np.sqrt(8 * EJ_f_plus * EC) - EC) / hbar
    derivative = (omega_01_plus - omega_01) / d_flux
    
    return EJ_f, omega_01, np.abs(derivative)

def bcs_integrand(E, omega_01):
    numerator = E * (E + hbar * omega_01) + Delta**2
    denominator = np.sqrt(E**2 - Delta**2) * np.sqrt((E + hbar * omega_01)**2 - Delta**2)
    return (numerator / denominator) * np.exp(-E / (kB * T_env))

def calculate_combined_error_rate(flux):
    """حساب احتمالية الخطأ المادي الكلي المدمج (QP + Flux Noise)"""
    EJ_f, omega_01, dispersion = get_qubit_dynamics(flux)
    
    if omega_01 <= 0: return 0.5
    
    integral_val, _ = integrate.quad(bcs_integrand, Delta, Delta + 10 * kB * T_env, args=(omega_01,))
    prefactor = np.cos(np.pi * flux / 2)**2 * (8 * EJ_f / (np.pi * hbar)) * (EC / (2 * hbar * omega_01))
    Gamma_1 = prefactor * x_qp * integral_val
    
    Gamma_phi = np.sqrt(2 * (A_flux)**2 * ln_factor) * dispersion
    Gamma_2 = (0.5 * Gamma_1) + Gamma_phi
    p_error = 1.0 - np.exp(-Gamma_2 * t_gate)
    return np.clip(p_error, 0, 0.5)

# ==========================================
# 2. محاكاة خوارزمية MWPM لشفرة السطح
# ==========================================
def simulate_surface_code_mwpm(d, p_error, trials=200):
    """محاكاة مجهرية لطوبولوجيا شفرة السطح وفك التشفير عبر MWPM"""
    if p_error <= 0: return 0.0
    if p_error >= 0.5: return 0.5
    
    failed_trials = 0
    num_syndromes = d - 1
    
    for _ in range(trials):
        errors = np.zeros((d, d), dtype=bool)
        for i in range(d):
            for j in range(d):
                if random.random() < p_error:
                    errors[i, j] = True
        
        syndromes = []
        for i in range(num_syndromes):
            for j in range(d):
                if errors[i, j] ^ errors[i+1, j]:
                    syndromes.append((i, j))
        
        if len(syndromes) == 0:
            continue  
            
        G = nx.Graph()
        for idx1, s1 in enumerate(syndromes):
            for idx2, s2 in enumerate(syndromes):
                if idx1 < idx2:
                    weight = abs(s1[0] - s2[0]) + abs(s1[1] - s2[1])
                    G.add_edge(idx1, idx2, weight=-weight)
                    
        matching = nx.algorithms.matching.max_weight_matching(G, maxcardinality=True)
        
        total_correction = np.zeros((d, d), dtype=bool)
        for idx1, idx2 in matching:
            s1, s2 = syndromes[idx1], syndromes[idx2]
            r_start, r_end = min(s1[0], s2[0]), max(s1[0], s2[0])
            for r in range(r_start + 1, r_end + 1):
                total_correction[r, s1[1]] ^= True
                
        residual_errors = errors ^ total_correction
        if np.any(np.all(residual_errors, axis=0)):
            failed_trials += 1
            
    return failed_trials / trials

if __name__ == "__main__":
    print("Transmon Quantum Noise & MWPM Simulation Suite Initialized.")
