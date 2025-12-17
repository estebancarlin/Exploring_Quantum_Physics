"""
Tests fonctions d'onde HO en représentation position (Extension D4).

Vérifie :
- Orthonormalité ⟨ψₙ|ψₘ⟩ = δₙₘ
- Cohérence avec base abstraite (énergies)
"""

import sys
from pathlib import Path
project_root = Path(__file__).parent.parent.parent.parent
sys.path.insert(0, str(project_root))

import pytest
import numpy as np
from quantum_simulation.systems.harmonic_oscillator import HarmonicOscillator


def test_orthonormality_wavefunctions():
    """
    Vérifier ⟨ψₙ|ψₘ⟩ = δₙₘ numériquement.
    """
    hbar = 1.054571817e-34
    mass = 9.1093837015e-31
    omega = 1e15
    
    ho = HarmonicOscillator(mass, omega, hbar, n_max=10)
    
    # Grille couvrant ±5σ
    x0 = np.sqrt(hbar / (mass * omega))
    x = np.linspace(-5*x0, 5*x0, 2048)
    
    # Tester n=0..4
    tolerance = 1e-6
    
    for n in range(5):
        psi_n = ho.wavefunction_position(n, x)
        
        # Auto-normalisation
        norm_n = psi_n.norm()
        assert abs(norm_n - 1.0) < tolerance, f"ψ{n} non normalisé : {norm_n}"
        
        # Orthogonalité avec m > n
        for m in range(n+1, 5):
            psi_m = ho.wavefunction_position(m, x)
            
            overlap = psi_n.inner_product(psi_m)
            
            assert abs(overlap) < tolerance, \
                f"⟨ψ{n}|ψ{m}⟩ = {overlap:.2e} ≠ 0"


def test_consistency_energy_spectrum():
    """
    Vérifier Eₙ = ℏω(n+1/2) cohérent entre base abstraite et représentation position.
    """
    hbar = 1.054571817e-34
    mass = 9.1093837015e-31
    omega = 1e15
    
    ho = HarmonicOscillator(mass, omega, hbar, n_max=10)
    
    x0 = np.sqrt(hbar / (mass * omega))
    x = np.linspace(-5*x0, 5*x0, 2048)
    
    from quantum_simulation.core.operators import Hamiltonian
    V_ho = lambda x_val: 0.5 * mass * omega**2 * x_val**2
    
    H = Hamiltonian(mass, hbar, potential=V_ho)
    
    # Tester n=0,1,2
    for n in range(3):
        psi_n = ho.wavefunction_position(n, x)
        
        # Énergie attendue
        E_expected = ho.energy_eigenvalue(n)
        
        # Énergie mesurée ⟨ψₙ|H|ψₙ⟩
        E_measured = H.expectation_value(psi_n)
        
        relative_error = abs(E_measured - E_expected) / E_expected
        
        assert relative_error < 1e-3, \
            f"Niveau n={n} : E_mesurée={E_measured:.2e}, E_attendue={E_expected:.2e}"


if __name__ == "__main__":
    pytest.main([__file__, "-v"])