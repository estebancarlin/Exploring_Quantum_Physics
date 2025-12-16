"""
Tests unitaires pour systèmes à potentiel.
"""

import sys
from pathlib import Path
# HACK: Ajouter racine projet au path (temporaire)
project_root = Path(__file__).parent.parent.parent
sys.path.insert(0, str(project_root))

import pytest
import numpy as np
from quantum_simulation.systems.potential_systems import (
    InfiniteWell, FiniteWell, PotentialBarrier
)


def test_infinite_well_energy_spectrum():
    """
    Énergies puits infini : Eₙ = n²π²ℏ²/(2mL²)
    """
    hbar = 1.054571817e-34
    mass = 9.1093837015e-31
    width = 1e-9
    
    well = InfiniteWell(width, mass, hbar)
    
    for n in range(1, 10):
        E_n = well.energy_eigenvalue(n)
        E_theory = (n**2 * np.pi**2 * hbar**2) / (2 * mass * width**2)
        
        assert abs(E_n - E_theory) < 1e-50, \
            f"E_{n} incorrect"


def test_infinite_well_wavefunction_normalization():
    """
    États propres puits infini doivent être normés.
    
    Note:
        Normalisation analytique √(2/L) donne erreur O(dx²) sur grille discrète.
        Pour n grand, oscillations rapides augmentent erreur d'intégration.
        Tolérance adaptée : 1e-6 (au lieu de 1e-8).
    """
    hbar = 1.054571817e-34
    mass = 9.1093837015e-31
    width = 1e-9
    
    x = np.linspace(-width, 2*width, 2048)
    
    well = InfiniteWell(width, mass, hbar)
    
    for n in [1, 2, 5]:
        state = well.eigenstate_wavefunction(n, x)
        # Tolérance 1e-6 : acceptable pour n ≤ 10 sur grille 2048 points
        assert state.is_normalized(tolerance=1e-6), \
            f"État n={n} non normé (norme={state.norm():.10f})"


def test_infinite_well_boundary_conditions():
    """
    ψ(0) = ψ(L) = 0 pour puits infini.
    
    Note:
        Grille doit contenir exactement x=0 et x=L pour tester conditions limites.
    """
    hbar = 1.054571817e-34
    mass = 9.1093837015e-31
    width = 1e-9
    
    # Grille alignée sur bords [0, L]
    # Inclut exactement x=0 (premier point) et x=L (dernier point)
    x = np.linspace(0, width, 1024)
    
    well = InfiniteWell(width, mass, hbar)
    
    for n in [1, 3, 7]:
        state = well.eigenstate_wavefunction(n, x)
        psi = state.wavefunction
        
        # Vérifier bords : indices 0 et -1 correspondent à x=0 et x=L
        assert abs(psi[0]) < 1e-10, \
            f"ψ_{n}(0) = {abs(psi[0]):.2e} devrait être nul"
        assert abs(psi[-1]) < 1e-10, \
            f"ψ_{n}(L) = {abs(psi[-1]):.2e} devrait être nul"


def test_infinite_well_orthonormality():
    """
    ⟨ψₙ|ψₘ⟩ = δₙₘ pour puits infini.
    """
    hbar = 1.054571817e-34
    mass = 9.1093837015e-31
    width = 1e-9
    
    x = np.linspace(-width/2, 1.5*width, 2048)
    
    well = InfiniteWell(width, mass, hbar)
    
    is_orthonormal = well.validate_orthonormality(n_states=5, x=x, tolerance=1e-6)
    assert is_orthonormal, "États propres devraient être orthonormés"


def test_finite_well_bound_states_estimate():
    """
    Nombre états liés puits fini doit être cohérent.
    """
    hbar = 1.054571817e-34
    mass = 9.1093837015e-31
    width = 1e-9
    depth = 1e-18  # ~6 eV
    
    well = FiniteWell(width, depth, mass, hbar)
    
    n_bound = well.estimate_bound_states_count()
    
    # Doit être au moins 1
    assert n_bound >= 1, "Puits fini doit avoir au moins 1 état lié"
    
    # Doit être fini (pas trop grand)
    assert n_bound < 100, "Nombre états liés suspects"


def test_barrier_transmission_classical():
    """
    Transmission sur-barrière (E > V₀) devrait être ~1.
    """
    hbar = 1.054571817e-34
    mass = 9.1093837015e-31
    width = 1e-9
    height = 1e-19
    
    barrier = PotentialBarrier(width, height, mass, hbar)
    
    # Énergie > hauteur
    E = 2 * height
    T = barrier.transmission_coefficient_approx(E)
    
    assert T > 0.9, "Transmission classique devrait être proche de 1"


def test_barrier_transmission_tunnel():
    """
    Transmission sous-barrière (effet tunnel) devrait être 0 < T < 1.
    """
    hbar = 1.054571817e-34
    mass = 9.1093837015e-31
    width = 1e-9
    height = 1e-18
    
    barrier = PotentialBarrier(width, height, mass, hbar)
    
    # Énergie < hauteur
    E = 0.5 * height
    T = barrier.transmission_coefficient_approx(E)
    
    assert 0 < T < 1, "Transmission tunnel devrait être 0 < T < 1"
    assert T < 0.5, "Transmission barrière opaque devrait être faible"


def test_potential_continuity():
    """
    Potentiels doivent être continus par morceaux.
    """
    hbar = 1.054571817e-34
    mass = 9.1093837015e-31
    width = 1e-9
    
    x = np.linspace(-2*width, 2*width, 1024)
    
    # Puits infini
    well = InfiniteWell(width, mass, hbar)
    V_well = well.potential(x)
    
    # Pas de NaN ou inf (sauf murs)
    assert not np.any(np.isnan(V_well))
    
    # Puits fini
    finite_well = FiniteWell(width, 1e-18, mass, hbar)
    V_finite = finite_well.potential(x)
    
    assert not np.any(np.isnan(V_finite))
    assert np.all(np.isfinite(V_finite))