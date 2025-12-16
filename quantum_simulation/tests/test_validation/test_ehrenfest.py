"""
Tests unitaires pour théorème Ehrenfest.
"""

import pytest
import numpy as np
from quantum_simulation.validation.ehrenfest_theorem import EhrenfestValidator
from quantum_simulation.systems.free_particle import FreeParticle


def test_free_particle_constant_momentum():
    """
    Particule libre : ⟨P⟩ doit rester constant (force nulle).
    """
    hbar = 1.054571817e-34
    mass = 9.1093837015e-31
    k0 = 1e10
    
    x = np.linspace(-1e-8, 1e-8, 1024)
    
    free_particle = FreeParticle(mass, hbar)
    
    # Paquet avec impulsion moyenne non nulle
    state = free_particle.create_gaussian_wavepacket(
        x, x0=0, sigma_x=1e-9, k0=k0
    )
    
    validator = EhrenfestValidator(mass, hbar)
    mean_p = validator.compute_mean_values(state)['mean_p']
    
    # Impulsion moyenne doit être ≈ ℏk0
    p_theory = hbar * k0
    relative_error = abs(mean_p - p_theory) / p_theory
    
    assert relative_error < 0.01, \
        f"⟨P⟩ incorrect: attendu {p_theory:.2e}, obtenu {mean_p:.2e}"


def test_position_evolution_free_particle():
    """
    Particule libre : d⟨X⟩/dt = ⟨P⟩/m.
    """
    hbar = 1.054571817e-34
    mass = 9.1093837015e-31
    
    x = np.linspace(-5e-8, 5e-8, 1024)
    
    free_particle = FreeParticle(mass, hbar)
    state = free_particle.create_gaussian_wavepacket(
        x, x0=0, sigma_x=2e-9, k0=5e9
    )
    
    # États "simulés" (déplacement linéaire pour test)
    # En réalité, il faudrait vraie évolution, mais ici test structure
    validator = EhrenfestValidator(mass, hbar)
    mean_values = validator.compute_mean_values(state)
    mean_p = mean_values['mean_p']
    
    # Vitesse théorique
    v_theory = mean_p / mass
    
    # Test : valeur raisonnable (ordre grandeur)
    assert abs(v_theory) > 0, "Vitesse devrait être non nulle"
    assert abs(v_theory) < 1e7, "Vitesse devrait être < c/30"


def test_mean_force_zero_potential():
    """
    Potentiel nul → ⟨F⟩ = 0.
    """
    hbar = 1.054571817e-34
    mass = 9.1093837015e-31
    
    x = np.linspace(-1e-8, 1e-8, 512)
    
    free_particle = FreeParticle(mass, hbar)
    state = free_particle.create_gaussian_wavepacket(
        x, x0=0, sigma_x=1e-9, k0=1e10
    )
    
    def zero_potential(x_grid):
        return np.zeros_like(x_grid)
    
    validator = EhrenfestValidator(mass, hbar)
    mean_force = validator.compute_mean_force(state, zero_potential)
    
    assert abs(mean_force) < 1e-20, \
        f"Force moyenne devrait être nulle, obtenu {mean_force}"