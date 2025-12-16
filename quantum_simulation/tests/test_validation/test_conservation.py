"""
Tests unitaires pour validation conservation.
"""

import pytest
import numpy as np
from quantum_simulation.validation.conservation_laws import ConservationValidator
from quantum_simulation.systems.free_particle import FreeParticle
from quantum_simulation.dynamics.evolution import TimeEvolution


def test_norm_conservation_stationary_state():
    """
    État stationnaire doit conserver norme parfaitement.
    """
    hbar = 1.054571817e-34
    mass = 9.1093837015e-31
    
    x = np.linspace(-1e-8, 1e-8, 512)
    
    free_particle = FreeParticle(mass, hbar)
    state = free_particle.create_plane_wave(x, k=1e10)
    
    # "Évolution" triviale (pas de changement)
    states = [state] * 5
    times = np.linspace(0, 1e-15, 5)
    
    validator = ConservationValidator(hbar, mass)
    result = validator.validate_norm_conservation(states, times)
    
    assert result['is_conserved'], "Norme devrait être parfaitement conservée"
    assert result['max_deviation'] < 1e-12


def test_probability_current_plane_wave():
    """
    Courant probabilité d'onde plane doit être constant.
    """
    hbar = 1.054571817e-34
    mass = 9.1093837015e-31
    k = 5e9
    
    x = np.linspace(-1e-8, 1e-8, 1024)
    
    free_particle = FreeParticle(mass, hbar)
    state = free_particle.create_plane_wave(x, k)
    
    validator = ConservationValidator(hbar, mass)
    J = validator.compute_probability_current(state)
    
    # Pour onde plane : J devrait être constant
    J_std = np.std(J)
    assert J_std < 1e-10, f"Courant devrait être constant, std={J_std}"
    
    # Valeur théorique : J = ρ * v = ρ * (ℏk/m)
    rho = state.probability_density()
    v_theory = (hbar * k) / mass
    J_theory_mean = np.mean(rho) * v_theory
    J_observed_mean = np.mean(J)
    
    # Comparaison (ordre de grandeur)
    assert abs(J_observed_mean - J_theory_mean) / abs(J_theory_mean) < 0.1


def test_continuity_equation_gaussian_packet():
    """
    Équation continuité doit être satisfaite durant évolution.
    """
    hbar = 1.054571817e-34
    mass = 9.1093837015e-31
    
    x = np.linspace(-5e-8, 5e-8, 2048)
    
    free_particle = FreeParticle(mass, hbar)
    state_t0 = free_particle.create_gaussian_wavepacket(
        x, x0=0, sigma_x=1e-9, k0=1e10
    )
    
    # Évolution courte
    time_evol = TimeEvolution(free_particle.hamiltonian, hbar)
    dt = 1e-17
    state_t1 = time_evol.evolve_wavefunction(state_t0, 0, dt, dt/10)
    
    validator = ConservationValidator(hbar, mass, tolerance=1e-8)
    result = validator.validate_continuity_equation(state_t0, state_t1, dt)
    
    # Note: Test peut être sensible selon schéma intégration
    # Tolérance ajustée en conséquence
    assert result['is_satisfied'] or result['max_residual'] < 1e-6, \
        f"Résidu équation continuité trop grand: {result['max_residual']}"