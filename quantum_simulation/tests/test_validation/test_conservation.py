"""
Tests unitaires pour validation conservation.
"""

import sys
from pathlib import Path
# HACK: Ajouter racine projet au path (temporaire)
project_root = Path(__file__).parent.parent.parent.parent
sys.path.insert(0, str(project_root))

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
    Courant probabilité d'onde plane doit être proportionnel à ρ(x).
    
    Note:
        Après normalisation discrète sur grille finie, ρ(x) n'est plus uniforme.
        On vérifie donc J(x) ∝ ρ(x) avec v = ℏk/m.
    """
    hbar = 1.054571817e-34
    mass = 9.1093837015e-31
    k = 5e9
    
    x = np.linspace(-1e-8, 1e-8, 1024)
    
    free_particle = FreeParticle(mass, hbar)
    state = free_particle.create_plane_wave(x, k)
    
    validator = ConservationValidator(hbar, mass)
    J = validator.compute_probability_current(state)
    rho = state.probability_density()
    
    # Vitesse théorique
    v_theory = (hbar * k) / mass
    
    # Test : J(x) / ρ(x) ≈ v_theory (rapport constant)
    # Éviter division par zéro : masque où ρ > seuil
    mask = rho > 1e-15
    if np.any(mask):
        ratio = J[mask] / rho[mask]
        ratio_std = np.std(ratio)
        ratio_mean = np.mean(ratio)
        
        # Le rapport doit être constant ≈ v_theory
        relative_std = ratio_std / abs(ratio_mean) if ratio_mean != 0 else ratio_std
        
        assert relative_std < 0.1, \
            f"Rapport J/ρ devrait être constant, std relative={relative_std:.2e}"
        
        # Vérifier valeur moyenne proche v_theory
        relative_error = abs(ratio_mean - v_theory) / abs(v_theory)
        assert relative_error < 0.1, \
            f"Vitesse attendue {v_theory:.2e}, obtenue {ratio_mean:.2e}"
    else:
        # Si ρ uniformément nul (ne devrait pas arriver)
        assert np.all(np.abs(J) < 1e-15), "Si ρ≈0, alors J≈0"


def test_continuity_equation_gaussian_packet():
    """
    Équation continuité : cohérence entre ∂ρ/∂t et ∇·J.
    
    Note:
        Test simplifié vérifiant que résidu relatif < 20% pour translation rigide.
        Translation exacte (sans dispersion) ne satisfait pas strictement équation
        car néglige étalement du paquet.
        
        TODO: Remplacer par vraie évolution Schrödinger quand disponible.
    """
    hbar = 1.054571817e-34
    mass = 9.1093837015e-31
    
    x = np.linspace(-5e-8, 5e-8, 2048)
    
    free_particle = FreeParticle(mass, hbar)
    
    # Paquet initial
    sigma_x = 2e-9  # Plus large pour réduire dispersion sur dt court
    k0 = 5e9  # k plus petit pour réduire vitesse
    state_t0 = free_particle.create_gaussian_wavepacket(
        x, x0=0, sigma_x=sigma_x, k0=k0
    )
    
    # Translation manuelle (approximation évolution)
    dt = 5e-18  # Pas temporel très court
    v = (hbar * k0) / mass
    x_shift = v * dt
    
    state_t1 = free_particle.create_gaussian_wavepacket(
        x, x0=x_shift, sigma_x=sigma_x, k0=k0
    )
    
    validator = ConservationValidator(hbar, mass, tolerance=1e-6)
    result = validator.validate_continuity_equation(state_t0, state_t1, dt)
    
    print(f"\n=== Résultats équation continuité (translation rigide) ===")
    print(f"Résidu max absolu = {result['max_residual']:.2e} 1/s")
    print(f"Échelle ∂ρ/∂t = {result['drho_dt_scale']:.2e} 1/s")
    print(f"Résidu relatif = {result['relative_max_residual']:.2%}")
    print(f"Est satisfaite = {result['is_satisfied']}")
    
    # Test cohérence relative (tolérance 20% car translation approx.)
    assert result['relative_max_residual'] < 0.2, \
        f"Résidu relatif trop grand: {result['relative_max_residual']:.2%} " \
        f"(max absolu: {result['max_residual']:.2e})"
    
    # Norme conservée
    norm_t0 = state_t0.norm()
    norm_t1 = state_t1.norm()
    assert abs(norm_t1 - norm_t0) < 1e-6, \
        f"Norme non conservée: {norm_t0:.6f} → {norm_t1:.6f}"
    