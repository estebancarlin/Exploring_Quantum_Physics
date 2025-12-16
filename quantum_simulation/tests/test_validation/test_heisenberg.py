"""
Tests unitaires pour validation relations Heisenberg.
"""

import pytest
import numpy as np
from quantum_simulation.validation.heisenberg_relations import HeisenbergValidator
from quantum_simulation.systems.free_particle import FreeParticle


def test_gaussian_packet_minimum_uncertainty():
    """
    Paquet gaussien doit satisfaire ΔX·ΔP = ℏ/2 (état minimum incertitude).
    """
    # Configuration
    hbar = 1.054571817e-34
    mass = 9.1093837015e-31
    
    # Grille
    x = np.linspace(-1e-8, 1e-8, 1024)
    
    # Système
    free_particle = FreeParticle(mass, hbar)
    
    # État gaussien (minimum incertitude par construction)
    sigma_x = 1e-9
    state = free_particle.create_gaussian_wavepacket(
        spatial_grid=x,
        x0=0.0,
        sigma_x=sigma_x,
        k0=1e10
    )
    
    # Validation Heisenberg
    validator = HeisenbergValidator(hbar)
    result = validator.validate_position_momentum(state)
    
    # Assertions
    assert result['is_valid'], "Inégalité Heisenberg violée"
    
    # Pour gaussienne : ΔX·ΔP devrait être très proche de ℏ/2
    excess = result['excess']
    assert abs(excess) < 0.01, f"Excès trop grand pour état minimum incertitude : {excess}"
    
    # Test qualité
    quality = validator.compute_minimum_uncertainty_state_quality(state)
    assert quality > 0.99, f"Qualité minimum incertitude insuffisante : {quality}"


def test_plane_wave_large_uncertainty():
    """
    Onde plane doit avoir ΔX → ∞ (délocalisée).
    """
    hbar = 1.054571817e-34
    mass = 9.1093837015e-31
    
    x = np.linspace(-1e-8, 1e-8, 2048)
    
    free_particle = FreeParticle(mass, hbar)
    state = free_particle.create_plane_wave(x, k=5e9)
    
    validator = HeisenbergValidator(hbar)
    result = validator.validate_position_momentum(state)
    
    # Onde plane : ΔX grand, ΔP petit → produit >> ℏ/2
    assert result['is_valid']
    assert result['product'] > 10 * result['heisenberg_bound'], \
        "Onde plane devrait avoir produit incertitudes >> ℏ/2"