"""
Tests unitaires pour module core.state.

Vérifie règles physiques fondamentales :
- R5.1 : Conservation normalisation ⟨ψ|ψ⟩ = 1
- Propriétés produit scalaire hermitien
"""

import sys
from pathlib import Path
# HACK: Ajouter racine projet au path (temporaire)
# print("PATH", Path(__file__).parent.parent.parent)
project_root = Path(__file__).parent.parent.parent
sys.path.insert(0, str(project_root))

import pytest
import numpy as np
from quantum_simulation.core.state import WaveFunctionState
from quantum_simulation.systems.free_particle import FreeParticle


def test_wavefunction_normalization():
    """
    Règle R5.1 : ⟨ψ|ψ⟩ = 1
    Source : [file:1, Chapitre III, § D-1-c]
    """
    hbar = 1.054571817e-34
    mass = 9.1093837015e-31
    
    x = np.linspace(-5e-8, 5e-8, 1024)
    
    free_particle = FreeParticle(mass, hbar)
    state = free_particle.create_gaussian_wavepacket(
        x, x0=0, sigma_x=2e-9, k0=5e9
    )
    
    # Vérification normalisation
    norm = state.norm()
    assert abs(norm - 1.0) < 1e-10, f"Norme devrait être 1, obtenu {norm}"
    assert state.is_normalized(tolerance=1e-10)


def test_inner_product_hermitian():
    """
    Propriété hermitienne : ⟨φ|ψ⟩* = ⟨ψ|φ⟩
    Source : [file:1, Chapitre II, § B-2-c]
    """
    hbar = 1.054571817e-34
    mass = 9.1093837015e-31
    
    x = np.linspace(-5e-8, 5e-8, 512)
    
    free_particle = FreeParticle(mass, hbar)
    psi = free_particle.create_gaussian_wavepacket(x, x0=-1e-8, sigma_x=1e-9, k0=3e9)
    phi = free_particle.create_gaussian_wavepacket(x, x0=1e-8, sigma_x=1.5e-9, k0=4e9)
    
    # Calcul produits scalaires
    phi_psi = phi.inner_product(psi)
    psi_phi = psi.inner_product(phi)
    
    # Vérification : ⟨φ|ψ⟩* = ⟨ψ|φ⟩
    assert abs(np.conj(phi_psi) - psi_phi) < 1e-12, \
        "Produit scalaire doit être hermitien"


def test_normalize_non_normalized_state():
    """
    Test méthode normalize() sur état non normé.
    """
    x = np.linspace(-1e-8, 1e-8, 256)
    dx = x[1] - x[0]
    
    # Créer fonction d'onde arbitraire non normée
    psi_unnormalized = np.exp(-((x)**2) / (2 * (1e-9)**2))
    
    state = WaveFunctionState(x, psi_unnormalized)
    
    # État initial non normé
    assert not state.is_normalized(tolerance=1e-10)
    
    # Normalisation
    state_normalized = state.normalize()
    
    # Vérifier normalisation
    assert state_normalized.is_normalized(tolerance=1e-10)
    
    # Vérifier conservation forme (seulement facteur multiplicatif)
    ratio = state_normalized.wavefunction / state.wavefunction
    assert np.allclose(ratio, ratio[0]), "Forme doit être préservée"


def test_probability_density_positive():
    """
    ρ(x) = |ψ(x)|² doit être positif partout.
    Règle R2.1
    """
    hbar = 1.054571817e-34
    mass = 9.1093837015e-31
    
    x = np.linspace(-3e-8, 3e-8, 512)
    
    free_particle = FreeParticle(mass, hbar)
    state = free_particle.create_gaussian_wavepacket(x, x0=0, sigma_x=1e-9, k0=0)
    
    rho = state.probability_density()
    
    # Positivité
    assert np.all(rho >= 0), "Densité probabilité doit être positive"
    
    # Intégrale = 1
    integral = np.trapz(rho, x)
    assert abs(integral - 1.0) < 1e-8, "Intégrale ρ doit être 1"


def test_probability_in_volume():
    """
    Test calcul probabilité dans intervalle.
    """
    hbar = 1.054571817e-34
    mass = 9.1093837015e-31
    
    x = np.linspace(-10e-8, 10e-8, 2048)
    
    free_particle = FreeParticle(mass, hbar)
    state = free_particle.create_gaussian_wavepacket(x, x0=0, sigma_x=2e-9, k0=0)
    
    # Probabilité dans [-σ, +σ] devrait être ~68% (gaussienne)
    sigma_x = 2e-9
    prob_1sigma = state.probability_in_volume(-sigma_x, sigma_x)
    
    # Tolérance large (discrétisation)
    assert 0.6 < prob_1sigma < 0.75, \
        f"Probabilité dans 1σ devrait être ~68%, obtenu {prob_1sigma:.2%}"
    
    # Probabilité totale
    prob_total = state.probability_in_volume(x[0], x[-1])
    assert abs(prob_total - 1.0) < 1e-6, "Probabilité totale doit être 1"


def test_inner_product_with_self():
    """
    ⟨ψ|ψ⟩ doit être réel et égal à norme²
    """
    hbar = 1.054571817e-34
    mass = 9.1093837015e-31
    
    x = np.linspace(-5e-8, 5e-8, 512)
    
    free_particle = FreeParticle(mass, hbar)
    state = free_particle.create_gaussian_wavepacket(x, x0=0, sigma_x=1.5e-9, k0=2e9)
    
    inner_self = state.inner_product(state)
    norm_squared = state.norm()**2
    
    # Doit être réel
    assert abs(inner_self.imag) < 1e-12, "⟨ψ|ψ⟩ doit être réel"
    
    # Égal à norme²
    assert abs(inner_self.real - norm_squared) < 1e-10


def test_incompatible_grids_raise_error():
    """
    Produit scalaire avec grilles différentes doit lever erreur.
    """
    x1 = np.linspace(-1e-8, 1e-8, 256)
    x2 = np.linspace(-2e-8, 2e-8, 512)  # Grille différente
    
    psi1 = np.exp(-((x1)**2) / (2 * (1e-9)**2))
    psi2 = np.exp(-((x2)**2) / (2 * (1e-9)**2))
    
    state1 = WaveFunctionState(x1, psi1).normalize()
    state2 = WaveFunctionState(x2, psi2).normalize()
    
    # Pattern regex corrigé pour matcher message d'erreur exact
    with pytest.raises(ValueError, match=r"Grilles spatiales incompatibles"):
        state1.inner_product(state2)