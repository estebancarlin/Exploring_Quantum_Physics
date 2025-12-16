"""
Système particule libre (V = 0).

Règle R3.2 avec V=0 : H = P²/2m = -ℏ²/2m Δ
Source : [file:1, Chapitre I, § C]
"""

import numpy as np
from typing import Tuple
from quantum_simulation.core.state import WaveFunctionState
from quantum_simulation.core.operators import Hamiltonian


class FreeParticle:
    """
    Particule quantique libre (pas de potentiel).
    
    États propres : ondes planes exp(ik·r)
    Spectre énergie : E(k) = ℏ²k²/2m (continu)
    """
    
    def __init__(self, mass: float, hbar: float):
        """
        Args:
            mass: Masse particule (kg)
            hbar: Constante de Planck réduite (J·s)
        """
        self.mass = mass
        self.hbar = hbar
        
        # Potentiel nul
        self.potential = lambda x: np.zeros_like(x)
        
        # Hamiltonien
        self.hamiltonian = Hamiltonian(
            mass=mass,
            potential=self.potential,
            hbar=hbar
        )
        
    def create_plane_wave(self, spatial_grid: np.ndarray, k: float) -> WaveFunctionState:
        """
        Crée onde plane exp(ikx) (état propre énergie).
        
        Args:
            spatial_grid: Grille spatiale
            k: Nombre d'onde (m⁻¹)
            
        Returns:
            État propre normalisé de H avec énergie E = ℏ²k²/2m
            
        Note: Sur grille finie, normalisation discrète appliquée.
        """
        psi = np.exp(1j * k * spatial_grid)
        
        # Normalisation sur grille
        dx = spatial_grid[1] - spatial_grid[0]
        norm = np.sqrt(np.sum(np.abs(psi)**2) * dx)
        psi /= norm
        
        return WaveFunctionState(spatial_grid, psi)
        
    def create_gaussian_wavepacket(self, spatial_grid: np.ndarray,
                                    x0: float, sigma_x: float, k0: float) -> WaveFunctionState:
        """
        Crée paquet d'ondes gaussien.
        
        ψ(x) = (2πσₓ²)^(-1/4) exp[-(x-x₀)²/(4σₓ²) + ik₀x]
        
        Args:
            spatial_grid: Grille spatiale
            x0: Position centrale (m)
            sigma_x: Largeur gaussienne (m)
            k0: Impulsion moyenne ℏk₀ (m⁻¹)
            
        Returns:
            État gaussien normalisé
            
        Propriétés:
            - ⟨X⟩ = x0
            - ΔX = σx
            - ⟨P⟩ = ℏk0
            - ΔP = ℏ/(2σx)
            - ΔX·ΔP = ℏ/2 (état minimum incertitude)
        """
        x = spatial_grid
        
        # Enveloppe gaussienne
        normalization = (2 * np.pi * sigma_x**2)**(-0.25)
        envelope = np.exp(-(x - x0)**2 / (4 * sigma_x**2))
        
        # Modulation onde plane
        phase = np.exp(1j * k0 * x)
        
        psi = normalization * envelope * phase
        
        # Vérification normalisation sur grille discrète
        state = WaveFunctionState(spatial_grid, psi)
        if not state.is_normalized(tolerance=1e-6):
            state = state.normalize()
            
        return state
        
    def energy_eigenvalue(self, k: float) -> float:
        """
        Calcule énergie onde plane de nombre d'onde k.
        
        E(k) = ℏ²k²/2m
        
        Returns:
            Énergie (J)
        """
        return (self.hbar**2 * k**2) / (2 * self.mass)
        
    def momentum_from_k(self, k: float) -> float:
        """
        Impulsion associée au nombre d'onde k.
        
        p = ℏk (Règle R1.2)
        
        Returns:
            Impulsion (kg·m/s)
        """
        return self.hbar * k
        
    def k_from_momentum(self, p: float) -> float:
        """
        Nombre d'onde associé à l'impulsion p.
        
        k = p/ℏ
        
        Returns:
            Nombre d'onde (m⁻¹)
        """
        return p / self.hbar
        
    def group_velocity(self, k: float) -> float:
        """
        Vitesse de groupe vg = dE/dℏk = ℏk/m.
        
        Pour particule libre, correspond à vitesse classique p/m.
        
        Returns:
            Vitesse (m/s)
        """
        return (self.hbar * k) / self.mass