from typing import Tuple
import numpy as np

class FreeParticle2D:
    """
    Particule libre en 2D (V=0).
    
    Hamiltonien : H = -(ℏ²/2m)(∂²/∂x² + ∂²/∂y²)
    
    Références:
        - Extension dimensionnelle Chapitre I, § C
        - Laplacien 2D : Document référence, § E1
    """
    
    def __init__(self, mass: float, hbar: float):
        self.mass = mass
        self.hbar = hbar
        
        # Hamiltonien 2D
        from quantum_simulation.core.operators import Hamiltonian2D
        self.hamiltonian = Hamiltonian2D(mass, hbar, potential=None)
        
    def create_gaussian_2d(self, X: np.ndarray, Y: np.ndarray,
                        x0: float, y0: float, sigma: float,
                        kx: float, ky: float) -> WaveFunctionState2D:
        """
        Paquet gaussien 2D.
        
        ψ(x,y,0) = N exp[-(x-x₀)²/4σ² - (y-y₀)²/4σ²] exp[i(kₓx + kᵧy)]
        
        Args:
            X, Y: Grilles meshgrid (nx × ny)
            x0, y0: Centre paquet (m)
            sigma: Largeur (identique x et y)
            kx, ky: Vecteur d'onde k⃗ = (kₓ, kᵧ)
            
        Returns:
            WaveFunctionState2D normalisé
            
        Validation:
            - ∫∫ |ψ|² dxdy = 1 ± 10⁻⁹
            - ⟨Pₓ⟩ = ℏkₓ, ⟨Pᵧ⟩ = ℏkᵧ
        """
        # Implémentation normalisée
        
    def expected_trajectory_2d(self, t: float, x0: float, y0: float,
                            kx: float, ky: float) -> Tuple[float, float]:
        """
        Trajectoire classique attendue (Ehrenfest 2D).
        
        Returns:
            (x(t), y(t)) = (x₀ + vₓt, y₀ + vᵧt) avec v⃗ = ℏk⃗/m
        """