from typing import Tuple
import numpy as np

class HarmonicOscillator2D:
    """
    Oscillateur harmonique 2D (découplé).
    
    Hamiltonien : H = ℏω₁(a₁†a₁ + 1/2) + ℏω₂(a₂†a₂ + 1/2)
    
    États propres : |n₁,n₂⟩ avec E_{n₁,n₂} = ℏω₁(n₁+1/2) + ℏω₂(n₂+1/2)
    
    Cas particuliers:
        - Isotrope : ω₁ = ω₂ → symétrie circulaire
        - Anisotrope : ω₁ ≠ ω₂ → ellipses
    
    Références:
        - Chapitre V, § B (1D) → Extension 2D
        - Complément BVI : Oscillateur 2D (si fourni)
    """
    
    def __init__(self, mass: float, omega_x: float, omega_y: float,
                hbar: float, n_max: int = 20):
        """
        Args:
            omega_x, omega_y: Pulsations (rad/s)
            n_max: Troncature base {|n₁,n₂⟩, 0 ≤ n₁,n₂ ≤ n_max}
        """
        
    def energy_eigenvalue(self, nx: int, ny: int) -> float:
        """E_{n₁,n₂} = ℏω₁(n₁+1/2) + ℏω₂(n₂+1/2)"""
        
    def wavefunction_2d(self, nx: int, ny: int, X: np.ndarray, 
                        Y: np.ndarray) -> WaveFunctionState2D:
        """
        Fonction d'onde ψ_{n₁,n₂}(x,y) = ψ_{n₁}(x) · ψ_{n₂}(y).
        
        Utilise scipy.special.eval_hermite (comme 1D).
        """