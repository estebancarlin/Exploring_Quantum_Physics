import numpy as np
from typing import Callable, Tuple

def gradient_1d(psi: np.ndarray, dx: float, order: int = 2) -> np.ndarray:
    """
    Calcule ∇ψ par différences finies.
    
    Implémente méthode numérique pour Règle R1.3 (P = -iℏ∇)
    
    Args:
        psi: Fonction d'onde sur grille 1D
        dx: Pas spatial
        order: Ordre différences finies (2 ou 4)
    
    Returns:
        Gradient ∂ψ/∂x
        
    Conditions aux bords : Dirichlet (ψ=0)
    """
    # À implémenter

def laplacian_1d(psi: np.ndarray, dx: float) -> np.ndarray:
    """
    Calcule Δψ = ∂²ψ/∂x² par différences finies ordre 2.
    
    Utilisé dans Règle R3.2 : H = -ℏ²/2m Δ + V
    
    Returns:
        Laplacien sur grille
    """
    # À implémenter

def integrate_1d(f: np.ndarray, dx: float) -> float:
    """
    Intégration numérique ∫f(x)dx par méthode trapèzes.
    
    Utilisé pour normalisation (Règle R2.1) : ∫|ψ|²dx = 1
    """
    return np.trapz(f, dx=dx)