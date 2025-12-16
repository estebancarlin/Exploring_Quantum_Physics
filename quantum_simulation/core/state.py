from abc import ABC, abstractmethod
import numpy as np
from quantum_simulation.utils.numerical import integrate_1d

class QuantumState(ABC):
    """
    État quantique |ψ⟩ d'un système.
    Sources : 
    - [file:1, Chapitre III, § B-1] : "vecteur d'état"
    - [file:1, Chapitre II, § B-2-a] : notations de Dirac
    
    Invariants physiques :
    - Normalisation : ⟨ψ|ψ⟩ = 1 (si normalisé)
    - Linéarité de l'espace des états
    """
    
    @abstractmethod
    def norm(self) -> float:
        """Calcule ⟨ψ|ψ⟩"""
        
    @abstractmethod
    def normalize(self) -> 'QuantumState':
        """Retourne état normé"""
        
    @abstractmethod
    def inner_product(self, other: 'QuantumState') -> complex:
        """
        Produit scalaire ⟨φ|ψ⟩
        Source : [file:1, Chapitre II, § B-2-c]
        """
        
    def is_normalized(self, tolerance: float = 1e-10) -> bool:
        """Vérifie si ⟨ψ|ψ⟩ ≈ 1"""

class WaveFunctionState(QuantumState):
    """
    Implémentation Règle R2.1 : ρ(r) = |ψ(r)|²
    Source : [file:1, Chapitre I, § B-2]
    """
    
    def __init__(self, spatial_grid: np.ndarray, wavefunction: np.ndarray):
        self.spatial_grid = spatial_grid
        self.wavefunction = wavefunction
        self.dx = spatial_grid[1] - spatial_grid[0]  # Grille uniforme
        
    def norm(self) -> float:
        """
        Calcule √⟨ψ|ψ⟩ via intégration numérique.
        Implémente Règle R5.1 : ⟨ψ|ψ⟩ doit rester = 1
        """
        rho = np.abs(self.wavefunction)**2
        return np.sqrt(integrate_1d(rho, self.dx))
        
    def normalize(self) -> 'WaveFunctionState':
        """Retourne état normé : |ψ⟩ / √⟨ψ|ψ⟩"""
        norm_value = self.norm()
        if norm_value < 1e-12:
            raise ValueError("État nul, normalisation impossible")
        return WaveFunctionState(
            self.spatial_grid,
            self.wavefunction / norm_value
        )
        
    def inner_product(self, other: 'WaveFunctionState') -> complex:
        """
        Calcule ⟨φ|ψ⟩ = ∫ φ*(x)ψ(x) dx
        Source : [file:1, Chapitre II, § B-2-c]
        """
        if not np.allclose(self.spatial_grid, other.spatial_grid):
            raise ValueError("Grilles spatiales incompatibles")
        integrand = np.conj(other.wavefunction) * self.wavefunction
        return integrate_1d(integrand, self.dx)
        
    def probability_density(self) -> np.ndarray:
        """Règle R2.1 : ρ(x) = |ψ(x)|²"""
        return np.abs(self.wavefunction)**2

class EigenStateBasis(QuantumState):
    """
    État décomposé sur base d'états propres {|un⟩}.
    Source : [file:1, Chapitre III, § D-2-a] : décomposition ψ = Σ cn|un⟩
    
    Attributs :
    - eigenstates : list[QuantumState]  # Base orthonormée
    - coefficients : np.ndarray          # Coefficients cn complexes
    - eigenvalues : np.ndarray           # Valeurs propres an associées
    """
    
    def __init__(self, eigenstates, coefficients, eigenvalues):
        """
        Hypothèses :
        - Base orthonormée : ⟨ui|uj⟩ = δij (Règle depuis [file:1, Chapitre II, § C-2-a])
        """
        
    def validate_orthonormality(self, tolerance: float = 1e-8) -> bool:
        """Vérifie ⟨ui|uj⟩ = δij"""
