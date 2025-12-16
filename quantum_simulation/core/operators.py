import numpy as np
from typing import Callable

from abc import ABC, abstractmethod

from quantum_simulation.core.state import WaveFunctionState, QuantumState
from quantum_simulation.utils.numerical import gradient_1d, laplacian_1d

class Observable(ABC):
    """
    Observable quantique (opérateur hermitique).
    Sources :
    - [file:1, Chapitre II, § D-1] : "observable = opérateur hermitique"
    - [file:1, Chapitre II, § D-2] : valeurs propres réelles
    
    Invariants :
    - Hermiticité : A† = A
    - Valeurs propres réelles
    - Vecteurs propres forment base orthonormée
    """
    
    @abstractmethod
    def apply(self, state: QuantumState) -> QuantumState:
        """Application A|ψ⟩"""
        
    @abstractmethod
    def expectation_value(self, state: QuantumState) -> float:
        """
        Calcule ⟨A⟩ = ⟨ψ|A|ψ⟩
        Source : [file:1, Chapitre III, § C-4] (Règle 1.4.1)
        """
        
    @abstractmethod
    def uncertainty(self, state: QuantumState) -> float:
        """
        Calcule ΔA = √(⟨A²⟩ - ⟨A⟩²)
        Source : [file:1, Chapitre III, § C-5] (Règle 1.4.2)
        """
        
    @abstractmethod
    def eigensystem(self) -> tuple[np.ndarray, list[QuantumState]]:
        """
        Retourne (valeurs_propres, vecteurs_propres)
        
        Limite : Méthode de diagonalisation non spécifiée dans synthèse.
        Implémentation numérique laissée libre (numpy.linalg, etc.)
        """
        
    def is_hermitian(self, tolerance: float = 1e-10) -> bool:
        """Vérifie A† = A"""
        
    def commutator(self, other: 'Observable') -> 'Observable':
        """
        Calcule [A,B] = AB - BA
        Source : [file:1, Chapitre III, § C-6-a] : compatibilité si [A,B]=0
        """

class PositionOperator(Observable):
    """
    Observable position R.
    Source : [file:1, Chapitre II, § E] : observable R(X,Y,Z)
    
    En représentation position : multiplication par r
    """
    
    def apply(self, state: WaveFunctionState) -> WaveFunctionState:
        """Applique R|ψ⟩ : multiplication r·ψ(r)"""

class MomentumOperator(Observable):
    """
    Règle R1.3 : P = -iℏ∇ en représentation position
    Source : [file:1, Chapitre II, § E-2]
    """
    
    def __init__(self, hbar: float, dimension: int = 1):
        self.hbar = hbar
        self.dimension = dimension
        
    def apply(self, state: WaveFunctionState) -> WaveFunctionState:
        """
        Applique P|ψ⟩ = -iℏ∇ψ
        Méthode : Différences finies ordre 2
        """
        grad_psi = gradient_1d(state.wavefunction, state.dx, order=2)
        return WaveFunctionState(
            state.spatial_grid,
            -1j * self.hbar * grad_psi
        )
        
    def expectation_value(self, state: WaveFunctionState) -> float:
        """
        ⟨P⟩ = ⟨ψ|P|ψ⟩
        Règle R4.1
        """
        P_psi = self.apply(state)
        return state.inner_product(P_psi).real
        
    def uncertainty(self, state: WaveFunctionState) -> float:
        """
        ΔP = √(⟨P²⟩ - ⟨P⟩²)
        Règle R4.2
        """
        P_psi = self.apply(state)
        P2_psi = self.apply(P_psi)
        
        exp_P = state.inner_product(P_psi).real
        exp_P2 = state.inner_product(P2_psi).real
        
        variance = exp_P2 - exp_P**2
        if variance < 0:
            variance = 0  # Erreurs numériques
        return np.sqrt(variance)

class Hamiltonian(Observable):
    """
    Règle R3.2 : H = -ℏ²/2m Δ + V(r)
    Source : [file:1, Chapitre III, § B-5-b]
    """
    
    def __init__(self, mass: float, potential: Callable[[np.ndarray], np.ndarray], hbar: float):
        self.mass = mass
        self.potential = potential
        self.hbar = hbar
        
    def apply(self, state: WaveFunctionState) -> WaveFunctionState:
        """
        Applique H|ψ⟩ = [-ℏ²/2m Δ + V(r)]ψ
        """
        # Terme cinétique
        laplacian_psi = laplacian_1d(state.wavefunction, state.dx)
        kinetic_term = -(self.hbar**2 / (2 * self.mass)) * laplacian_psi
        
        # Terme potentiel
        V_values = self.potential(state.spatial_grid)
        potential_term = V_values * state.wavefunction
        
        return WaveFunctionState(
            state.spatial_grid,
            kinetic_term + potential_term
        )